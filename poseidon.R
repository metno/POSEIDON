#!/usr/bin/env Rscript
# + POSEIDON - Precipitation Spatial Data Quality Control
# mailto: cristianl@met.no
# # https://github.com/metno/POSEIDON
#
# command line:
#  >poseidon.R input_file output_file [options]
# to list available options:
#  >poseidon.R --help 
#-----------------------------------------------------------------------------
#  This file is free software: you may copy, redistribute and/or modify it  
#  under the terms of the GNU General Public License as published by the  
#  Free Software Foundation, either version 2 of the License, or (at your  
#  option) any later version.  
#  
#  This file is distributed in the hope that it will be useful, but  
#  WITHOUT ANY WARRANTY; without even the implied warranty of  
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  
#  General Public License for more details.  
#  
#  You should have received a copy of the GNU General Public License  
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#-----------------------------------------------------------------------------
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("sp"))
suppressPackageStartupMessages(library("raster"))
suppressPackageStartupMessages(library("rgdal"))
suppressPackageStartupMessages(library("tripack"))
options(warn = 2, scipen = 999)
#------------------------------------------------------------------------------
# FUNCTIONS

# auxiliary function to keep/blacklist observations
setCode_lonlat<-function(lonlat,code) {
# lonlat. vector. 1=lon; 2=lat
  ix<-datatmp$lon==lonlat[1] & datatmp$lat==lonlat[2]
  if (length(ix)>0)  {
    aux[ix]<-code
    assign("aux",aux,envir=.GlobalEnv)
  }
  return(length(ix))
}

#+ number of stations within "drmin" units from location "xy" 
nstat<-function(xy,drmin) {
# input
# xy=vector. 2=x;3=y
# output
# 1=nobs
#------------------------------------------------------------------------------
  if (any(is.na(xy))) return(NA)
  i<-which(abs(xy[1]-xtot)<=drmin & 
           abs(xy[2]-ytot)<=drmin)
  return(length(i))
}

#+ summary statistics in a box of /pm "drmin" centered on location "ixyz" 
statSpat<-function(ixyzt,drmin,tcor.flag=T,gamma=-0.0065) {
# input
# ixyz=vector. 1=index;2=x;3=y;4=z;5=t
# output
# 1=nobs;2=maxVertDist[m];3=mean(temp);4=sd(temp)
# NOTE: temperature is adjusted for elevation differences (gamma=-0.0065 K/m)
#------------------------------------------------------------------------------
  if (any(is.na(ixyzt))) return(c(NA,NA,NA,NA))
  i<-which(abs(ixyzt[2]-xtot)<=drmin & 
           abs(ixyzt[3]-ytot)<=drmin)
  if (length(i)==1) return(c(1,NA,NA,NA))
  dz<-ixyzt[4]-ztot[i]
  if (tcor.flag) {
    tcor<-ttot[i]+gamma*dz
  } else {
    tcor<-ttot[i]
  }
  tmean<-mean(tcor) 
  tsd<-sd(tcor)
  dz_mx<-max(abs(dz))
  return(c(length(i),round(dz_mx,0),round(tmean,1),round(tsd,3)))
}

#+ Box-Cox transformation
boxcox<-function(x,lambda) {
  if (lambda==0) {
    return(log(x))
  } else {
    return((x**lambda-1)/lambda)
  }
}

#+ Box-Cox inverse transformation
tboxcox<-function(x,lambda) {
  if (lambda==0) {
    return(exp(x))
  } else {
    return((1+lambda*x)**(1./lambda))
  }
}

#+ function to plot points
plotp<-function(x,y,val,br,col,
                map=NULL,map.br=NULL,map.col=NULL,
                xl=NULL,yl=NULL) {
  if (is.null(xl)) {
    plot(x,y)
  } else {
    plot(x,y,xlim=xl,ylim=yl)
  }
  if (!is.null(map)) {
    if (!is.null(map.br)) {
      image(map,add=T,breaks=map.br,col=map.col)
    } else {
      image(map,add=T)
    }
  } 
  for (i in 1:length(col)) {
    ix<-which(val>=br[i] & val<br[i+1])
    if (length(ix)>0) 
      points(x[ix],y[ix],col=col[i],pch=19)
    if (i==1) {
      legstr<-paste("<=",br[1],sep="")
    } else if (i==length(col)) {
      legstr<-c(legstr,paste(">=",br[length(col)],sep=""))
    } else {
      legstr<-c(legstr,br[i])
    }
  }
  legend(x="bottomright",fill=rev(col),legend=rev(legstr))
}

#+ SCT - spatial consistency test
sct<-function(ixynp,
              nmin=50,dzmin=30,
              Dhmin=10000,Dz=200,eps2=0.5,
              T2=16,sus.code=4) {
# ref:
#  Lussana, C., Uboldi, F., & Salvati, M. R. (2010). A spatial consistency 
#   test for surface observations from mesoscale meteorological networks.
#   Quarterly Journal of the Royal Meteorological Society, 136(649), 1075-1088.
# input
#  ixynp= vector(4). 1=box identifier on the grid; 
#                   2/3=easting/northing coord (center of the box);
#                   4=number of stations within the box
#  NOTE: stations are associated to a box before running this function
#  nmin= numeric. minimum number of stations to fit a vertical profile
#  dzmin= numeric. minimum elevation range to fit a vertical profile [m]
#  Dhmin= numeric. minimum value for OI horizontal decorellation length [m]
#  Dz= numeric. OI vertical decorellation length [m]
#  eps2= numeric. OI ratio between obs_err_variance/backg_err_variance
#  T2=numeric. SCT threshold. (obs-pred)^2/(varObs+varPred)^2 > T2, suspect!
#  sus.code=numeric. identifier code for suspect observation
# output
#  number of rejected stations. (special cases: (i) NA if the function has not
#   been applied; (ii) -1 if just one station in the domain
#  
#  NOTE: "dqcflag" global variable is also updated
#------------------------------------------------------------------------------
  # something strange with the number of stations
  if (is.na(ixynp[4]) | is.null(ixynp[4]) | !is.finite(ixynp[4]) ) return(NA)
  # j, index for the stations in the box
  j<-which(itot==ixynp[1])
  # case of just one station
  if (ixynp[4]==1) {
    sctpog[ix[j]]<--1
    assign("sctpog",sctpog,envir=.GlobalEnv)
    return(-1)
  }
  # "zopt" and "topt" are set as global variables, so that they can be used by 
  #   other functions in the optimizaiton procedure
#  assign("zopt",ztot[j],envir=.GlobalEnv)
#  dz<-as.numeric(quantile(zopt,probs=0.95))-
#      as.numeric(quantile(zopt,probs=0.05))
#  assign("topt",ttot[j],envir=.GlobalEnv)
  to<-ttot[j]
  ta<-to
  tav<-to
  ta[]<-NA
  tav[]<-NA
  #
  tb<-tbtot[j]
  # OI for SCT (Lussana et al., 2010)
  # initialize variables
  #  dqc flags
  dqctmp<-dqcflag[ix[j]]
  #  probability of gross error (pog)
  pog<-dqctmp
  pog[]<-NA
  # distance matrices
  disth<-(outer(xtot[j],xtot[j],FUN="-")**2.+
          outer(ytot[j],ytot[j],FUN="-")**2.)**0.5
  distz<-abs(outer(ztot[j],ztot[j],FUN="-"))
  # set to optimal Dh
  Dh<-max(Dhmin,
          as.numeric(quantile(disth,probs=0.1)))
  # background error correlation matrix
  S<-exp(-0.5*(disth/Dh)**2.-0.5*(distz/Dz)**2.)
  if (argv$laf.sct) {
    S<-S * (1-(1-argv$lafmin.sct)*abs(outer(laftot[j],laftot[j],FUN="-")))
  }
  # S+eps2I
  diag(S)<-diag(S)+eps2
  # innvoation
  d<-to-tb
  # select observations to test 
  sel<-which(is.na(dqctmp) | dqctmp==keep.code)
  sel2check<-which(is.na(dqctmp))
  first<-T
  # loop over SCT iterations 
  # NOTE: SCT flags at most one observation, iterate until no observations fail
  while (length(sel)>1) { 
#    # update selection
#    sel<-which(is.na(dqctmp) | dqctmp==keep.code)
#    sel2check<-which(is.na(dqctmp))
#    if (length(sel2check)==0) break
    # first iteration, inver the matrix
    if (first) {
      SRinv<-chol2inv(chol(S))
      # from S+R go back to S
      diag(S)<-diag(S)-eps2
      first<-F
    } else {
      # Update inverse matrix (Uboldi et al 2008, Appendix AND erratum!)
      aux<-SRinv
      SRinv<-aux[-indx,-indx]-
             (tcrossprod(aux[indx,-indx],aux[-indx,indx]))*Zinv[indx]
      S<-S[-indx,-indx]
      rm(aux)
    }
    # next tree lines: compute cvres=(obs - CrossValidation_prediction)
    Zinv<-1/diag(SRinv)
    SRinv.d<-crossprod(SRinv,d[sel])
    ares<-crossprod(S,SRinv.d)-d[sel] #   a-Obs
    cvres<--Zinv*SRinv.d              # CVa-Obs
    ta[sel]<-to[sel]+ares
    tav[sel]<-to[sel]+cvres
    sig2o<-mean(d[sel]*(-ares))       # Lussana et al 2010, Eq(32)
    if (sig2o<0.01) sig2o<-0.01       # safe threshold  
    # pog=cvres/(sig2obs+sig2CVpred), Lussana et al 2010 Eq(20)
    pog[sel]<-(ares*cvres)/sig2o
    sctpog[ix[j[sel]]]<-pog[sel]
    assign("sctpog",sctpog,envir=.GlobalEnv)
    if (length(sel2check)==0) break
    # check if any obs fails the test
    if (any(pog[sel2check]>T2)) {
      # flag as suspect only the observation with the largest cvres 
      indx<-which.max(pog[sel2check])
      dqctmp[sel2check[indx]]<-sus.code
      # update global variable with flags
      dqcflag[ix[j[sel2check[indx]]]]<-sus.code
      assign("dqcflag",dqcflag,envir=.GlobalEnv)
      # update corep (useful from 2nd iteration onwards, if an obs fails the
      # sct, then no need for representativeness
      corep[ix[j[sel2check[indx]]]]<-NA
      assign("corep",corep,envir=.GlobalEnv)
      # update selection
      sel<-which(is.na(dqctmp) | dqctmp==keep.code)
      sel2check<-which(is.na(dqctmp))
    } else {
      break
    }
  } # end cycle SCT model
  # coefficient of observation representativeness
  # this call to ecdf(x)(x) should be the same as rank(x)/length(x)
  corep[ix[j[sel]]]<-(d[sel]*(-ares))/sig2o
  assign("corep",corep,envir=.GlobalEnv)
  # debug: begin
  if (argv$debug) {
    susi<-which(!is.na(dqctmp) & dqctmp!=keep.code)
    gold<-which(dqctmp==keep.code)
    if (!dir.exists(argv$debug.dir)) 
      dir.create(argv$debug.dir,showWarnings=F,recursive=T)
    f<-file.path(argv$debug.dir,
         paste("poseidon_sctm_it_",formatC(i,width=2,flag="0"),
               "_subd",formatC(ixynp[1],width=5,flag="0"),".png",sep=""))
    png(file=f,width=800,height=800)
    aux<-which(is.na(dqctmp))
    plot(to[aux],tb[aux],
         xlim=c(-2,max(c(to,tb,ta,tav,na.rm=T))),
         ylim=c(-2,max(c(to,tb,ta,tav,na.rm=T))),
         pch=19,col="black",cex=2)
    points(to[aux],tav[aux],pch=19,col="blue",cex=2)
    points(to[aux],ta[aux],pch=19,col="cyan",cex=2)
    lines(-100:100,-100:100,col="gray",lty=1)
    abline(h=seq(-100,100,by=0.25),lty=2,col="gray")
    abline(v=seq(-100,100,by=0.25),lty=2,col="gray")
    aux<-which(dqctmp==sct.code)
    points(to[aux],tb[aux],pch=17,col="black",cex=2.2)
    points(to[aux],tav[aux],pch=17,col="red",cex=2.2)
    points(to[aux],ta[aux],pch=17,col="pink",cex=2.2)
    dev.off()
    f<-file.path(argv$debug.dir,
         paste("poseidon_scthorz_it_",formatC(i,width=2,flag="0"),
               "_subd",formatC(ixynp[1],width=5,flag="0"),".png",sep=""))
    png(file=f,width=800,height=800)
    plot(xtot[j],ytot[j])
    xmna<-ixynp[2]-100000
    xmxa<-ixynp[2]+100000
    ymna<-ixynp[3]-100000
    ymxa<-ixynp[3]+100000
    if ( (argv$dem | argv$dem.fill) &
        (xmna>=rdem.xmn & xmxa<=rdem.xmx & 
         ymna>=rdem.ymn & ymxa<=rdem.ymx) ) {
      if (!exists("dem1"))
        dem1<-crop(rdem,extent(c(xmna,xmxa,ymna,ymxa)))
      image(dem1,add=T,
            breaks=c(0,10,25,50,100,250,500,750,1000,1250,1500,1750,2000,2500,3000),
            col=gray.colors(14))
    }
    brt<-c(0,0.1,0.5,1,2,3,4,5,6,7,8,9,10,3000)
    br<-boxcox(brt,
               lambda=argv$boxcox.lambda)
    col<-c("beige",rev(rainbow((length(br)-2))))
    for (c in 1:length(col)) {
      if (c==1) {
        legstr<-paste("<",brt[2],"mm",sep="")
      } else if (c==length(col)) {
        legstr<-c(legstr,paste(">=",brt[c],sep=""))
      } else {
        legstr<-c(legstr,paste("[",brt[c],",",brt[c+1],")",sep=""))
      }
      aux<-which(to>=br[c] & 
                 to<br[c+1] & 
                 is.na(dqctmp))
      if (length(aux)>0) {
        points(xtot[j][aux],ytot[j][aux],pch=19,col=col[c],cex=2)
        points(xtot[j][aux],ytot[j][aux],pch=1,col="black",cex=2)
      }
      aux<-which(to>=br[c] & 
                 to<br[c+1] & 
                 dqctmp==sct.code)
      if (length(aux)>0) {
        points(xtot[j][aux],ytot[j][aux],pch=17,col=col[c],cex=2)
        points(xtot[j][aux],ytot[j][aux],pch=2,col="black",cex=2)
      }
    }
    legend(x="bottomright",fill=rev(col),legend=rev(legstr),cex=1.5)
    dev.off()
  }
  # debug: end
  return(length(which(dqctmp==sus.code)))
}

# used for debug
plotSCTgrid<-function() {
  png(file="domain.png",width=800,height=800)
  par(mar=c(1,1,1,1))
  image(rlaf,breaks=c(0,0.000001,1),col=c("cornflowerblue","beige"),axes=F,
        main="",xlab="",ylab="")
  plot(e,add=T)
  #points(xr,yr,pch=19,col="black")
  xrv<-unique(xr)
  yrv<-unique(yr)
  yrmn<-min(yr,na.rm=T)-res(r)[2]/2
  yrmx<-max(yr,na.rm=T)+res(r)[2]/2
  xrmn<-min(xr,na.rm=T)-res(r)[1]/2
  xrmx<-max(xr,na.rm=T)+res(r)[1]/2
  for (iii in 1:(length(yrv))) {
    lines(c(xrmn,xrmx),c(yrv[iii]+res(r)[2]/2,yrv[iii]+res(r)[2]/2),lwd=2.5)
  }
  for (iii in 1:(length(xrv))) {
    lines(c(xrv[iii]+res(r)[1]/2,xrv[iii]+res(r)[1]/2),c(yrmn,yrmx),lwd=2.5)
  }
  box()
  dev.off()
}

# +
RR1_dqc_isolatedDry<-function(obs,
                              rrinf=1,
                              pmax=20,
                              dmax=150000,
                              n.sector=16,
                              dmin.dry=150000) {
#------------------------------------------------------------------------------
# [] Contiguous NO-Rain (nor) Area identification
  p<-length(obs$x)
  pmax<-min(p,pmax)
  psub<-vector(mode="numeric",length=p)
  psub[]<-NA
  sector.angle<-360/n.sector
  ydqc.flag<-vector(mode="numeric",length=p)
  ydqc.flag[]<-0
# First Guess
  dry_ix<-which(obs$yo<rrinf)
  if (length(dry_ix)==0) return(ydqc.flag)
  # vectors:
  # +nnor.vec>number of no-rain areas
  # +Lnor.vec[i]>(i=1,..,nnor.vec) number of stns in i-th no-rain area 
  # +nor.vec[i,n]-> n-th stn (i.e. pointers to VecS) in i-th no-rain area
  #  matrix(i,j) i=1,nnor.vec j=1,Lnor.vec[i]
  #  note: an no-rain areas could contain both wet and dry stations
  nor.vec<-matrix(data=NA,ncol=p,nrow=p)
  Lnor.vec<-vector(mode="numeric",length=p)
  Lnor.vec[]<-NA
  nnor.vec<-0
  nor.aux<-matrix(data=NA,ncol=p,nrow=p)
  Lnor.aux<-vector(mode="numeric",length=p)
  Lnor.aux[]<-NA
  nnor.aux<-0
  for (b in dry_ix) {  # START: Cycle over dry stations 
    # a. identify the closest stations to the b-th dry station
    #  outputs: -close2b-> pointer to VecS of the closest stations.
    #            constraints: distance must be less than Lsubsample.DHmax 
    #                         max number of stations allowed is Lsubsample.max
    #           -psub[b]-> actual number of stations used
    Disth<-sqrt((obs$y[b]-obs$y)**2.+(obs$x[b]-obs$x)**2.)
    if (any(Disth<=dmax)) {
      ix<-which(Disth<=dmax)
      psub[b]<-min(pmax,length(ix))
      c2b<-order(Disth,decreasing=F)[1:psub[b]]
    } else {
      psub[b]<-0
    }
    if (psub[b]<3) next
    close2b<-c2b[1:psub[b]]
    rm(c2b)
    # c. Establish (local-) triangulation (Delauney)
    #    note: all stations are used (wet and dry)
    aux<-cbind(obs$x[close2b],obs$y[close2b])
    nokloop<-1
    while (anyDuplicated(aux)) {
      if (nokloop%%2==0) {
        obs$x[close2b][which(duplicated(aux))]<-obs$x[close2b][which(duplicated(aux))]+1
      } else {
        obs$y[close2b][which(duplicated(aux))]<-obs$y[close2b][which(duplicated(aux))]+1
      }
      aux<-cbind(obs$x[close2b],obs$y[close2b])
      nokloop<-nokloop+1
      if (nokloop>10) {
        print("ERROR: check the input data for too many duplicated stations")
        quit(status=1)
      }
    }
    tri.rr<-tri.mesh(obs$x[close2b],obs$y[close2b])
    # d. identify all the stations (wheter wet or dry) which belongs
    #    to adjacent nodes (respect to the b-th dry station node)
    #  procedure: tri.find-> returns the three vertex indexes of triangle
    #   containing the specified point. To find all the neighbouring 
    #   triangles just span the area surrounding the b-th dry station
    #   (circle of radius=1Km)
    #  output: nodes-> station position respect to VecS
    nodes<-vector(mode="numeric")
    tri.loc<-tri.find(tri.rr,obs$x[b],obs$y[b])
    # note: due to the fact that (obs$x[b],obs$y[b]) belongs to the mesh
    # used to establish the triangulation, ...i1,i2,i3 must be >0
    nodes<-close2b[c(tri.loc$i1,tri.loc$i2,tri.loc$i3)]
    for (cont.clock in 1:n.sector) {
      x.aux<-obs$x[b]+sin(sector.angle*(cont.clock-1)*pi/180.)
      y.aux<-obs$y[b]+cos(sector.angle*(cont.clock-1)*pi/180.)
      tri.loc<-tri.find(tri.rr,x.aux,y.aux)
      nodes.aux<-close2b[c(tri.loc$i1,tri.loc$i2,tri.loc$i3)]
      nodes<-c(nodes,
               nodes.aux[which( (!(nodes.aux %in% nodes)) &
                                  (nodes.aux>0) )])
    }
    rm(x.aux,y.aux,tri.loc,nodes.aux)
    lnodes<-length(nodes)
    # e. update the temporary structure used to identify nor 
    #    merge in a (new) nor all the (old) temporary nor (if any) 
    #    which contain at least one station present in nodes (avoiding
    #    repetition); erase (i.e. set to NAs) the (old) temporary nor merged
    aux<-vector()
    if (nnor.aux>0) {
      for (nor in 1:nnor.aux) {
        if (Lnor.aux[nor]==0) next
        if ( any(nodes %in% nor.aux[nor,1:Lnor.aux[nor]]) ) {
          aux<-c(aux[which(!(aux %in% nor.aux[nor,1:Lnor.aux[nor]]))],nor.aux[nor,1:Lnor.aux[nor]])
          nor.aux[nor,]<-NA
          Lnor.aux[nor]<-0
        }
      }
    }
    nnor.aux<-nnor.aux+1
    Lnor.aux[nnor.aux]<-length(c(aux,nodes[which(!(nodes%in%aux))]))
    nor.aux[nnor.aux,1:Lnor.aux[nnor.aux]]<-c(aux,nodes[which(!(nodes%in%aux))])
  } # END: Cycle over the dry stations
# Reorganise no-rain area labels
  y.nor<-vector(length=p)
  y.nor[1:p]<-NA
  nnor.vec<-0
  if (nnor.aux>0) {
    for (nor in 1:nnor.aux) {
      if (Lnor.aux[nor]==0) next
      nnor.vec<-nnor.vec+1
      Lnor.vec[nnor.vec]<-Lnor.aux[nor]
      nor.vec[nnor.vec,1:Lnor.vec[nnor.vec]]<-nor.aux[nor,1:Lnor.aux[nor]]
      y.nor[nor.vec[nnor.vec,1:Lnor.vec[nnor.vec]]]<-nnor.vec
    }
  }
  rm(nor.aux,Lnor.aux,nnor.aux)
# DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC 
# identify dry-stations surrounded only by wet-stations (close enough)
  if (nnor.vec>0) {
    for (j in 1:nnor.vec) {
      nor.j<-nor.vec[j,1:Lnor.vec[j]]
      indx.dry.stns<-which(nor.j %in% dry_ix)
      n.dry<-length(indx.dry.stns)
      if (n.dry==1 | (n.dry<4 & (n.dry/Lnor.vec[j])<0.2) ) {
        nor.j.dry<-nor.j[indx.dry.stns]
        aux.dist<-(outer(obs$y[nor.j.dry],obs$y[nor.j],FUN="-")**2.+
                   outer(obs$x[nor.j.dry],obs$x[nor.j],FUN="-")**2.)**0.5
        min.dist.from.dry.stn<-min(aux.dist[aux.dist>0])
        if (min.dist.from.dry.stn<dmin.dry) 
          ydqc.flag[nor.j.dry]<-200
      }
    }
    rm(nor.j,indx.dry.stns)
  }
  if (exists("aux.dist")) rm(aux.dist,min.dist.from.dry.stn,nor.j.dry)
  ydqc.flag
}

RR1_dqc_isolatedWet<-function(obs,
                              rrinf=1,
                              pmax=20,
                              dmax=150000,
                              n.sector=16,
                              dmin.wet=100000) {
#------------------------------------------------------------------------------
# [] Precipitation events identification (eve) (contiguous rain areas)
  p<-length(obs$x)
  pmax<-min(p,pmax)
  psub<-vector(mode="numeric",length=p)
  psub[]<-NA
  sector.angle<-360/n.sector
  ydqc.flag<-vector(mode="numeric",length=p)
  ydqc.flag[]<-0
# First Guess
  wet_ix<-which(obs$yo>=rrinf)
  if (length(wet_ix)==0) return(ydqc.flag) 
  # vectors:
  # +neve.vec>number of events 
  # +Leve.vec[i]>(i=1,..,neve.vec) number of stns in i-th event
  # +eve.vec[i,n]-> n-th stn (i.e. pointers to VecS) in i-th event
  #  matrix(i,j) i=1,neve.vec j=1,Leve.vec[i]
  #  note: an event could contain both wet and dry stations
  eve.vec<-matrix(data=NA,ncol=p,nrow=p)
  Leve.vec<-vector(mode="numeric",length=p)
  Leve.vec[]<-NA
  neve.vec<-0
  eve.aux<-matrix(data=NA,ncol=p,nrow=p)
  Leve.aux<-vector(mode="numeric",length=p)
  Leve.aux[]<-NA
  neve.aux<-0
  for (b in wet_ix) {  # START: Cycle over wet stations 
    # a. identify the closest stations to the b-th wet station
    #  outputs: -close2b-> position (i.e. pointer to VecS) of the closest stations.
    #            constraints: distance must be less than Lsubsample.DHmax 
    #                         max number of stations allowed is Lsubsample.max
    #           -Lsubsample.vec[b]-> actual number of stations used
    Disth<-sqrt((obs$y[b]-obs$y)**2.+(obs$x[b]-obs$x)**2.)
    if (any(Disth<=dmax)) {
      ix<-which(Disth<=dmax)
      psub[b]<-min(pmax,length(ix))
      c2b<-order(Disth,decreasing=F)[1:psub[b]]
    } else {
      psub[b]<-0
    }
    if (psub[b]<3) next 
    close2b<-c2b[1:psub[b]]
    rm(c2b)
    # c. Establish (local-) triangulation (Delauney)
    #    note: all stations are used (wet and dry)
    aux<-cbind(obs$x[close2b],obs$y[close2b])
    nokloop<-1
    while (anyDuplicated(aux)) {
      if (nokloop%%2==0) {
        obs$x[close2b][which(duplicated(aux))]<-obs$x[close2b][which(duplicated(aux))]+1
      } else {
        obs$y[close2b][which(duplicated(aux))]<-obs$y[close2b][which(duplicated(aux))]+1
      }
      aux<-cbind(obs$x[close2b],obs$y[close2b])
      nokloop<-nokloop+1
      if (nokloop>10) {
        print("ERROR: check the input data for too many duplicated stations")
        quit(status=1)
      }
    }
    tri.rr<-tri.mesh(obs$x[close2b],obs$y[close2b])
    # d. identify all the stations (wheter wet or dry) which belongs
    #    to adjacent nodes (respect to the b-th wet station node)
    #  procedure: tri.find-> returns the three vertex indexes of triangle
    #   containing the specified point. To find all the neighbouring 
    #   triangles just span the area surrounding the b-th wet station
    #   (circle of radius=1Km)
    #  output: nodes-> station position respect to VecS
    nodes<-vector(mode="numeric")
    tri.loc<-tri.find(tri.rr,obs$x[b],obs$y[b])
    # note: due to the fact that (obs$x[b],obs$y[b]) belongs to the mesh
    # used to establish the triangulation, ...i1,i2,i3 must be >0
    nodes<-close2b[c(tri.loc$i1,tri.loc$i2,tri.loc$i3)]
    for (cont.clock in 1:n.sector) {
      x.aux<-obs$x[b]+sin(sector.angle*(cont.clock-1)*pi/180.)
      y.aux<-obs$y[b]+cos(sector.angle*(cont.clock-1)*pi/180.)
      tri.loc<-tri.find(tri.rr,x.aux,y.aux)
      nodes.aux<-close2b[c(tri.loc$i1,tri.loc$i2,tri.loc$i3)]
      nodes<-c(nodes,nodes.aux[which( (!(nodes.aux %in% nodes)) & (nodes.aux>0) )])
    }
    rm(x.aux,y.aux,tri.loc,nodes.aux)
    lnodes<-length(nodes)
    # e. update the temporary structure used to identify eve 
    #    merge in a (new) eve all the (old) temporary eve (if any) 
    #    which contain at least one station present in nodes (avoiding
    #    repetition); erase (i.e. set to NAs) the (old) temporary eve merged
    aux<-vector()
    if (neve.aux>0) {
      for (eve in 1:neve.aux) {
        if (Leve.aux[eve]==0) next
        if ( any(nodes %in% eve.aux[eve,1:Leve.aux[eve]]) ) {
          aux<-c(aux[which(!(aux %in% eve.aux[eve,1:Leve.aux[eve]]))],eve.aux[eve,1:Leve.aux[eve]])
          eve.aux[eve,]<-NA
          Leve.aux[eve]<-0
        }
      }
    }
    neve.aux<-neve.aux+1
    Leve.aux[neve.aux]<-length(c(aux,nodes[which(!(nodes%in%aux))]))
    eve.aux[neve.aux,1:Leve.aux[neve.aux]]<-c(aux,nodes[which(!(nodes%in%aux))])
  } # END: Cycle over the wet stations
# reorganise event labels
  y.eve<-vector(length=p)
  y.eve[1:p]<-NA
  neve.vec<-0
  if (neve.aux>0) {
    for (eve in 1:neve.aux) {
      if (Leve.aux[eve]==0) next
      neve.vec<-neve.vec+1
      Leve.vec[neve.vec]<-Leve.aux[eve]
      eve.vec[neve.vec,1:Leve.vec[neve.vec]]<-eve.aux[eve,1:Leve.aux[eve]]
      y.eve[eve.vec[neve.vec,1:Leve.vec[neve.vec]]]<-neve.vec
    }
  }
  rm(eve.aux,Leve.aux,neve.aux)
# DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC DQC 
# wet-stations surrounded only by dry-stations (close enough)
  if (neve.vec>0) {
    for (j in 1:neve.vec) {
      eve.j<-eve.vec[j,1:Leve.vec[j]]
      indx.wet.stns<-which(eve.j %in% wet_ix)
      n.wet<-length(indx.wet.stns)
      if (n.wet==1) {
        eve.j.wet<-eve.j[indx.wet.stns]
        aux.dist<-(outer(obs$y[eve.j.wet],obs$y[eve.j],FUN="-")**2.+
                   outer(obs$x[eve.j.wet],obs$x[eve.j],FUN="-")**2.)**0.5
        min.dist.from.wet.stn<-min(aux.dist[aux.dist>0])
        if (min.dist.from.wet.stn<dmin.wet) 
          ydqc.flag[eve.j.wet]<-300
      }
    }
  }
#  rm(eve.j,indx.wet.stns)
  if (exists("aux.dist")) rm(aux.dist,min.dist.from.wet.stn,eve.j.wet)
  ydqc.flag
}

#+ plot summary figures 
plotsummary<-function(ixynp) {
#------------------------------------------------------------------------------
  # something strange with the number of stations
  if (is.na(ixynp[4]) | is.null(ixynp[4]) | !is.finite(ixynp[4]) ) return(NA)
  xmna<-ixynp[2]-res(r)[1]/2
  xmxa<-ixynp[2]+res(r)[1]/2
  ymna<-ixynp[3]-res(r)[2]/2
  ymxa<-ixynp[3]+res(r)[2]/2
  # j, index for the stations in the box
  j<-which(itot==ixynp[1])
  to<-data$value[j]
  dqctmp<-dqcflag[j]
  susj<-which(dqctmp!=0)
  good<-which(dqctmp==0)
  if (!dir.exists(argv$debug.dir)) 
    dir.create(argv$debug.dir,showWarnings=F,recursive=T)
  f<-file.path(argv$debug.dir,
       paste("poseidon_summary_subd",
             formatC(ixynp[1],width=5,flag="0"),
             ".png",sep=""))
  png(file=f,width=800,height=800)
  plot(x[j],y[j],col="white",xlim=c(xmna,xmxa),ylim=c(ymna,ymxa))
  if ( (argv$dem | argv$dem.fill) &
       (xmna>=rdem.xmn & xmxa<=rdem.xmx & 
        ymna>=rdem.ymn & ymxa<=rdem.ymx) ) {
    if (!exists("dem1"))
      dem1<-crop(rdem,extent(c(xmna,xmxa,ymna,ymxa)))
    image(dem1,add=T,
          breaks=c(0,10,25,50,100,250,500,750,
                   1000,1250,1500,1750,2000,2500,3000),
            col=gray.colors(14))
  }
  br<-c(0,0.1,0.5,1,2,3,4,5,6,7,8,9,10,3000)
  col<-c("beige",rev(rainbow((length(br)-2))))
  for (c in 1:length(col)) {
    if (c==1) {
      legstr<-paste("<",br[2],"mm",sep="")
    } else if (c==length(col)) {
      legstr<-c(legstr,paste(">=",br[c],sep=""))
    } else {
      legstr<-c(legstr,paste("[",br[c],",",br[c+1],")",sep=""))
    }
    aux<-which(to>=br[c] & 
               to<br[c+1] & 
               dqctmp==0)
    if (length(aux)>0) {
      points(x[j][aux],y[j][aux],pch=19,col=col[c],cex=2)
      points(x[j][aux],y[j][aux],pch=1,col="black",cex=2)
    }
    aux<-which(to>=br[c] & 
               to<br[c+1] & 
               dqctmp!=0)
    if (length(aux)>0) {
      points(x[j][aux],y[j][aux],pch=17,col=col[c],cex=2)
      points(x[j][aux],y[j][aux],pch=2,col="black",cex=2)
    }
  }
  legend(x="bottomright",fill=rev(col),legend=rev(legstr),cex=1.5)
  dev.off()
  # debug: end
  return(0)
}


#==============================================================================
#  MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN
#==============================================================================
t0<-Sys.time()
#
# create parser object
p <- arg_parser("poseidon")
# specify our desired options 
# by default ArgumentParser will add an help option 
p <- add_argument(p, "input",help="input file",type="character")
p <- add_argument(p, "output",help="output file",type="character",
                  default="output.txt")
# more than one provider
p <- add_argument(p, "--input.files",help="additional input files (provider2 provider3 ...)",
                  type="character",default=NULL,nargs=Inf,short="-i")
#
p <- add_argument(p, "--debug",help="debug mode",flag=T,short="-dbg")
p <- add_argument(p, "--debug.dir",help="directory for debug output",
                  type="character",default=".",short="-dbgd")
p <- add_argument(p, "--verbose",help="debug mode",flag=T,short="-v")
# NOTE: lat-lon setup to have Oslo in a single box
p <- add_argument(p, "--lonmin",help="longitude of south-eastern domain corner",
                  type="numeric",default=5,short="-lon")
p <- add_argument(p, "--lonmax",help="longitude of south-western domain corner",
                  type="numeric",default=28,short="-lox")
p <- add_argument(p, "--latmin",help="latitude of south-eastern domain corner",
                  type="numeric",default=53.25,short="-lan")
p <- add_argument(p, "--latmax",help="latitude of north-western domain corner",
                  type="numeric",default=71.8,short="-lax")
# variable names
p <- add_argument(p, "--separator",
                  help="separator character",
                  type="character",default=";")
p <- add_argument(p, "--varname.lat",
                  help="name for the latitude variable (in/out)",
                  type="character",default="lat",short="-vlat")
p <- add_argument(p, "--varname.lon",
                  help="name for the longitude variable (in/out)",
                  type="character",default="lon",short="-vlon")
p <- add_argument(p, "--varname.elev",
                  help="name for the elevation variable (in/out)",
                  type="character",default="elev",short="-vele")
p <- add_argument(p, "--varname.value",
                  help="name for the temperature variable (in/out)",
                  type="character",default="value",short="-vval")
p <- add_argument(p, "--varname.opt",
     help="additional optional variables to be written on the output (out)",
                  type="character",default=NA,nargs=Inf,short="-vopt")
p<- add_argument(p, "--varname.prid",
                 help="name for the provider identifier (out)",
                 type="character",default="prid",short="-vprid")
p<- add_argument(p, "--varname.dqc",
                 help="name for the data quality control flag (out)",
                 type="character",default="dqc",short="-vdqc")
p<- add_argument(p, "--varname.sct",
            help="name for the spatial consistency test returned value (out)",
                 type="character",default="sct",short="-vsct")
p<- add_argument(p, "--varname.rep",
            help="name for the coefficient of representativeness (out)",
                  type="character",default="rep",short="-vrep")
# geographical parameters
p <- add_argument(p, "--spatconv",
                  help="flag for conversion of spatial coordinates before running the data quality checks",
                  flag=T,short="-c")
p <- add_argument(p, "--proj4from",
                  help="proj4 string for the original coordinate reference system",
                  type="character",
                  default="+proj=longlat +datum=WGS84",short="-pf")
p <- add_argument(p, "--proj4to",
                  help="proj4 string for the coordinate reference system where the DQC is performed",
                  type="character",
                  default="+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0",
#                  default="+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06",
                  short="-pt")
# parameter for the Box-Cox transformation
p <- add_argument(p, "--boxcox.lambda",
                  help="parameter used in the Box-Cox transformation",
                  type="numeric",default=0.5,short="-l")
# metadata check
p <- add_argument(p, "--zmin",help="minimum allowed elevation in the domain [m amsl]",
                  type="numeric",default=0,short="-z")
p <- add_argument(p, "--zmax",help="maximum allowed elevation in the domain [m amsl]",
                  type="numeric",default=2500,short="-Z")
# Plausibility check
p <- add_argument(p, "--tmax",help="maximum allowed value [mm]",
                  type="numeric",default=90,short="-TP")
# Cliamtological check
# default based on Norwegian hourly precipitation from 2010-2017
# (threshold=seasonal_maximum * 1.5 (with some subjective adjustments)
p <- add_argument(p, "--tmax.clim",
                  help="maximum allowed value [mm]",
                  type="numeric",nargs=12,short="-TC",
                  default=c(20,20,30,30,50,70,70,70,40,40,40,20))
p <- add_argument(p, "--month.clim",help="month (number 1-12)",
                  type="numeric",short="-mC",
#                  default=as.numeric(format(Sys.time(), "%m")))
                  default=NA)
# Check for isolated dry observations
p<-add_argument(p, "--i.isodry",
                help="number of iterations",
                type="numeric",default=1,short="-iI")
p<-add_argument(p, "--rr.isodry",
                help="precipitation/no-precipitation threshold [mm]",
                type="numeric",default=1,short="-rrI")
p<-add_argument(p, "--pmax.isodry",
    help="number of observations defining the neighbourhood to consider",
                type="numeric",default=20,short="-pI")
p<-add_argument(p, "--dmax.isodry",
    help="distance used to define the neighbourhood to consider [m]",
                type="numeric",default=150000,short="-dxI")
p<-add_argument(p, "--dmin.isodry",
    help="minimum distance between a dry observations surrounded only by wet observations and the closest observations to be considered good [m]",
                type="numeric",default=150000,short="-dI")
# check for dry observations by comparing against a first-guess field
p<- add_argument(p, "--file.fg",
                 help="first-guess file",
                 type="character",default=NA,short="-fF")
p<- add_argument(p, "--proj4.fg",
                 help="first-guess variable in netCDF file",
                 type="character",
                 default="+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0",
                 short="-pF")
p<- add_argument(p, "--varname.fg",
                 help="first-guess variable in netCDF file",
                 type="character",
                 default="lwe_precipitation_rate",short="-fF")
p<-add_argument(p, "--dr.fg",
                help="distance for the check against the first-guess field [m]",
                type="numeric",default=1500,short="-dF")
p<-add_argument(p, "--rr.fg",
                help="precipitation/no-precipitation threshold [mm]",
                type="numeric",default=.1,short="-rrF")
# Check for isolated wet observations
p<-add_argument(p, "--i.isowet",
                help="number of iterations",
                type="numeric",default=1,short="-iW")
p<-add_argument(p, "--rr.isowet",
                help="precipitation/no-precipitation threshold [mm]",
                type="numeric",default=1,short="-rrW")
p<-add_argument(p, "--pmax.isowet",
    help="number of observations defining the neighbourhood to consider",
                type="numeric",default=20,short="-pW")
p<-add_argument(p, "--dmax.isowet",
    help="distance used to define the neighbourhood to consider [m]",
                type="numeric",default=150000,short="-dxW")
p<-add_argument(p, "--dmin.isowet",
    help="minimum distance between a wet observations surrounded only by dry observations and the closest observations to be considered good [m]",
                type="numeric",default=100000,short="-dW")
# Buddy-check
p <- add_argument(p, "--dr.buddy",
                  help="perform the buddy-check in a dr-by-dr square-box around each observation [m]",
                  type="numeric",default=3000,short="-dB")
p <- add_argument(p, "--i.buddy",
                  help="number of buddy-check iterations",
                  type="integer",default=1,short="-iB")
p <- add_argument(p, "--thr.buddy",
                  help="buddy-check threshold. flag observation if: abs(obs-pred)/st_dev > thr.buddy",
                  type="numeric",default=3,short="-thB")
p <- add_argument(p, "--n.buddy",
                  help="minimum number of neighbouring observations to perform the buddy-check",
                  type="integer",default=5,short="-nB")
p <- add_argument(p, "--dz.buddy",
                  help="maximum allowed range of elevation in a square-box to perform the buddy-check (i.e. no check if elevation > dz.buddy)",
                  type="numeric",default=1500,short="-zB")
# isolated stations
p <- add_argument(p, "--dr.isol",
                  help="check for the number of observation in a dr-by-dr square-box around each observation [m]",
                  type="numeric",default=25000,short="-dI")
p <- add_argument(p, "--n.isol",
                  help="threshold (number of neighbouring observations) for the identification of isolated observations.",
                  type="integer",default=10,short="-nI")
# spatial consistency test
p <- add_argument(p, "--grid.sct",
                  help="nrow ncol (i.e. number_of_rows number_of_columns). used to define grid of boxes where the SCT is performed. SCT in each box is independent from the others",
                  type="integer",nargs=2,default=c(20,20),short="-gS")
p <- add_argument(p, "--i.sct",
                  help="number of SCT iterations",
                  type="integer",default=1,short="-iS")
p <- add_argument(p, "--n.sct",help="minimum number of stations in a box to run SCT",
                  type="integer",default=2,short="-nS")
p <- add_argument(p, "--dr.sct",
                  help="distance defining the neighbourhood used to obtain the background by averaging observations [m]",
                  type="numeric",default=20000,short="-dS")
p <- add_argument(p, "--nbg.sct",
                  help="minimum number of stations to compute the background value",
                  type="numeric",default=10,short="-nbgS")
p <- add_argument(p, "--rr.sct",
                  help="precipitation/no-precipitation threshold [mm]",
                  type="numeric",default=.1,short="-rrS")
p <- add_argument(p, "--dz.sct",
                  help="minimum range of elevation in a box to run SCT [m]",
                  type="numeric",default=0,short="-zS")
p <- add_argument(p, "--DhorMin.sct",
                  help="OI, minimum allowed value for the horizontal de-correlation lenght (of the background error correlation) [m]",
                  type="numeric",default=10000,short="-hS")
p <- add_argument(p, "--Dver.sct",
                  help="OI, vertical de-correlation lenght  (of the background error correlation) [m]",
                  type="numeric",default=1000,short="-vS")
p <- add_argument(p, "--eps2.sct",
                  help="OI, ratio between observation error variance and background error variance",
                  type="numeric",default=0.1,short="-eS")
p <- add_argument(p, "--thr.sct",
                  help="SCT threshold. flag observation if: (obs-Cross_Validation_pred)^2/(varObs+varCVpred) > thr.sct",
                  type="numeric",default=16,short="-tS")
p <- add_argument(p, "--laf.sct",
                  help="use land area fraction in the OI correlation function (0-100%)",
                  flag=T,short="-lS")
p <- add_argument(p, "--lafmin.sct",
                  help="land area fraction influence in the OI correlation function",
                  type="numeric",default=0.5,short="-lmS")
p <- add_argument(p, "--laf.file",
                  help="land area fraction file (netCDF in kilometric coordinates)",
                  type="character",default=NULL,short="-lfS")
p <- add_argument(p, "--proj4laf",
                  help="proj4 string for the laf",
                  type="character",
                  default="+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0",
#                  default="+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06",
                  short="-pl")
# observation representativeness
p <- add_argument(p, "--mean.corep",
                  help="average coefficient for the observation representativeness",
                  type="numeric",default=NA,nargs=Inf,short="-avC")
p <- add_argument(p, "--min.corep",
     help="minimum value for the coefficient for the observation representativeness",
                  type="numeric",default=NA,nargs=Inf,short="-mnC")
p <- add_argument(p, "--max.corep",
     help="maximum value for the coefficient for the observation representativeness",
                  type="numeric",default=NA,nargs=Inf,short="-mxC")
# check elevation against dem
p <- add_argument(p, "--dem",
     help="check elevation against digital elevation model (dem)",
                  flag=T,short="-dm")
p <- add_argument(p, "--dz.dem",
     help="maximum allowed deviation between observation and dem elevations [m]",
                  type="numeric",default=500,short="-zD")
p <- add_argument(p, "--dem.fill",help="fill missing elevation with data from dem",
                  flag=T,short="-df")
p <- add_argument(p, "--dem.file",
     help="land area fraction file (netCDF in kilometric coordinates)",
                  type="character",default=NULL,short="-dmf")
p <- add_argument(p, "--proj4dem",help="proj4 string for the dem",
                  type="character",
                  default="+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0",
#     default="+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06",
                  short="-pd")
# blacklist
# specified by triple/pairs of numbers: either (lat,lon,IDprovider) OR (index,IDprovider)
p <- add_argument(p, "--blacklist.lat",
                  help="observation blacklist (latitude)",
                  type="numeric",default=NA,nargs=Inf,short="-bla")
p <- add_argument(p, "--blacklist.lon",
                  help="observation blacklist (longitude)",
                  type="numeric",default=NA,nargs=Inf,short="-blo")
p <- add_argument(p, "--blacklist.fll",
                  help="observation blacklist (ID provider)",
                  type="numeric",default=NA,nargs=Inf,short="-bfll")
p <- add_argument(p, "--blacklist.idx",
                  help="observation blacklist (position in input file)",
                  type="numeric",default=NA,nargs=Inf,short="-bix")
p <- add_argument(p, "--blacklist.fidx",
                  help="observation blacklist (ID provider)",
                  type="numeric",default=NA,nargs=Inf,short="-bfix")
# keep (keep them)
# specified by triple/pairs of numbers: either (lat,lon,IDprovider) OR (index,IDprovider)
p <- add_argument(p, "--keeplist.lat",
                  help="observation keeplist (latitude)",
                  type="numeric",default=NA,nargs=Inf,short="-kla")
p <- add_argument(p, "--keeplist.lon",
                  help="observation keeplist (longitude)",
                  type="numeric",default=NA,nargs=Inf,short="-klo")
p <- add_argument(p, "--keeplist.fll",
                  help="observation keeplist (ID provider)",
                  type="numeric",default=NA,nargs=Inf,short="-kfll")
p <- add_argument(p, "--keeplist.idx",
                  help="observation keeplist (position in input file)",
                  type="numeric",default=NA,nargs=Inf,short="-kix")
p <- add_argument(p, "--keeplist.fidx",
                  help="observation keeplist (ID provider)",
                  type="numeric",default=NA,nargs=Inf,short="-kfix")
#
argv <- parse_args(p)
#
#-----------------------------------------------------------------------------
# checks on input arguments
if (!file.exists(argv$input)) {
  print("ERROR: input file not found")
  print(argv$input)
  quit(status=1)
}
# more than oen input file
if (any(!is.na(argv$input.files))) {
  for (j in 1:length(argv$input.files)) {
    if (!file.exists(argv$input.files[j])) {
      print("ERROR: input file not found")
      print(argv$input.files[j])
      quit(status=1)
    }
  }
  argv$input.files<-c(argv$input,argv$input.files)
} else {
  argv$input.files<-argv$input
}
nfin<-length(argv$input.files)
if (argv$dem | argv$dem.fill) {
  if (!file.exists(argv$dem.file)) {
    print("ERROR: dem file not found")
    print(argv$dem.file)
    quit(status=1)
  }
}
if (argv$laf.sct) {
  if (!file.exists(argv$laf.file)) {
    print("ERROR: laf file not found")
    print(argv$laf.file)
    quit(status=1)
  }
}
# TODO: adapt the procedure for input data others than lat-lon
if (!argv$spatconv) {
  print("ERROR: \"--spatconv\" (-c) option must be used onthe command line")
  print("input metadata are expected to be lat-lon coordinates")
  print(" some DQC tests takes place in kilometric coordinates specified by the user")
  print("output is in lat-lon coordinates")
  quit(status=1)
}
if ( argv$laf.sct & argv$proj4to!=argv$proj4laf ) {
  print("ERROR: anomalies found in the proj4 strings:")
  print(paste("proj4 laf=",argv$proj4laf))
  print(paste("proj4  to=",argv$proj4to))
  print("they must be the same")
  quit(status=1)
}
if ( (argv$dem | argv$dem.fill) & argv$proj4to!=argv$proj4dem ) {
  print("ERROR: anomalies found in the proj4 strings:")
  print(paste("proj4 dem=",argv$proj4dem))
  print(paste("proj4  to=",argv$proj4to))
  print("they must be the same")
  quit(status=1)
}
if (argv$laf.sct & (argv$dem | argv$dem.fill) &
    (argv$proj4laf!=argv$proj4dem) |
     !is.character(argv$proj4laf)  |
     !is.character(argv$proj4dem) ) {
  print("ERROR: anomalies found in the proj4 strings:")
  print(paste("proj4 laf=",argv$proj4laf))
  print(paste("proj4 dem=",argv$proj4dem))
  quit(status=1)
}
if (argv$laf.sct & (argv$dem | argv$dem.fill))
  suppressPackageStartupMessages(library("ncdf4")) 
#
if (!is.na(argv$month.clim) & (argv$month.clim<1 | argv$month.clim>12)) {
  print("ERROR: month number is wrong:")
  print(paste("month number=",argv$month.clim))
  quit(status=1)
}
# blacklist
if (any(!is.na(argv$blacklist.lat)) | 
    any(!is.na(argv$blacklist.lon)) |
    any(!is.na(argv$blacklist.fll)) ) {
  if ( (length(argv$blacklist.lat)!=length(argv$blacklist.lon))  |
       (length(argv$blacklist.lat)!=length(argv$blacklist.fll))  |
       (any(is.na(argv$blacklist.fll))) | 
       (any(is.na(argv$blacklist.lat))) | 
       (any(is.na(argv$blacklist.lon))) ) {
    print("ERROR in the blacklist definition, must have same number of lat,lon,IDprovider points")
    print(paste("lat number=",argv$blacklist.lat))
    print(paste("lon number=",argv$blacklist.lon))
    print(paste("ID provider number=",argv$blacklist.fll))
    quit(status=1)
  }
}
if (any(!is.na(argv$blacklist.idx)) | 
    any(!is.na(argv$blacklist.fidx)) ) {
  if ( (length(argv$blacklist.idx)!=length(argv$blacklist.fidx))  |
       (any(is.na(argv$blacklist.idx))) | 
       (any(is.na(argv$blacklist.fidx))) ) {
    print("ERROR in the blacklist definition, must have same number of index and IDprovider points")
    print(paste("index number=",argv$blacklist.idx))
    print(paste("ID provider number=",argv$blacklist.fidx))
    quit(status=1)
  }
}
# keeplist
if (any(!is.na(argv$keeplist.lat)) | 
    any(!is.na(argv$keeplist.lon)) |
    any(!is.na(argv$keeplist.fll)) ) {
  if ( (length(argv$keeplist.lat)!=length(argv$keeplist.lon))  |
       (length(argv$keeplist.lat)!=length(argv$keeplist.fll))  |
       (any(is.na(argv$keeplist.fll))) | 
       (any(is.na(argv$keeplist.lat))) | 
       (any(is.na(argv$keeplist.lon))) ) {
    print("ERROR in the keeplist definition, must have same number of lat,lon,IDprovider points")
    print(paste("lat number=",argv$keeplist.lat))
    print(paste("lon number=",argv$keeplist.lon))
    print(paste("ID provider number=",argv$keeplist.fll))
    quit(status=1)
  }
}
if (any(!is.na(argv$keeplist.idx)) | 
    any(!is.na(argv$keeplist.fidx)) ) {
  if ( (length(argv$keeplist.idx)!=length(argv$keeplist.fidx))  |
       (any(is.na(argv$keeplist.idx))) | 
       (any(is.na(argv$keeplist.fidx))) ) {
    print("ERROR in the keeplist definition, must have same number of index and IDprovider points")
    print(paste("index number=",argv$keeplist.idx))
    print(paste("ID provider number=",argv$keeplist.fidx))
    quit(status=1)
  }
}
# observation representativeness
if (any(is.na(argv$mean.corep)) | 
    any(is.na(argv$min.corep))  | 
    any(is.na(argv$max.corep)) ) {
  argv$mean.corep<-rep(1,nfin)
  argv$min.corep<-rep(.5,nfin)
  argv$max.corep<-rep(2,nfin)
}
if ( (length(argv$mean.corep)!=nfin) |
     (length(argv$min.corep)!=nfin)  |
     (length(argv$max.corep)!=nfin) ) {
  print("ERROR in the definitions of the coefficient for observation representativeness, there must be one coefficient for each input file (i.e. provider):")
  print(paste("number of input files=",nfin))
  print(paste("lenght of vector for mean value=",length(argv$mean.corep)))
  print(paste("lenght of vector for  min value=",length(argv$min.corep)))
  print(paste("lenght of vector for  max value=",length(argv$max.corep)))
  quit(status=1)
}
#
#-----------------------------------------------------------------------------
if (argv$verbose | argv$debug) print(">> POSEIDON <<")
#
#-----------------------------------------------------------------------------
# constants
nometa.code<-1
p.code<-2
clim.code<-3
isodry.code<-4
fg.code<-4
isowet.code<-5
buddy.code<-6
sct.code<-7
dem.code<-8
isol.code<-9
black.code<-100
keep.code<-200
#
#-----------------------------------------------------------------------------
# read data
first<-T
for (f in 1:nfin) {
  datain<-read.table(file=argv$input.files[f],header=T,sep=argv$separator,
                      stringsAsFactors=F,strip.white=T)
  varidx<-match(c(argv$varname.lat,
                  argv$varname.lon,
                  argv$varname.elev,
                  argv$varname.val),
                names(datain))
  if (any(is.na(varidx))) {
    print("ERROR in the specification of the variable names")
    print(paste("  latitutde=",argv$varname.lat))
    print(paste("  longitude=",argv$varname.lon))
    print(paste("  elevation=",argv$varname.elev))
    print(paste("temperature=",argv$varname.val))
    print("header of input file:")
    print(argv$input.files[f])
    print(names(datain))
    quit(status=1)
  }
  datatmp<-data.frame(datain[,varidx])
  names(datatmp)<-c("lat","lon","elev","value")
  datatmp$lat<-suppressWarnings(as.numeric(datatmp$lat))
  datatmp$lon<-suppressWarnings(as.numeric(datatmp$lon))
  datatmp$elev<-suppressWarnings(as.numeric(datatmp$elev))
  auxz<-datatmp$elev
  datatmp$value<-suppressWarnings(as.numeric(datatmp$value))
  ndatatmp<-length(datatmp$lat)
  if (ndatatmp==0) next
  # set provider id
  datatmp$prid<-rep(f,ndatatmp)
  aux<-rep(NA,length=ndatatmp)
  if (any(!is.na(argv$blacklist.idx)) & any(argv$blacklist.fidx==f)) {
    aux[argv$blacklist.idx[which(argv$blacklist.fidx==f)]]<-black.code  
  }
  if (any(!is.na(argv$blacklist.lat)) & any(argv$blacklist.fll==f)) {
    out<-apply(cbind(argv$blacklist.lon[argv$blacklist.fll==f],
                     argv$blacklist.lat[argv$blacklist.fll==f])
               ,FUN=setCode_lonlat,MARGIN=1,code=black.code)
  }
  if (any(!is.na(argv$keeplist.idx)) & any(argv$keeplist.fidx==f)) {
    aux[argv$keeplist.idx[which(argv$keeplist.fidx==f)]]<-keep.code  
  }
  if (any(!is.na(argv$keeplist.lat)) & any(argv$keeplist.fll==f)) {
    out<-apply(cbind(argv$keeplist.lon[argv$keeplist.fll==f],
                     argv$keeplist.lat[argv$keeplist.fll==f])
               ,FUN=setCode_lonlat,MARGIN=1,code=keep.code)
  }
  if (first) {
    data<-datatmp
    first<-F
    z <- auxz
    dqcflag<-aux
    sctpog<-rep(NA,length=ndatatmp)
    corep<-rep(NA,length=ndatatmp)
    if (any(!is.na(argv$varname.opt))) {
      varidx.opt<-match(argv$varname.opt,
                        names(datain))
      dataopt<-array(data=NA,dim=c(ndatatmp,
                                   length(argv$varname.opt)))
      dataopt<-datain[,varidx.opt[which(!is.na(varidx.opt))],drop=F]
    }
  } else {
    data<-rbind(data,datatmp)
    dqcflag<-c(dqcflag,aux)
    z <- c(z, auxz)
    sctpog<-c(sctpog,rep(NA,length=ndatatmp))
    corep<-c(corep,rep(NA,length=ndatatmp))
    if (any(!is.na(argv$varname.opt))) {
      varidx.opt.check<-match(argv$varname.opt,
                        names(datain))
      if (any(varidx.opt.check!=varidx.opt | 
              (is.na(varidx.opt.check) & !is.na(varidx.opt)) |
              (is.na(varidx.opt) & !is.na(varidx.opt.check)) )) {
        print("ERROR the header of file")
        print(argv$input.files[f])
        print("is different from the header of the first file")
        print(argv$input.files[1])
        quit(status=1)
      }
      dataopt<-rbind(dataopt,
        datain[,varidx.opt[which(!is.na(varidx.opt))],drop=F])
    }
  }
}
rm(datatmp)
ndata<-length(data$lat)
if (ndata==0) {
  print("input file is empty")
  quit(status=0)
}
if (argv$verbose | argv$debug) {
  print(paste("number of observations=",ndata))
  if (any(!is.na(argv$blacklist.idx)) | any(!is.na(argv$blacklist.lat)))
    print(paste("number of blacklisted observations=",
          length(which(dqcflag==black.code))) )
  if (any(!is.na(argv$keeplist.idx)) | any(!is.na(argv$keeplist.lat)))
    print(paste("number of keeplisted  observations=",
          length(which(dqcflag==keep.code))) )
  if (nfin>1) {
    for (f in 1:nfin) { 
      print(paste("  number of observations provider",f,"=",
            length(which(data$prid==f))))
      if (any(!is.na(argv$blacklist.idx)) | any(!is.na(argv$blacklist.lat)))
        print(paste("  number of blacklisted observations provider",f,"=",
              length(which(data$prid==f & dqcflag==black.code))) )
      if (any(!is.na(argv$keeplist.idx)) | any(!is.na(argv$keeplist.lat)))
        print(paste("  number of keeplisted  observations provider",f,"=",
              length(which(data$prid==f & dqcflag==keep.code))) )
    }
  }
  print("+---------------------------------+")
}
#
#-----------------------------------------------------------------------------
# coordinate transformation
if (argv$spatconv) {
  if (argv$debug) print("spatial conversion required")
  coord<-SpatialPoints(cbind(data$lon,data$lat),
                       proj4string=CRS(argv$proj4from))
  coord.new<-spTransform(coord,CRS(argv$proj4to))
  xy.new<-coordinates(coord.new)
  x<-round(xy.new[,1],0)
  y<-round(xy.new[,2],0)
  rm(xy.new)
  xp<-expand.grid(c(argv$lonmin,argv$lonmax),c(argv$latmin,argv$latmax))
  coord<-SpatialPoints(xp,
                       proj4string=CRS(argv$proj4from))
  coord.new<-spTransform(coord,CRS(argv$proj4to))
  # define the extent for the SCT grid
  e<-extent(coord.new)
  xl<-e[1:2]
  yl<-e[3:4]
} else {
  x<-data$lon
  y<-data$lat
  xl<-c(argv$lonmin,argv$lonmax)
  yl<-c(argv$latmin,argv$latmax)
  e<-extent(c(xl,yl))
}
#
#-----------------------------------------------------------------------------
# Read (optional) geographical information
if (argv$dem | argv$dem.fill) {
  # read netCDF and load into memory
  rdem0<-raster(argv$dem.file)
  crs(rdem0)<-argv$proj4dem
  rdem<-raster(rdem0)
  rdem[]<-getValues(rdem0)
  rm(rdem0)
  zdem<-extract(rdem,cbind(x,y))
  rdem.xmn<-xmin(rdem)
  rdem.xmx<-xmax(rdem)
  rdem.ymn<-ymin(rdem)
  rdem.ymx<-ymax(rdem)
  # fill missing elevation with dem
  if (argv$dem.fill) {
    iz<-which(is.na(z) & !is.na(zdem))
    z[iz]<-zdem[iz]
    rm(iz)
  }  
}
if (argv$laf.sct) {
  rlaf0<-raster(argv$laf.file)
  crs(rlaf0)<-argv$proj4laf
  rlaf<-raster(rlaf0)
  rlaf[]<-getValues(rlaf0)
  rm(rlaf0)
  laf<-extract(rlaf,cbind(x,y))/100.
} else {
  # use a fake laf
  laf<-rep(1,ndata)
}
#
#-----------------------------------------------------------------------------
# test for no metadata 
# use only (probably) good observations
ix<-which(is.na(dqcflag) | dqcflag==keep.code)
if (length(ix)>0) {
  meta<-!is.na(data$lat[ix]) & 
        !is.na(data$lon[ix]) &
        !is.na(z[ix]) & z[ix]>=argv$zmin & z[ix]<=argv$zmax &
        !is.na(data$value[ix]) &
        !is.na(laf[ix])
  if (any(!meta)) dqcflag[ix[which(!meta)]]<-nometa.code
} else {
  print("no valid observations left, no metadata check")
}
if (argv$verbose | argv$debug) {
  print("test for no metdata")
#  print(paste(data$lat[which(!meta)],data$lon[which(!meta)],
#              z[which(!meta)],data$value[which(!meta)]))
  print(paste("# observations lacking metadata=",length(which(!meta))))
  print("+---------------------------------+")
}
#
#-----------------------------------------------------------------------------
# plausibility test
ix<-which( (is.na(dqcflag) | dqcflag==keep.code) &
           data$value<0 &
           data$value>argv$tmax)
if (length(ix)>0) dqcflag[ix]<-p.code
if (argv$verbose | argv$debug) {
  print("plausibility test")
  print(paste("# suspect observations=",length(ix)))
  print("+---------------------------------+")
}
#
#-----------------------------------------------------------------------------
# climatological check 
# use only (probably) good observations
if (!is.na(argv$month.clim)) {
  ix<-which(is.na(dqcflag))
  if (length(ix)>0) {
    sus<-which(data$value[ix]<0 | 
               data$value[ix]>argv$tmax.clim[argv$month.clim] )
    # set dqcflag
    if (length(sus)>0) dqcflag[ix[sus]]<-clim.code
  } else {
    print("no valid observations left, no climatological check")
  }
  if (argv$verbose | argv$debug) {
    print(paste("climatological test (month=",argv$month.clim,")",sep=""))
    print(paste("# suspect observations=",length(which(dqcflag==clim.code))))
    print("+---------------------------------+")
  }
}
#
#-----------------------------------------------------------------------------
# check for isolated dry observations 
options(warn=1,scipen = 999)
if (argv$verbose | argv$debug) nprev<-0
for (i in 1:argv$i.isodry) {
  # use only (probably) good observations
  ix<-which(is.na(dqcflag) | dqcflag==keep.code)
  if (length(ix)>0) {
    t0a<-Sys.time()
    obs<-data.frame(x[ix],y[ix],data$value[ix])
    names(obs)<-c("x","y","yo")
    aux<-RR1_dqc_isolatedDry(obs,
                             rrinf=argv$rr.isodry,
                             pmax=argv$pmax.isodry,
                             dmax=argv$dmax.isodry,
                             n.sector=16,
                             dmin.dry=argv$dmin.isodry)
    sus<-which(aux!=0)
    # set dqcflag
    if (length(sus)>0) dqcflag[ix[sus]]<-isodry.code
  } else {
    print("no valid observations left, no check for isolated dry observations")
  }
  if (argv$verbose | argv$debug) {
    t1a<-Sys.time()
    print(paste("check for isolated dry observations, iteration=",i,
                "/time",round(t1a-t0a,1),attr(t1a-t0a,"unit")))
    ncur<-length(which(dqcflag==isodry.code))
    print(paste("# suspect observations=",ncur-nprev))
    nprev<-length(which(dqcflag==isodry.code))
  }
}
#
if (argv$debug) {
  if (!dir.exists(argv$debug.dir)) 
    dir.create(argv$debug.dir,showWarnings=F,recursive=T)
  sus<-which(dqcflag==isodry.code)
  for (j in 1:length(sus)) {
    i<-sus[j]
    f<-file.path(argv$debug.dir,
         paste("poseidon_dry_",
               formatC(i,width=2,flag="0"),
               ".png",sep=""))
    susi<-which(dqcflag==isodry.code)
    png(file=f,width=800,height=800)
    xmnj<-x[i]-30000
    xmxj<-x[i]+30000
    ymnj<-y[i]-30000
    ymxj<-y[i]+30000
    ee<-extent(xmnj,xmxj,ymnj,ymxj)
    plot(x,y,
         xlim=c(xmnj,xmxj),
         ylim=c(ymnj,ymxj),
         main="",
         xlab="",
         ylab="",cex=2, col="white")
    if (file.exists(argv$dem.file) &
        xmnj>=rdem.xmn & xmxj<=rdem.xmx & 
        ymnj>=rdem.ymn & ymxj<=rdem.ymx) {
      dem0<-crop(rdem,ee)
      image(dem0,add=T,
            breaks=c(0,10,25,50,100,250,500,750,1000,
                     1250,1500,1750,2000,2500,3000),
            col=gray.colors(14))
    }
    susi<-which(dqcflag==isodry.code)
    wet<-which(data$value>=argv$rr.isodry & is.na(dqcflag))
    dry<-which(data$value<argv$rr.isodry & is.na(dqcflag))
    points(x[wet],y[wet],pch=19,col="blue",cex=2)
    points(x[dry],y[dry],pch=15,col="orange",cex=2)
    points(x[susi],y[susi],pch=17,col="red",cex=2)
    dev.off()
  }
}
if (argv$verbose | argv$debug) 
  print("+---------------------------------+")
options(warn = 2, scipen = 999)
#
#-----------------------------------------------------------------------------
# check for dry observations by comparing against a first-guess field
if (!is.na(argv$file.fg)) {
  t0a<-Sys.time()
  ix<-which(is.na(dqcflag) | dqcflag==keep.code)
  if (length(ix)>0) {
    if (!file.exists(argv$file.fg)) {
      print("first-guess file not found")
      print(argv$file.fg)
    } else {
      xtot<-x[ix]
      ytot<-y[ix]
      ttot<-data$value[ix]
      rfg0<-raster(argv$file.fg)
      crs(rfg0)<-argv$proj4.fg
      rfg<-raster(rfg0)
      rfg[]<-getValues(rfg0)
      rm(rfg0)
      rfg.xmn<-xmin(rfg)
      rfg.xmx<-xmax(rfg)
      rfg.ymn<-ymin(rfg)
      rfg.ymx<-ymax(rfg)
      rrfg<-extract(rfg,cbind(xtot,ytot),buffer=argv$dr.fg,fun=max,na.rm=T)
      sus<-which(rrfg>=argv$rr.fg & ttot<argv$rr.fg)
      if (length(sus)>0) dqcflag[ix[sus]]<-fg.code
    }
  } else {
    print("no valid observations left, no check for dry observations against first-guess")
  }
  #
  if (argv$verbose | argv$debug) {
    t1a<-Sys.time()
    print(paste("check for dry observations against a first-guess",
                "/time",round(t1a-t0a,1),attr(t1a-t0a,"unit")))
    ncur<-length(which(dqcflag==fg.code))
    print(paste("# suspect observations=",ncur))
    if (argv$debug) {
      if (!dir.exists(argv$debug.dir)) 
      dir.create(argv$debug.dir,showWarnings=F,recursive=T)
      f<-file.path(argv$debug.dir,
           paste("poseidon_fg.png",sep=""))
      susi<-which(dqcflag[ix]==fg.code)
      png(file=f,width=800,height=800)
      plot(ttot,rrfg,
#           xlim=c(min(c(ttot,rrfg),na.rm=T),
#                  max(c(ttot,rrfg),na.rm=T)),
#           ylim=c(min(c(ttot,rrfg),na.rm=T),
#                  max(c(ttot,rrfg),na.rm=T)),
           xlim=c(0,6),
           ylim=c(0,6),
#           main=paste("/ #sus=",length(susi)),
           main="",
           xlab="Observations (mm)",
           ylab="First guess (mm)",col="white" )
      points(ttot,rrfg,pch=19,col="black")
      lines(-100:100,-100:100,col="gray")
      points(ttot[susi],rrfg[susi],pch=19,col="red")
      abline(h=seq(-1000,10000,by=2.5),col="gray",lty=2)
      abline(v=seq(-1000,10000,by=2.5),col="gray",lty=2)
      abline(h=0,lwd=2,col="gray",lty=1)
      abline(v=0,lwd=2,col="gray",lty=1)
      dev.off()
      sus<-which(dqcflag==fg.code)
      for (j in 1:length(sus)) {
        i<-sus[j]
        f<-file.path(argv$debug.dir,
             paste("poseidon_fg_",
                   formatC(i,width=2,flag="0"),
                   ".png",sep=""))
        susi<-which(dqcflag==fg.code)
        png(file=f,width=800,height=800)
        xmnj<-x[i]-30000
        xmxj<-x[i]+30000
        ymnj<-y[i]-30000
        ymxj<-y[i]+30000
        ee<-extent(xmnj,xmxj,ymnj,ymxj)
        plot(x,y,
             xlim=c(xmnj,xmxj),
             ylim=c(ymnj,ymxj),
             main="",
             xlab="",
             ylab="",cex=2,col="white")
        br<-c(0,.1,.5,1,2,3,4,5,6,7,8,9,10,11,3000)
        col<-c("beige",rev(rainbow((length(br)-2))))
        if (xmnj>=rfg.xmn & xmxj<=rfg.xmx & 
            ymnj>=rfg.ymn & ymxj<=rfg.ymx ) {
          fg0<-crop(rfg,ee)
          image(fg0,add=T,
                breaks=br,
                col=col)
        }
        susi<-which(dqcflag==fg.code)
        wet<-which(data$value>=argv$rr.fg & is.na(dqcflag))
        dry<-which(data$value<argv$rr.fg & is.na(dqcflag))
        points(x[wet],y[wet],pch=19,col="blue",cex=2)
        points(x[wet],y[wet],pch=1,col="black",cex=2)
        points(x[dry],y[dry],pch=15,col="orange",cex=2)
        points(x[dry],y[dry],pch=0,col="black",cex=2)
        points(x[susi],y[susi],pch=17,col="red",cex=2)
        points(x[susi],y[susi],pch=2,col="black",cex=2)
        dev.off()
      }
    }
  }
  if (argv$verbose | argv$debug) 
    print("+---------------------------------+")
}
#
#-----------------------------------------------------------------------------
# Check for isolated wet observations
options(warn=1,scipen = 999)
if (argv$verbose | argv$debug) nprev<-0
for (i in 1:argv$i.isowet) {
  # use only (probably) good observations
  ix<-which(is.na(dqcflag) | dqcflag==keep.code)
  if (length(ix)>0) {
    t0a<-Sys.time()
    obs<-data.frame(x[ix],y[ix],data$value[ix])
    names(obs)<-c("x","y","yo")
    aux<-RR1_dqc_isolatedWet(obs,
                             rrinf=argv$rr.isowet,
                             pmax=argv$pmax.isowet,
                             dmax=argv$dmax.isowet,
                             n.sector=16,
                             dmin.wet=argv$dmin.isowet)
    sus<-which(aux!=0)
    # set dqcflag
    if (length(sus)>0) dqcflag[ix[sus]]<-isowet.code
  } else {
    print("no valid observations left, no check for isolated wet observations")
  }
  if (argv$verbose | argv$debug) {
    t1a<-Sys.time()
    print(paste("check for isolated wet observations, iteration=",i,
                "/time",round(t1a-t0a,1),attr(t1a-t0a,"unit")))
    ncur<-length(which(dqcflag==isowet.code))
    print(paste("# suspect observations=",ncur-nprev))
    nprev<-length(which(dqcflag==isowet.code))
  }
}
#
if (argv$debug) {
  if (!dir.exists(argv$debug.dir)) 
    dir.create(argv$debug.dir,showWarnings=F,recursive=T)
  sus<-which(dqcflag==isowet.code)
  for (j in 1:length(sus)) {
    i<-sus[j]
    f<-file.path(argv$debug.dir,
         paste("poseidon_wet_",
               formatC(i,width=2,flag="0"),
               ".png",sep=""))
    susi<-which(dqcflag==isowet.code)
    png(file=f,width=800,height=800)
    xmnj<-x[i]-30000
    xmxj<-x[i]+30000
    ymnj<-y[i]-30000
    ymxj<-y[i]+30000
    ee<-extent(xmnj,xmxj,ymnj,ymxj)
    plot(x,y,
         xlim=c(xmnj,xmxj),
         ylim=c(ymnj,ymxj),
         main="",
         xlab="",
         ylab="",cex=2,col="white")
    if (file.exists(argv$dem.file) & 
        xmnj>=rdem.xmn & xmxj<=rdem.xmx & 
        ymnj>=rdem.ymn & ymxj<=rdem.ymx ) {
      dem0<-crop(rdem,ee)
      image(dem0,add=T,
            breaks=c(0,10,25,50,100,250,500,750,1000,
                     1250,1500,1750,2000,2500,3000),
            col=gray.colors(14))
    }
    susi<-which(dqcflag==isowet.code)
    wet<-which(data$value>=argv$rr.isowet & is.na(dqcflag))
    dry<-which(data$value<argv$rr.isowet & is.na(dqcflag))
    points(x[wet],y[wet],pch=19,col="blue",cex=2)
    points(x[dry],y[dry],pch=15,col="orange",cex=2)
    points(x[susi],y[susi],pch=17,col="cyan",cex=2)
    dev.off()
  }
}
if (argv$verbose | argv$debug) 
  print("+---------------------------------+")
options(warn = 2, scipen = 999)
#
#-----------------------------------------------------------------------------
# buddy check 
#  compare each observation against the average of neighbouring observations 
if (argv$verbose | argv$debug) nprev<-0
for (i in 1:argv$i.buddy) {
  # use only (probably) good observations
  ix<-which(is.na(dqcflag) | dqcflag==keep.code)
  if (length(ix)>0) {
    t0a<-Sys.time()
    # define global 1D vector used in statSpat (1D for fast access)
    xtot<-x[ix]
    ytot<-y[ix]
    ztot<-z[ix]
    ttot<-boxcox(x=data$value[ix],lambda=argv$boxcox.lambda)
    # apply will loop over this 4D array
    ixyzt_tot<-cbind(1:length(xtot),xtot,ytot,ztot,ttot)
    stSp<-apply(ixyzt_tot,FUN=statSpat,MARGIN=1,drmin=argv$dr.buddy,tcor.flag=F)
    # probability of gross error
    stSp[4,][which(stSp[4,]<0.1)]<-0.1
    pog<-abs(ttot-stSp[3,])/stSp[4,]
    # suspect if: 
    sus<-which(pog>argv$thr.buddy & 
               stSp[1,]>argv$n.buddy & 
               stSp[2,]<argv$dz.buddy &
               is.na(dqcflag[ix]))
    # set dqcflag
    if (length(sus)>0) dqcflag[ix[sus]]<-buddy.code
  } else {
    print("no valid observations left, no buddy check")
  }
  if (argv$verbose | argv$debug) {
    t1a<-Sys.time()
    print(paste("buddy-check, iteration=",i,
                "/time",round(t1a-t0a,1),attr(t1a-t0a,"unit")))
    ncur<-length(which(dqcflag==buddy.code))
    print(paste("# suspect observations=",ncur-nprev))
    nprev<-length(which(dqcflag==buddy.code))
  }
}
if (argv$verbose | argv$debug) { 
  if (argv$debug) {
    if (!dir.exists(argv$debug.dir)) 
    dir.create(argv$debug.dir,showWarnings=F,recursive=T)
    f<-file.path(argv$debug.dir,
         paste("poseidon_buddy_mean.png",sep=""))
    susi<-which(dqcflag[ix]==buddy.code)
    png(file=f,width=800,height=800)
    plot(ttot,stSp[3,],
         xlim=c(min(c(ttot,stSp[3,]),na.rm=T),
                max(c(ttot,stSp[3,]),na.rm=T)),
         ylim=c(min(c(ttot,stSp[3,]),na.rm=T),
                max(c(ttot,stSp[3,]),na.rm=T)),
#         xlim=c(0,6),
#         ylim=c(0,6),
#           main=paste("/ #sus=",length(susi)),
         main="",
         xlab="Observations (mm)",
         ylab="Average (mm)",col="white" )
    points(ttot,stSp[3,],pch=19,col="black")
    lines(-100:100,-100:100,col="gray")
    points(ttot[susi],stSp[3,][susi],pch=19,col="red")
    abline(h=seq(-1000,10000,by=2.5),col="gray",lty=2)
    abline(v=seq(-1000,10000,by=2.5),col="gray",lty=2)
    abline(h=0,lwd=2,col="gray",lty=1)
    abline(v=0,lwd=2,col="gray",lty=1)
    dev.off()
    f<-file.path(argv$debug.dir,
         paste("poseidon_buddy_sd.png",sep=""))
    susi<-which(dqcflag[ix]==buddy.code)
    png(file=f,width=800,height=800)
    plot(ttot,stSp[4,],
#         xlim=c(min(c(ttot,stSp[4,]),na.rm=T),
#                max(c(ttot,stSp[4,]),na.rm=T)),
#         ylim=c(min(c(ttot,stSp[4,]),na.rm=T),
#                max(c(ttot,stSp[4,]),na.rm=T)),
         xlim=c(-2,4),
         ylim=c(0,3),
#           main=paste("/ #sus=",length(susi)),
         main="",
         xlab="Observations (mm)",
         ylab="Standard deviation (mm)",col="white" )
    points(ttot,stSp[4,],pch=19,col="black")
#    lines(-100:100,-100:100,col="gray")
    points(ttot[susi],stSp[4,susi],pch=19,col="red")
    abline(h=seq(-1000,10000,by=.25),col="gray",lty=3)
    abline(h=seq(-1000,10000,by=1),col="gray",lty=2)
    abline(v=seq(-1000,10000,by=1),col="gray",lty=2)
    abline(h=0,lwd=2,col="blue",lty=1)
    abline(v=-2,lwd=2,col="blue",lty=1)
    dev.off()
    #
    susi<-which(dqcflag[ix]==buddy.code)
    for (j in 1:length(susi)) {
      i<-susi[j]
      f<-file.path(argv$debug.dir,
           paste("poseidon_buddy_",
                 formatC(i,width=5,flag="0"),
                 ".png",sep=""))
      png(file=f,width=800,height=800)
      xmnj<-xtot[i]-argv$dr.buddy
      xmxj<-xtot[i]+argv$dr.buddy
      ymnj<-ytot[i]-argv$dr.buddy
      ymxj<-ytot[i]+argv$dr.buddy
      ee<-extent(xmnj,xmxj,ymnj,ymxj)
      plot(xtot,xtot,
           xlim=c(xmnj-1*argv$dr.buddy,xmxj+1*argv$dr.buddy),
           ylim=c(ymnj-1*argv$dr.buddy,ymxj+1*argv$dr.buddy),
           main="",
           xlab="",
           ylab="",cex=2,col="white")
      brt<-c(0,0.1,0.5,1,2,3,4,5,6,7,8,9,10,3000)
      br<-boxcox(brt,
                 lambda=argv$boxcox.lambda)
      col<-c("beige",rev(rainbow((length(br)-2))))
      xmnj<-xtot[i]-2*argv$dr.buddy
      xmxj<-xtot[i]+2*argv$dr.buddy
      ymnj<-ytot[i]-2*argv$dr.buddy
      ymxj<-ytot[i]+2*argv$dr.buddy
      if (file.exists(argv$dem.file) & 
          xmnj>=rdem.xmn & xmxj<=rdem.xmx & 
          ymnj>=rdem.ymn & ymxj<=rdem.ymx ) {
        ed<-extent(xmnj,xmxj,ymnj,ymxj)
        dem0<-crop(rdem,ed)
        image(dem0,add=T,
              breaks=c(0,10,25,50,100,250,500,750,1000,
                       1250,1500,1750,2000,2500,3000),
              col=gray.colors(14))
      }
      for (c in 1:length(col)) {
        if (c==1) {
          legstr<-paste("<",brt[2],"mm",sep="")
        } else if (c==length(col)) {
          legstr<-c(legstr,paste(">=",brt[c],sep=""))
        } else {
          legstr<-c(legstr,paste("[",brt[c],",",brt[c+1],")",sep=""))
        }
        aux<-which(ttot>=br[c] & 
                   ttot<br[c+1] & 
                   is.na(dqcflag[ix]))
        if (length(aux)>0) {
          points(xtot[aux],ytot[aux],pch=19,col=col[c],cex=2)
          points(xtot[aux],ytot[aux],pch=1,col="black",cex=2)
        }
        aux<-which(ttot>=br[c] & 
                   ttot<br[c+1] & 
                   dqcflag[ix]==buddy.code)
        if (length(aux)>0) {
          points(xtot[aux],ytot[aux],pch=17,col=col[c],cex=2)
          points(xtot[aux],ytot[aux],pch=2,col="black",cex=2)
        }
      }
      plot(ee,add=T,lwd=3)
#      legend(x="bottomright",fill=rev(col),legend=rev(legstr),cex=1.5)
      dev.off()
    }
  }
  print("+---------------------------------+")
}
#
#-----------------------------------------------------------------------------
# SCT - Spatial Consistency Test
if (argv$verbose | argv$debug) nprev<-0
for (i in 1:argv$i.sct) {
  # use only (probably) good observations
  ix<-which(is.na(dqcflag) | dqcflag==keep.code)
  if (length(ix)>0) {
    t0a<-Sys.time()
    # create the grid for SCT. SCT is done independently in each box
    # NOTE: box size around 100Km should be ok
    if (i==1) {
      r<-raster(e,ncol=argv$grid.sct[2],nrow=argv$grid.sct[1])
      xy<-xyFromCell(r,1:ncell(r))
      xr<-xy[,1]
      yr<-xy[,2]
      ir<-1:ncell(r)
      r[]<-1:ncell(r)
    }
    # define global 1D vector used in statSpat (1D for fast access)
    xtot<-x[ix]
    ytot<-y[ix]
    ztot<-z[ix]
    ttot<-boxcox(x=data$value[ix],lambda=argv$boxcox.lambda)
    # apply will loop over this 4D array
    ixyzt_tot<-cbind(1:length(xtot),xtot,ytot,ztot,ttot)
    stSp<-apply(ixyzt_tot,
                FUN=statSpat,
                MARGIN=1,
                drmin=argv$dr.sct,
                tcor.flag=F)
    ixok<-which(stSp[1,]>argv$nbg.sct & 
                data$value[ix]>argv$rr.sct)
    if (length(ixok)>3) {
      # linear regression such that x=1/argv$boxcox.lambda corresponds 
      #  to y=1/argv$boxcox.lambda 
      yyy<-stSp[3,ixok]+1/argv$boxcox.lambda
      xxx<-ttot[ixok]+1/argv$boxcox.lambda
      m<-as.numeric(lm(yyy~xxx+0)$coefficients)
      rm(yyy,xxx)
      # background (large scale) value
      tbtot<-m*ttot+1/argv$boxcox.lambda*m-1/argv$boxcox.lambda
      tbtot[which(!is.na(stSp[3,]))]<-stSp[3,][which(!is.na(stSp[3,]))]
      # assign each station to the corresponding box
      itot<-extract(r,cbind(xtot,ytot))
      # count the number of observations in each box
      rnobs<-rasterize(cbind(xtot,ytot),r,ttot,fun=function(x,...)length(x))
      nr<-getValues(rnobs)
      # create the 4D array for the function call via apply
      ixyn<-cbind(ir,xr,yr,nr)
      # SCT within each (grid) box
      # NOTE: dqcflag is updated in "sct" function
      out<-apply(ixyn,
                 FUN=sct,MARGIN=1,
                 nmin=argv$n.sct,
                 dzmin=argv$dz.sct,
                 Dhmin=argv$DhorMin.sct,
                 Dz=argv$Dver.sct,
                 eps2=argv$eps2.sct,
                 T2=argv$thr.sct,
                 sus.code=sct.code)
    } else {
      print("not enough suitable observations left for the SCT")
    }
  } else {
    print("no valid observations left, no SCT")
  }
  if (argv$verbose | argv$debug) {
    t1a<-Sys.time()
    print(paste("SCT, iteration=",i,
                "/time",round(t1a-t0a,1),attr(t1a-t0a,"unit")))
    ncur<-length(which(dqcflag==sct.code))
    print(paste("# suspect observations=",ncur-nprev))
    nprev<-length(which(dqcflag==sct.code))
    if (argv$debug) {
      if (!dir.exists(argv$debug.dir)) 
      dir.create(argv$debug.dir,showWarnings=F,recursive=T)
      f<-file.path(argv$debug.dir,
           paste("poseidon_sctm_it_",formatC(i,width=2,flag="0"),".png",sep=""))
      susi<-which(dqcflag[ix]==sct.code)
      png(file=f,width=800,height=800)
      plot(ttot[ixok],stSp[3,ixok],
           xlim=c(min(c(ttot,stSp[3,ixok])),
                  max(c(ttot,stSp[3,ixok]))),
           ylim=c(min(c(ttot,stSp[3,ixok])),
                  max(c(ttot,stSp[3,ixok]))),
           main=paste("ang coeff=",round(m,3),
                      "/ #sus=",length(susi)),
           xlab="Observations (Box-Cox transformed)",
           ylab="Large-Scale value / Background" )
      aaa<-m*-100:100+1/argv$boxcox.lambda*m-1/argv$boxcox.lambda
      lines(-100:100,aaa,col="darkblue")
      lines(-100:100,-100:100,col="gray")
      points(ttot,m*ttot+1/argv$boxcox.lambda*m-1/argv$boxcox.lambda,col="blue")
      points(ttot[susi],m*ttot[susi]+1/argv$boxcox.lambda*m-1/argv$boxcox.lambda,pch=19,col="red")
      abline(h=seq(-1000,10000,by=10),col="gray",lty=2)
      abline(v=seq(-1000,10000,by=10),col="gray",lty=2)
      abline(h=0,lwd=2,col="gray",lty=1)
      abline(v=0,lwd=2,col="gray",lty=1)
      abline(h=-1/argv$boxcox.lambda,lwd=2,col="black",lty=1)
      abline(v=-1/argv$boxcox.lambda,lwd=2,col="black",lty=1)
      abline(v=ttot[susi],col="red")
      dev.off()
    }
  }
} # end SCT loop
if (argv$verbose | argv$debug) 
  print("+---------------------------------+")
#
# coefficient of observation representativeness
#-----------------------------------------------------------------------------
# corep has been set by function sct to the observation error variance
qmn<-0.25
qmx<-0.75
qav<-0.5
ix<-which(!is.na(corep) & (is.na(dqcflag) | dqcflag==keep.code)) 
if (length(ix)>0) {
  acorep<-abs(corep[ix])
  # ecdf(x)(x) here should give us something similar to rank(x)/length(x)
  qcorep<-ecdf(acorep)(acorep)
  if (any(qcorep<qmn)) qcorep[which(qcorep<qmn)]<-qmn
  if (any(qcorep>qmx)) qcorep[which(qcorep>qmx)]<-qmx
  for (f in 1:nfin) {
    if (any(data$prid[ix]==f)) {
      ip<-which(data$prid[ix]==f & qcorep<=qav)
      if (length(ip)>0)      
        corep[ix[ip]]<-argv$min.corep[f]+
          (argv$mean.corep[f]-argv$min.corep[f])*(qcorep[ip]-qmn)/(qav-qmn)
      ip<-which(data$prid[ix]==f & qcorep>qav)
      if (length(ip)>0)      
        corep[ix[ip]]<-argv$mean.corep[f]+
          (argv$max.corep[f]-argv$mean.corep[f])*(qcorep[ip]-qav)/(qmx-qav)
    } else {
      print(paste("provider ",f,": no valid first guess for obs-err-var",sep=""))
    }
  }
} else {
  print("no valid first guess for the observation error variance found")
}
#
#-----------------------------------------------------------------------------
# check for isolated dry observations 
options(warn=1,scipen = 999)
if (argv$verbose | argv$debug) nprev<-0
for (i in 1:argv$i.isodry) {
  # use only (probably) good observations
  ix<-which(is.na(dqcflag) | dqcflag==keep.code)
  if (length(ix)>0) {
    t0a<-Sys.time()
    obs<-data.frame(x[ix],y[ix],data$value[ix])
    names(obs)<-c("x","y","yo")
    aux<-RR1_dqc_isolatedDry(obs,
                             rrinf=argv$rr.isodry,
                             pmax=argv$pmax.isodry,
                             dmax=argv$dmax.isodry,
                             n.sector=16,
                             dmin.dry=argv$dmin.isodry)
    sus<-which(aux!=0)
    # set dqcflag
    if (length(sus)>0) dqcflag[ix[sus]]<-isodry.code
  } else {
    print("no valid observations left, no check for isolated dry observations")
  }
  if (argv$verbose | argv$debug) {
    t1a<-Sys.time()
    print(paste("check for isolated dry observations, iteration=",i,
                "/time",round(t1a-t0a,1),attr(t1a-t0a,"unit")))
    ncur<-length(which(dqcflag==isodry.code))
    print(paste("# suspect observations=",ncur-nprev))
    nprev<-length(which(dqcflag==isodry.code))
  }
}
#
if (argv$debug) {
  if (!dir.exists(argv$debug.dir)) 
    dir.create(argv$debug.dir,showWarnings=F,recursive=T)
  sus<-which(dqcflag==isodry.code)
  for (j in 1:length(sus)) {
    i<-sus[j]
    f<-file.path(argv$debug.dir,
         paste("poseidon_dry_",
               formatC(i,width=2,flag="0"),
               ".png",sep=""))
    susi<-which(dqcflag==isodry.code)
    png(file=f,width=800,height=800)
    xmnj<-x[i]-30000
    xmxj<-x[i]+30000
    ymnj<-y[i]-30000
    ymxj<-y[i]+30000
    ee<-extent(xmnj,xmxj,ymnj,ymxj)
    plot(x,y,
         xlim=c(xmnj,xmxj),
         ylim=c(ymnj,ymxj),
         main="",
         xlab="",
         ylab="",cex=2, col="white")
    if (file.exists(argv$dem.file) &
        xmnj>=rdem.xmn & xmxj<=rdem.xmx & 
        ymnj>=rdem.ymn & ymxj<=rdem.ymx) {
      dem0<-crop(rdem,ee)
      image(dem0,add=T,
            breaks=c(0,10,25,50,100,250,500,750,1000,
                     1250,1500,1750,2000,2500,3000),
            col=gray.colors(14))
    }
    susi<-which(dqcflag==isodry.code)
    wet<-which(data$value>=argv$rr.isodry & is.na(dqcflag))
    dry<-which(data$value<argv$rr.isodry & is.na(dqcflag))
    points(x[wet],y[wet],pch=19,col="blue",cex=2)
    points(x[dry],y[dry],pch=15,col="orange",cex=2)
    points(x[susi],y[susi],pch=17,col="red",cex=2)
    dev.off()
  }
}
if (argv$verbose | argv$debug) 
  print("+---------------------------------+")
options(warn = 2, scipen = 999)
#
#-----------------------------------------------------------------------------
# Check for isolated wet observations
options(warn=1,scipen = 999)
if (argv$verbose | argv$debug) nprev<-0
for (i in 1:argv$i.isowet) {
  # use only (probably) good observations
  ix<-which(is.na(dqcflag) | dqcflag==keep.code)
  if (length(ix)>0) {
    t0a<-Sys.time()
    obs<-data.frame(x[ix],y[ix],data$value[ix])
    names(obs)<-c("x","y","yo")
    aux<-RR1_dqc_isolatedWet(obs,
                             rrinf=argv$rr.isowet,
                             pmax=argv$pmax.isowet,
                             dmax=argv$dmax.isowet,
                             n.sector=16,
                             dmin.wet=argv$dmin.isowet)
    sus<-which(aux!=0)
    # set dqcflag
    if (length(sus)>0) dqcflag[ix[sus]]<-isowet.code
  } else {
    print("no valid observations left, no check for isolated wet observations")
  }
  if (argv$verbose | argv$debug) {
    t1a<-Sys.time()
    print(paste("check for isolated wet observations, iteration=",i,
                "/time",round(t1a-t0a,1),attr(t1a-t0a,"unit")))
    ncur<-length(which(dqcflag==isowet.code))
    print(paste("# suspect observations=",ncur-nprev))
    nprev<-length(which(dqcflag==isowet.code))
  }
}
#
if (argv$debug) {
  if (!dir.exists(argv$debug.dir)) 
    dir.create(argv$debug.dir,showWarnings=F,recursive=T)
  sus<-which(dqcflag==isowet.code)
  for (j in 1:length(sus)) {
    i<-sus[j]
    f<-file.path(argv$debug.dir,
         paste("poseidon_wet_",
               formatC(i,width=2,flag="0"),
               ".png",sep=""))
    susi<-which(dqcflag==isowet.code)
    png(file=f,width=800,height=800)
    xmnj<-x[i]-30000
    xmxj<-x[i]+30000
    ymnj<-y[i]-30000
    ymxj<-y[i]+30000
    ee<-extent(xmnj,xmxj,ymnj,ymxj)
    plot(x,y,
         xlim=c(xmnj,xmxj),
         ylim=c(ymnj,ymxj),
         main="",
         xlab="",
         ylab="",cex=2,col="white")
    if (file.exists(argv$dem.file) & 
        xmnj>=rdem.xmn & xmxj<=rdem.xmx & 
        ymnj>=rdem.ymn & ymxj<=rdem.ymx ) {
      dem0<-crop(rdem,ee)
      image(dem0,add=T,
            breaks=c(0,10,25,50,100,250,500,750,1000,
                     1250,1500,1750,2000,2500,3000),
            col=gray.colors(14))
    }
    susi<-which(dqcflag==isowet.code)
    wet<-which(data$value>=argv$rr.isowet & is.na(dqcflag))
    dry<-which(data$value<argv$rr.isowet & is.na(dqcflag))
    points(x[wet],y[wet],pch=19,col="blue",cex=2)
    points(x[dry],y[dry],pch=15,col="orange",cex=2)
    points(x[susi],y[susi],pch=17,col="cyan",cex=2)
    dev.off()
  }
}
if (argv$verbose | argv$debug) 
  print("+---------------------------------+")
options(warn = 2, scipen = 999)
#
#-----------------------------------------------------------------------------
# check elevation against dem 
if (argv$dem) {
  # use only (probably) good observations
  ix<-which(is.na(dqcflag))
  if (length(ix)>0) {
    ixna<-which(!is.na(z) & !is.na(zdem) & is.na(dqcflag))
    sus<-which(abs(z[ixna]-zdem[ixna])>argv$dz.dem)
    # set dqcflag
    if (length(sus)>0) dqcflag[ixna[sus]]<-dem.code
  }  else {
    print("no valid observations left, no dem check")
  }
  if (argv$verbose | argv$debug) {
    print(paste("# observations too far from dem=",length(which(dqcflag==dem.code))))
    print("+---------------------------------+")
  }
}
#
#-----------------------------------------------------------------------------
# check for isolated stations
# use only (probably) good observations
ix<-which(is.na(dqcflag))
if (length(ix)>0) {
  # define global 1D vector used in statSpat (1D for fast access)
  xtot<-x[ix]
  ytot<-y[ix]
  xy<-cbind(xtot,ytot)
  ns<-apply(xy,FUN=nstat,MARGIN=1,drmin=argv$dr.isol)
  sus<-which(ns<argv$n.isol)
  # set dqcflag
  if (length(sus)>0) dqcflag[ix[sus]]<-isol.code
} else {
  print("no valid observations left, no check for isolated observations")
}
if (argv$verbose | argv$debug) {
  print(paste("# isolated observations=",length(which(dqcflag==isol.code))))
  print("+---------------------------------+")
}
#
#-----------------------------------------------------------------------------
# observations not flagged are assumed to be good observaitons 
dqcflag[is.na(dqcflag)]<-0
if (argv$verbose | argv$debug) {
  print("summary")
  print(paste("# total suspect observations=",
              length(which(dqcflag!=0))," [",
              round(100*length(which(dqcflag!=0))/ndata,0),
              "%]",sep="") )
}
#
#-----------------------------------------------------------------------------
# Debug: create Figures 
if (argv$debug) {
  r<-raster(e,ncol=argv$grid.sct[2],nrow=argv$grid.sct[1])
  xy<-xyFromCell(r,1:ncell(r))
  xr<-xy[,1]
  yr<-xy[,2]
  ir<-1:ncell(r)
  r[]<-1:ncell(r)
  # apply will loop over this 4D array
  ixyzt<-cbind(1:length(x),x,y,z,data$value)
  # assign each station to the corresponding box
  itot<-extract(r,cbind(x,y))
  # count the number of observations in each box
  rnobs<-rasterize(cbind(x,y),r,data$value,fun=function(x,...)length(x))
  nr<-getValues(rnobs)
  # create the 4D array for the function call via apply
  ixyn<-cbind(ir,xr,yr,nr)
  out<-apply(ixyn,
             FUN=plotsummary,
             MARGIN=1)

}
#
#-----------------------------------------------------------------------------
# write the output file
varidx.out<-varidx
if (any(!is.na(argv$varname.opt))) 
  varidx.out<-c(varidx,varidx.opt[which(!is.na(varidx.opt))]) 
dataout<-array(data=NA,
               dim=c(length(data$lat),(length(varidx.out)+4)))
ord.varidx.out<-order(varidx.out)
str<-vector()
for (s in 1:length(ord.varidx.out)) {
  varidx.out.s<-varidx.out[ord.varidx.out[s]]
  pos.s<-which(varidx.out==varidx.out.s)
  if (pos.s>4) {
    posopt.s<-which(varidx.opt==varidx.out.s & !is.na(varidx.opt))
    posopt.nona.s<-which(varidx.opt[which(!is.na(varidx.opt))]==varidx.out.s)
    str[s]<-argv$varname.opt[posopt.s]
    dataout[,s]<-dataopt[,posopt.nona.s]
  } else if (pos.s==1) {
    str[s]<-"lat"
    dataout[,s]<-round(data$lat,5)
  } else if (pos.s==2) {
    str[s]<-"lon"
    dataout[,s]<-round(data$lon,5)
  } else if (pos.s==3) {
    str[s]<-"elev"
    dataout[,s]<-round(z,1)
  } else if (pos.s==4) {
    str[s]<-"value"
    dataout[,s]<-round(data$value,1)
  }
}
str[s+1]<-argv$varname.prid
dataout[,(s+1)]<-data$prid
str[s+2]<-argv$varname.dqc
dataout[,(s+2)]<-dqcflag
str[s+3]<-argv$varname.sct
dataout[,(s+3)]<-round(sctpog,2)
str[s+4]<-argv$varname.rep
dataout[,(s+4)]<-round(corep,5)
dataout<-as.data.frame(dataout,stringsAsFactors=F)
names(dataout)<-str
write.table(file=argv$output,
            dataout,
            quote=F,
            col.names=str,
            row.names=F,
            dec=".",
            sep=";")
#
#-----------------------------------------------------------------------------
# Normal exit
t1<-Sys.time()
if (argv$verbose | argv$debug) 
 print(paste("normal exit /time",round(t1-t0,1),attr(t1-t0,"unit")))
q(status=0)
