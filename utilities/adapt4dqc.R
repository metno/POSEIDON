#!/usr/bin/env Rscript
#
# + <adapt4dqc.R> adapt file format for Data Quality Control routines
#
#------------------------------------------------------------------------------
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("sp"))
suppressPackageStartupMessages(library("rgdal"))
#options(warn = 2, scipen = 999)
#
#------------------------------------------------------------------------------
# create parser object
p <- arg_parser("adapt4dqc")
# specify our desired options 
# by default ArgumentParser will add an help option 
p <- add_argument(p, "input",
                  help="input file",type="character")
p <- add_argument(p, "output",
                  help="output file",type="character",
                  default="output.txt")
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
p <- add_argument(p, "--proj4from",
                  help="proj4 string for the original coordinate reference system",
                  type="character",
                  default="+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0",
                  short="-pf")
p <- add_argument(p, "--proj4to",
                  help="proj4 string for the coordinate reference system where the DQC is performed",
                  type="character",
                  default="+proj=longlat +datum=WGS84",
#                  default="+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06",
                  short="-pt")
#
argv <- parse_args(p)
#
#------------------------------------------------------------------------------
# checks
if (!file.exists(argv$input)) {
  print("ERROR: input file not found")
  print(argv$input)
  quit(status=1)
}
#
#------------------------------------------------------------------------------
# 
datain<-read.table(file=argv$input,
                   header=T,sep=argv$separator,
                   stringsAsFactors=F,strip.white=T)
varidx<-match(c(argv$varname.lat,
                argv$varname.lon,
                argv$varname.elev,
                argv$varname.value),
              names(datain))
if (any(is.na(varidx))) {
  print("ERROR in the specification of the variable names")
  print(paste("  latitutde=",argv$varname.lat))
  print(paste("  longitude=",argv$varname.lon))
  print(paste("  elevation=",argv$varname.elev))
  print(paste("      value=",argv$varname.value))
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

if (any(!is.na(argv$varname.opt))) {
  varidx.opt<-match(argv$varname.opt,
                    names(datain))
  dataopttmp<-datain[,varidx.opt[which(!is.na(varidx.opt))],drop=F]
}
coord<-SpatialPoints(cbind(datatmp$lon,datatmp$lat),
                     proj4string=CRS(argv$proj4from))
coord.new<-spTransform(coord,CRS(argv$proj4to))
xy.new<-coordinates(coord.new)
x<-round(xy.new[,1],6)
y<-round(xy.new[,2],6)
rm(xy.new)
datatmp$lat<-y
datatmp$lon<-x
#
#------------------------------------------------------------------------------
# write the output file
ix<-which(!is.na(datatmp$lat) & !is.na(datatmp$lon) & 
          !is.na(datatmp$elev) & !is.na(datatmp$value) )
data<-data.frame(datatmp$lat[ix],
                 datatmp$lon[ix],
                 datatmp$elev[ix],
                 datatmp$value[ix])
names(data)<-c("lat","lon","elev","value")
dataopt<-dataopttmp[ix,,drop=F]
#
#------------------------------------------------------------------------------
# write the output file
varidx.out<-varidx
if (any(!is.na(argv$varname.opt))) 
  varidx.out<-c(varidx,varidx.opt[which(!is.na(varidx.opt))]) 
dataout<-array(data=NA,
               dim=c(length(data$lat),(length(varidx.out))))
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
    dataout[,s]<-round(data$elev,1)
  } else if (pos.s==4) {
    str[s]<-"value"
    dataout[,s]<-round(data$value,1)
  }
}
#str[s+1]<-argv$varname.prid
#dataout[,(s+1)]<-data$prid
#str[s+2]<-argv$varname.dqc
#dataout[,(s+2)]<-dqcflag
#str[s+3]<-argv$varname.sct
#dataout[,(s+3)]<-round(sctpog,2)
#str[s+4]<-argv$varname.rep
#dataout[,(s+4)]<-round(corep,5)
print(str)
dataout<-as.data.frame(dataout,stringsAsFactors=F)
names(dataout)<-str
write.table(file=argv$output,
            dataout,
            quote=F,
            col.names=str,
            row.names=F,
            dec=".",
            sep=";")
#write.table(file=argv$output,
#            data,
#            quote=F,
#            col.names=names(datatmp),
#            row.names=F,
#            dec=".",
#            sep=";")
#
#------------------------------------------------------------------------------
# Exit
print("Normal Exit")
quit(status=0)
