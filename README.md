# POSEIDON
(almost) Acronym of _PrecipitatiOn SpatIal Data quality CONtrol_

Available checks are (applied sequentially as in this list):

* Plausibility check

* Climatological check, predefined range for each month (optional)

* check for isolated dry observations

* check for isolated wet observations

* Buddy-check

* Spatial Consistency Test (SCT)

* check elevations against digital elevation model (optional)

* detect isolated observations

Possibility to have observation black-list and keep(-it-no-matter-what)-list.

Installation Instructions
-------------------------
Ensure the following R-libraries (and their dependencies) are installed:

   * argparser
   * sp
   * raster
   * rgdal
   * ncdf4 (optional, used only if additional geographical information are required)


Running the program
-------------------
To see program options, run:

```
   poseidon.R --help
```

run a test case with:

```
./poseidon.R test/TA_2017072112.txt test/dqc_2017072112.txt --spatconv --month.clim 7 --i.sct 3 --i.buddy 3 -v
```

Copyright and license
---------------------
Copyright (C) 2017 MET Norway. POSEIDON is licensed under [GPL
version 3](https://github.com/cristianlussana/POSEIDON/blob/master/LICENSE) or (at
your option) any later version.

Contact
-------
E-mail: cristianl@met.no
