* Replication Code for Koch data
use "C:\R\match\docs\papers\koch\genharvard.dta", clear

*Table 2 
* Democrat
reg pdcanid dviswom dvisman demcan1 dempty rideo dproj demft aware
* Democrat (voters only)
reg pdcanid dviswom dvisman demcan1 dempty rideo dproj demft aware if(voter==1)
* Republican
reg prcanid rviswom rvisman repcan1 goppty rideo rproj repft aware
* Republican (voters only)
reg prcanid rviswom rvisman repcan1 goppty rideo rproj repft aware if(voter==1)
