This directory contains physicochemical measurements, to be contextualized with eDNA data. The procedure how to connect ASV and metadata is described in script [DataLoad.R](../DataLoad.R). For this, place all files in a single R working directory of your choice.

Physicochemical data includes:
- [complete sample info](./sample_info.txt)
- [oceanographic data](./CTD.txt): temperature, salinity, depth, oxygen concentration + saturation, sensor-chlorophyll, pH, pCO2, Atlantic/Polar Water proportions
- [Water column stratification](./Strat.txt): mixed layer depth and related variables
- [Nutrient concentrations](./Nutrients.txt): nitrate + nitrite, nitrite, phosphate, silicate
- [Sea-ice concentration](./IceConc.txt)
- [Distance to sea-ice edge](./IceDist.txt)
- [satellite chlorophyll](./Chl_sat.txt)

In case samples from several timepoints were pooled, a "mean date" is calculated, and the corresponding environmental parameters averaged as well.


