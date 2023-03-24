This directory contains physicochemical measurements, to be contextualized with eDNA data. The procedure how to connect ASV and metadata is described in script [DataLoad.R](../DataLoad.R). For this, place all files in a single R working directory of your choice.

Physicochemical data includes:
- [complete sample info](./sample_info.txt)
- [oceanographic variables](./CTD.txt)
- [Water column stratification](./Strat.txt)
- [Nutrient concentrations](./Nutrients.txt)
- [Sea-ice concentration](./IceConc.txt)
- [Distance to sea-ice edge](./IceDist.txt)
- [satellite chlorophyll](./Chl_sat.txt)

In case samples from several timepoints were pooled, a "mean date" is calculated, and the corresponding environmental parameters averaged as well.


