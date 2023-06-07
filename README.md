# METR-5673

RadarCode.m and RainGaugeCode.m are Matlab scripts used to read and QC KHGX radar and Harris County Rain Gauge data between 25-29 August 2017 during Hurricane Harvey


RainGaugeCode.m can be run once you have downloaded the .xlsx files that hold the rain gauge data.

RadarCode.m cannot be run since the Z and ZDR files are not included (too large). I reccommend downloading all the gridded .nc files from 25-29 August 2017 at the lowest elevation angle for both Z and ZDR from KHGX if you are interested in running that script.


It is NOT needed to run the RainGaugeCode.m and RadarCode.m since I have provided the .mat files from each script after it runs.


Once you load RainGaugeData.mat and Z&ZDRData.mat, you are able to run QauntifyingPrecipCode.m. There are lots of comments so please go through them to understand it all. It takes a few minutes to run the full script < 3 min, so please be patient :)


If you have any questions, please message me.


Acknowledgements:
Harris County Flood Control District for providing the rain gauge data
NOAA Climate Weather Toolkit for providing the KHGX files
