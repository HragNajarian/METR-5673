%% Load external arrays

load("RainGaugedata.mat");
SiteLat = GaugeLocations.latitude;
SiteLon = GaugeLocations.longitude;
%% Preliminary
FRAME_M = 584;              % 5 minute interval to analyze
FRAME_H = 49;               % 60 minute interval to analyze
                            % Highest intensity occurs at 584/49
I_PTS = length(SiteLat);   % Number of interpolation points

RGDH = reshape(transpose(sum(reshape(RainGaugedata,12,[]))),[120 150]);

%% Plotting rain gauge data in 5 minute intervals

figure("Name","Rain Gauge Density Plot")
geodensityplot(SiteLat,SiteLon,RainGaugedata(FRAME_M,:))
geolimits([29.5 30],[-96 -94.5])
geobasemap streets
title("Rain Gauge Density Plot at Frame "+FRAME_M);
movegui('northwest');
figure("Name","Rain Gauge Bubble Plot")
geobubble(SiteLat,SiteLon,RainGaugedata(FRAME_M,:))
geolimits([29.5 30],[-96 -94.5])
geobasemap streets
title("Rain Gauge Density Plot at Frame "+FRAME_M);
movegui('southwest');

%% Plotting rain gauge data in 60 minute intervals

figure("Name","Rain Gauge Density Plot, Hourly")
geodensityplot(SiteLat,SiteLon,RGDH(FRAME_H,:))
geolimits([29.5 30],[-96 -94.5])
geobasemap streets
title("Rain Gauge Density Plot at Frame "+FRAME_H);
movegui('northwest');
figure("Name","Rain Gauge Bubble Plot, Hourly")
geobubble(SiteLat,SiteLon,RGDH(FRAME_H,:))
geolimits([29.5 30],[-96 -94.5])
geobasemap streets
title("Rain Gauge Density Plot at Frame "+FRAME_H);
movegui('southwest');

%% Interpolation
%  vq = griddata(x,y,v,xq,yq) is used to fit a surface defined by (x,y,v)
%                             at query points (xq, yq).


[xR,yR] = meshgrid(linspace(min(SiteLat),max(SiteLat),I_PTS), ...
    linspace(min(SiteLon),max(SiteLon),I_PTS));

iRGD1 = griddata(SiteLat,SiteLon,RGDH(FRAME_H,:),xR,yR);

RxR = reshape(xR,[1 I_PTS^2]);
RyR = reshape(yR,[1 I_PTS^2]);
RiD = reshape(iRGD1,[1 I_PTS^2]);

figure("Name","Interploated Rain Gauge Density Plot")
geodensityplot(RxR,RyR,RiD)
geolimits([29.5 30],[-96 -94.5])
geobasemap streets
title("Interpolated Rain Gauge Data at Frame "+FRAME_H)

%% Define Harris County borders

HarrisCountyP01 = [29.580138771715720 -95.42413999039637];
HarrisCountyP02 = [29.788625856004195 -95.82594752936723];
HarrisCountyP03 = [30.054674035807650 -95.92144314890938];
HarrisCountyP04 = [30.163577521244953 -95.96071237498661];
HarrisCountyP05 = [30.081711418469734 -95.82836297631081];
HarrisCountyP06 = [30.103219047502340 -95.66777899236332];
HarrisCountyP07 = [30.170369966756630 -95.54157555354173];
HarrisCountyP08 = [30.103939929561548 -95.47329317575064];
HarrisCountyP09 = [30.115211701348695 -95.42399755115744];
HarrisCountyP10 = [30.027813222004085 -95.29288114984554];
HarrisCountyP11 = [30.032386725851133 -95.26371447998915];
HarrisCountyP12 = [30.167148927317560 -95.09670830343960];
HarrisCountyP13 = [29.993335388312110 -95.03118639884012];
HarrisCountyP14 = [29.826490534313580 -94.90950942574400];
HarrisCountyP15 = [29.673725883432630 -94.93020197150076];
HarrisCountyP16 = [29.548806498530528 -95.01780158735795];
HarrisCountyP17 = [29.497675718901960 -95.16314627579786];
HarrisCountyP18 = [29.553978674619206 -95.25396597773938];
HarrisCountyP19 = [29.595149445465353 -95.27483841576833];

HarrisCountyL01Lat = [HarrisCountyP01(1) HarrisCountyP02(1)];
HarrisCountyL01Lon = [HarrisCountyP01(2) HarrisCountyP02(2)];
HarrisCountyL02Lat = [HarrisCountyP02(1) HarrisCountyP03(1)];
HarrisCountyL02Lon = [HarrisCountyP02(2) HarrisCountyP03(2)];
HarrisCountyL03Lat = [HarrisCountyP03(1) HarrisCountyP04(1)];
HarrisCountyL03Lon = [HarrisCountyP03(2) HarrisCountyP04(2)];
HarrisCountyL04Lat = [HarrisCountyP04(1) HarrisCountyP05(1)];
HarrisCountyL04Lon = [HarrisCountyP04(2) HarrisCountyP05(2)];
HarrisCountyL05Lat = [HarrisCountyP05(1) HarrisCountyP06(1)];
HarrisCountyL05Lon = [HarrisCountyP05(2) HarrisCountyP06(2)];
HarrisCountyL06Lat = [HarrisCountyP06(1) HarrisCountyP07(1)];
HarrisCountyL06Lon = [HarrisCountyP06(2) HarrisCountyP07(2)];
HarrisCountyL07Lat = [HarrisCountyP07(1) HarrisCountyP08(1)];
HarrisCountyL07Lon = [HarrisCountyP07(2) HarrisCountyP08(2)];
HarrisCountyL08Lat = [HarrisCountyP08(1) HarrisCountyP09(1)];
HarrisCountyL08Lon = [HarrisCountyP08(2) HarrisCountyP09(2)];
HarrisCountyL09Lat = [HarrisCountyP09(1) HarrisCountyP10(1)];
HarrisCountyL09Lon = [HarrisCountyP09(2) HarrisCountyP10(2)];
HarrisCountyL10Lat = [HarrisCountyP10(1) HarrisCountyP11(1)];
HarrisCountyL10Lon = [HarrisCountyP10(2) HarrisCountyP11(2)];
HarrisCountyL11Lat = [HarrisCountyP11(1) HarrisCountyP12(1)];
HarrisCountyL11Lon = [HarrisCountyP11(2) HarrisCountyP12(2)];
HarrisCountyL12Lat = [HarrisCountyP12(1) HarrisCountyP13(1)];
HarrisCountyL12Lon = [HarrisCountyP12(2) HarrisCountyP13(2)];
HarrisCountyL13Lat = [HarrisCountyP13(1) HarrisCountyP14(1)];
HarrisCountyL13Lon = [HarrisCountyP13(2) HarrisCountyP14(2)];
HarrisCountyL14Lat = [HarrisCountyP14(1) HarrisCountyP15(1)];
HarrisCountyL14Lon = [HarrisCountyP14(2) HarrisCountyP15(2)];
HarrisCountyL15Lat = [HarrisCountyP15(1) HarrisCountyP16(1)];
HarrisCountyL15Lon = [HarrisCountyP15(2) HarrisCountyP16(2)];
HarrisCountyL16Lat = [HarrisCountyP16(1) HarrisCountyP17(1)];
HarrisCountyL16Lon = [HarrisCountyP16(2) HarrisCountyP17(2)];
HarrisCountyL17Lat = [HarrisCountyP17(1) HarrisCountyP18(1)];
HarrisCountyL17Lon = [HarrisCountyP17(2) HarrisCountyP18(2)];
HarrisCountyL18Lat = [HarrisCountyP18(1) HarrisCountyP19(1)];
HarrisCountyL18Lon = [HarrisCountyP18(2) HarrisCountyP19(2)];
HarrisCountyL19Lat = [HarrisCountyP19(1) HarrisCountyP01(1)];
HarrisCountyL19Lon = [HarrisCountyP19(2) HarrisCountyP01(2)];

%% Plot the location of each rain gauge, of KHGX, and county border

figure("Name","Location Data")
geoplot(SiteLat,SiteLon,'*','Color','blue'); hold on;
geoplot(29.4719,-95.0788,'*','Color','red'); 
geoplot(HarrisCountyL01Lat,HarrisCountyL01Lon,'--','Color','black')
geoplot(HarrisCountyL02Lat,HarrisCountyL02Lon,'--','Color','black')
geoplot(HarrisCountyL03Lat,HarrisCountyL03Lon,'--','Color','black')
geoplot(HarrisCountyL04Lat,HarrisCountyL04Lon,'--','Color','black')
geoplot(HarrisCountyL05Lat,HarrisCountyL05Lon,'--','Color','black')
geoplot(HarrisCountyL06Lat,HarrisCountyL06Lon,'--','Color','black')
geoplot(HarrisCountyL07Lat,HarrisCountyL07Lon,'--','Color','black')
geoplot(HarrisCountyL08Lat,HarrisCountyL08Lon,'--','Color','black')
geoplot(HarrisCountyL09Lat,HarrisCountyL09Lon,'--','Color','black')
geoplot(HarrisCountyL10Lat,HarrisCountyL10Lon,'--','Color','black')
geoplot(HarrisCountyL11Lat,HarrisCountyL11Lon,'--','Color','black')
geoplot(HarrisCountyL12Lat,HarrisCountyL12Lon,'--','Color','black')
geoplot(HarrisCountyL13Lat,HarrisCountyL13Lon,'--','Color','black')
geoplot(HarrisCountyL14Lat,HarrisCountyL14Lon,'--','Color','black')
geoplot(HarrisCountyL15Lat,HarrisCountyL15Lon,'--','Color','black')
geoplot(HarrisCountyL16Lat,HarrisCountyL16Lon,'--','Color','black')
geoplot(HarrisCountyL17Lat,HarrisCountyL17Lon,'--','Color','black')
geoplot(HarrisCountyL18Lat,HarrisCountyL18Lon,'--','Color','black')
geoplot(HarrisCountyL19Lat,HarrisCountyL19Lon,'--','Color','black')
hold off;
geolimits([29.65 30],[-96.1 -94.6])
geobasemap streets
legend(["Rain Gauge","KHGX"]);
title('Locations of HCFWS Rain Gauges and KHGX')
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 7.25], ...
    'PaperUnits', 'Inches', 'PaperSize', [7.25, 7.25])
