%% Important variables to be familiar with

% Radar

    % ZData: Reflectivity Data (dBZ)
    % ZDRData: Differential Reflectivity (dBZ)
    % latZ and latZDR: Gridded latitude depending on Z and ZDR data (they
        % are different because of issues out of our control)
    % lonZ and lonZDR: Gridded latitude depending on Z and ZDR data (they
        % are different because of issues out of our control)
    % timeRadar: SECONDS since jan 0, 0000. Divide by 86400 before using
        % datetime
    % baseTime: Radar time but SECONDS since 1970-1-1
    % fileInfoZ/ZDR: Information on the .nc files
    
% Rain Gauge
    
    % RainGaugedata [mm/hr]: Rain rates from rain 150 gauges from Aug 25-29.
        % [5-minute intervals X number of rain gauges]
        % The order of rain gauges corresponds to GaugeLocations below.
        
    % GaugeLocations: Rain gauge locations(table) 
        % siteNum: Identifies the sites
        % latitude and longitude: Identifies the location of the respective sites
    % endTime: The end of the 5-minute collection of rain rate data.
    
       
%% Load .mat files

% WSR-88D KHGX data
load Z&ZDRData.mat

% Harris County Flood Control District
% Houston Rain Gauge Network data
load RainGaugeData.mat

%% Match up the lat and lon between gauges and Z/ZDR
set(groot, 'DefaultAxesFontName', 'Arial')
% matchCoord represents the [lon,lat] for each site (150) at each time step (1337)
    % for Z (1) and ZDR (2)
matchCoord = cell(150,2);
% Basically we need to find the closest pixel that matches closest to the
    % rain gauages lon and lat. To do this, we need to find the nearest
    % neighbor radar pixel for each radar site.

% matchCoord will provide the index that corresponds to the radars lon and 
    % lat that is the nearest neighbor for each site.
% Go to the end of this section for an example to make sure that the
    % matchCoord is doing its job correctly

for i = 1:length(matchCoord)
    
    coordZ = zeros(length(ZData),2);
    coordZDR = zeros(length(ZDRData),2);
    
    for j = 1:length(ZData)
        
        [~,coordZ(j,1)] = min(abs(GaugeLocations.longitude(i)-lonZ{j}));
        [~,coordZ(j,2)] = min(abs(GaugeLocations.latitude(i)-latZ{j}));
        
        [~,coordZDR(j,1)] = min(abs(GaugeLocations.longitude(i)-lonZDR{j}));
        [~,coordZDR(j,2)] = min(abs(GaugeLocations.latitude(i)-latZDR{j}));
        
    end
    
    matchCoord{i,1} = coordZ;
    matchCoord{i,2} = coordZDR;
    
end
clear coordZ coordZDR i j

% i.e. if you want to know the radar lon and lat that best matched the 
    % first site's lon and lat, for the first time step run: 
    % i and j picks random times and site to compare respectively
    i = randi(length(ZData));
    j = randi(length(matchCoord));
   sprintf('Radar lon,lat: %0.4f, %0.4f',lonZ{i}(matchCoord{j,1}(i,1)),latZ{i}(matchCoord{j,1}(i,2)))
   sprintf('Rain Gauge lon,lat: %0.4f, %0.4f',GaugeLocations.longitude(j),GaugeLocations.latitude(j))
clear i j

%% Match up the times between gauges and radar

% Times between Z and ZDR are the same, and are the same all across the
    % times
% This function below is required since the radar scans are over 5 minutes
    % long which means the rain gauge data is higher frequency than radar.
    % Since they don't match in temporal resolution, we need to find the
    % times that best match.
    
% matchTime will hold the indicies for RainGauges that best match times of
    % radar data
matchTime = zeros(length(ZData),1);

for i = 1:length(ZData)    % loop over each site
    
    [~,matchTime(i)] = min(abs(timeRadar(i)-endTime));
    
end

%% Find the site distance from radar

radarDistance = zeros(length(matchCoord),1);

for i = 1:length(radarDistance)
    
    dx = abs(lonZ{1}(centerCoord(1,1)) - GaugeLocations.longitude(i));
    dy = abs(latZ{1}(centerCoord(1,2)) - GaugeLocations.latitude(i));
    
    % a^2 + b^2 = c^2 to get radius from radar
    radarDistance(i) = sqrt((dx*111)^2 +(dy*111)^2);    % [km]
    
end

% Find the time it takes for the retreival above each site to be detected
% Assmume 8 m/s fall rate
fallRate = 8;   
t = (tand(0.46)*radarDistance*1000)/fallRate;   % [seconds]

% Add radarDistance to GaugeLocations
T = table(radarDistance,t,'VariableNames',{'radarDistance','time2ground'});
GaugeLocations = [GaugeLocations,T];

clear T t fallRate dx dy i

%% Match the times of rain being detected from the received radar signal

matchTimeAdj = repmat(matchTime,1,150);

for i = 1:length(matchCoord)
    
    for j = 1:length(ZData)
        
        [~,matchTimeAdj(j,i)] = min(abs((timeRadar(j)+GaugeLocations.time2ground(i))-endTime));
    
    end
    
end

clear i

%% Find a relationship between Z and ZDR and rain rate

% Reflectivity over the site (no smoothing/just 1 pixel)
matchZ = zeros(size(ZData,1),size(RainGaugedata,2));
matchZDR = zeros(size(ZDRData,1),size(RainGaugedata,2),1);

% This loop's purpose is to collect the Z and ZDR values at each site over
    % time
for i = 1:size(RainGaugedata,2) % Loop through all the sites
    
    for j = 1:length(ZData) % Loop through the radar time steps
        
        lon = matchCoord{i,1}(j,1);
        lat = matchCoord{i,1}(j,2);
        matchZ(j,i) = ZData{j}(lon,lat);

        lon = matchCoord{i,2}(j,1);
        lat = matchCoord{i,2}(j,2);
        matchZDR(j,i) = ZDRData{j}(lon,lat);
        
    end
    
end

%%
% Reflectivity over the site (Smoothing over 3x3 grid centered over rain gauge)
matchZSmooth = zeros(size(ZData,1),size(RainGaugedata,2));
matchZDRSmooth = zeros(size(ZDRData,1),size(RainGaugedata,2),1);

% This loop's purpose is to collect the Z and ZDR values at each site over
    % time
for i = 1:size(RainGaugedata,2) % Loop through all the sites
    
    for j = 1:length(ZData) % Loop through the radar time steps
        
        lon = matchCoord{i,1}(j,1);
        lat = matchCoord{i,1}(j,2);
        matchZSmooth(j,i) = mean(ZData{j}(lon-1:lon+1,lat-1:lat+1),'all','omitnan');

        lon = matchCoord{i,2}(j,1);
        lat = matchCoord{i,2}(j,2);
        matchZDRSmooth(j,i) = mean(ZDRData{j}(lon-1:lon+1,lat-1:lat+1),'all','omitnan');
        
    end
    
end

clear i j lon lat

%% Compare and contrast between the different methods Z

f = figure('Position',[-210 1060 1047 842]);

% Non-time adjusted + non-smooth
subplot(2,2,1)

x = matchZ;
y = RainGaugedata(matchTime,:); % only look at rain gauge data that match the times

% This portion will take rain rates and average them over their respective
    % reflectivity bin
xavg = -20:5:70;   % min Z is -17.5, max is 65.5 over the whole time
yavg = zeros(length(xavg),1);
n = zeros(length(xavg),4);    % This will provide the number of samples

for i = 1:(length(xavg)-1)
    
    % The averages the rain rate values within a reflectivity bin
    yavg(i) = mean(y(x<xavg(i+1) & x>=xavg(i)),'all','omitnan');
    
    % Counts the number of samples to check the potential of outliers
    n(i,1) = sum((x<xavg(i+1) & x>=xavg(i)),'all');
    
end
yavg(isnan(yavg)) = 0;

% Significance/Variance testing through the Monte Carlo method (itter=1000, w/ repetition)
itter = 1000;
bound = .10/2;  % 90% significance
yavgSig = zeros(length(xavg),2);
temp = zeros(itter,2);

for i = 1:(length(xavg)-1)
    
    % The averages the rain rate values within a reflectivity bin
    a = y(x<xavg(i+1) & x>=xavg(i));
    if isempty(a)
        continue
    else
        for j = 1:itter

            randIndex = randi(length(a),[length(a) 1]);
            temp(j) = mean(a(randIndex),'all','omitnan');
        end
    end
    
    temp = sort(temp);
    yavgSig(i,:) =  [temp(itter*bound) temp(itter*(1-bound))];
    
end

bar(xavg,yavg,.75)
hold on
% Plot the confidence intervals at each bin
er = errorbar(xavg(n(:,1)~=0),yavg(n(:,1)~=0),yavgSig((n(:,1)~=0),1)-yavg((n(:,1)~=0)),yavgSig((n(:,1)~=0),2)-yavg((n(:,1)~=0)));
er.Color = [0,0,0];
er.LineWidth = 1.5;
er.LineStyle = 'none';

% Plot the exponential fit
xtemp = xavg(n(:,1)~=0);
ytemp = yavg(n(:,1)~=0);
expfit = fit(xtemp',ytemp,'exp1');
p = plot(xavg,expfit(xavg));
p.LineWidth = 3;
p.Color = '#D95319';
temp = sqrt((1/length(ytemp))*sum(abs(ytemp-expfit(xtemp)).^2));

a = sprintf('Total RMSE: %1.4f',mean(temp)); % this is the RMSE from ALL the observations
text(20,40,a,'fontsize',14,'HorizontalAlignment','left','VerticalAlignment','bottom');

legend(p,sprintf('%1.4f*exp(%1.4fx)',expfit.a,expfit.b),'Location','northeast')

ax = gca;
ax.FontSize = 16;
ax.Box = 'on';
ax.XGrid = 'on';
ax.YGrid = 'on';
xticks(-20:10:70)
xlabel('Reflectivity [dBZ]')
ylabel('Rain Gauge Rain Rate [mm/hr]')
ylim([0,50])
xlim([-25,70])
title('Rain Rate and Reflectivity Relationship')
subtitle('Non-Time Adjusted, Non-Smooth')


% Time-adjusted + non-smooth
subplot(2,2,2)
x = matchZ;
y = RainGaugedata(matchTimeAdj); % only look at rain gauge data that match the times

% This portion will take rain rates and average them over their respective
    % reflectivity bin
    
for i = 1:(length(xavg)-1)
    
    % The averages the rain rate values within a reflectivity bin
    yavg(i) = mean(y(x<xavg(i+1) & x>=xavg(i)),'all','omitnan');
    n(i,2) = sum((x<xavg(i+1) & x>=xavg(i)),'all');
    
end
yavg(isnan(yavg)) = 0;

% Significance/Variance testing (itter=1000, w/ repeated)
yavgSig = zeros(length(xavg),2);
temp = zeros(itter,2);

for i = 1:(length(xavg)-1)
    
    % The averages the rain rate values within a reflectivity bin
    a = y(x<xavg(i+1) & x>=xavg(i));
    if isempty(a)
        continue
    else
        for j = 1:itter

            randIndex = randi(length(a),[length(a) 1]);
            temp(j) = mean(a(randIndex),'all','omitnan'); 
        end
    end
    
    temp = sort(temp);
    yavgSig(i,:) =  [temp(itter*bound) temp(itter*(1-bound))];
    
end

bar(xavg,yavg,.75)
hold on

er = errorbar(xavg(n(:,2)~=0),yavg(n(:,2)~=0),yavgSig((n(:,2)~=0),1)-yavg((n(:,2)~=0)),yavgSig((n(:,2)~=0),2)-yavg((n(:,2)~=0)));
er.Color = [0,0,0];
er.LineWidth = 2;
er.LineStyle = 'none';

% Plot the exponential fit
xtemp = xavg(n(:,2)~=0);
ytemp = yavg(n(:,2)~=0);
expfit = fit(xtemp',ytemp,'exp1');
p = plot(xavg,expfit(xavg));
p.LineWidth = 3;
p.Color = '#D95319';
temp = sqrt((1/length(ytemp))*sum(abs(ytemp-expfit(xtemp)).^2));

a = sprintf('Total RMSE: %1.4f',mean(temp)); % this is the RMSE from ALL the observations
text(-23,40,a,'fontsize',14,'HorizontalAlignment','left','VerticalAlignment','bottom');

legend(p,sprintf('%1.4f*exp(%1.4fx)',expfit.a,expfit.b),'Location','northwest')


ax = gca;
ax.FontSize = 16;
ax.Box = 'on';
ax.XGrid = 'on';
ax.YGrid = 'on';
xticks(-20:10:70)
xlabel('Reflectivity [dBZ]')
ylabel('Rain Gauge Rain Rate [mm/hr]')
ylim([0,50])
xlim([-25,70])
title('Rain Rate and Reflectivity Relationship')
subtitle('Time Adjusted, Non-Smooth')


% Non-time adjusted + Smoothed
subplot(2,2,3)

x = matchZSmooth;
y = RainGaugedata(matchTime,:); % only look at rain gauge data that match the times

% This portion will take rain rates and average them over their respective
    % reflectivity bin

for i = 1:(length(xavg)-1)
    
    % The averages the rain rate values within a reflectivity bin
    yavg(i) = mean(y(x<xavg(i+1) & x>=xavg(i)),'all','omitnan');
    n(i,3) = sum((x<xavg(i+1) & x>=xavg(i)),'all');
    
end
yavg(isnan(yavg)) = 0;
% Significance/Variance testing (itter=1000, w/ repeated)
yavgSig = zeros(length(xavg),2);
temp = zeros(itter,2);

for i = 1:(length(xavg)-1)
    
    % The averages the rain rate values within a reflectivity bin
    a = y(x<xavg(i+1) & x>=xavg(i));
    
    if isempty(a)
        continue
    else
        for j = 1:itter

            randIndex = randi(length(a),[length(a) 1]);
            temp(j) = mean(a(randIndex),'all','omitnan');

        end
    end
    
    temp = sort(temp);
    yavgSig(i,:) =  [temp(itter*bound) temp(itter*(1-bound))];
    
end

bar(xavg,yavg,.75)
hold on
er = errorbar(xavg(n(:,3)~=0),yavg(n(:,3)~=0),yavgSig((n(:,3)~=0),1)-yavg((n(:,3)~=0)),yavgSig((n(:,3)~=0),2)-yavg((n(:,3)~=0)));
er.Color = [0,0,0];
er.LineWidth = 2;
er.LineStyle = 'none';

% Plot the exponential fit
xtemp = xavg(n(:,3)~=0);
ytemp = yavg(n(:,3)~=0);
expfit = fit(xtemp',ytemp,'exp1');
p = plot(xavg,expfit(xavg));
p.LineWidth = 3;
p.Color = '#D95319';
temp = sqrt((1/length(ytemp))*sum(abs(ytemp-expfit(xtemp)).^2));

a = sprintf('Total RMSE: %1.4f',mean(temp)); % this is the RMSE from ALL the observations
text(20,40,a,'fontsize',14,'HorizontalAlignment','left','VerticalAlignment','bottom');

legend(p,sprintf('%1.4f*exp(%1.4fx)',expfit.a,expfit.b),'Location','northeast')

ax = gca;
ax.FontSize = 16;
ax.Box = 'on';
ax.XGrid = 'on';
ax.YGrid = 'on';
xticks(-20:10:70)
xlabel('Reflectivity [dBZ]')
ylabel('Rain Gauge Rain Rate [mm/hr]')
ylim([0,50])
xlim([-25,70])
title('Rain Rate and Reflectivity Relationship')
subtitle('Non-Time Adjusted, Smoothed')


% Time adjusted + Smoothed
subplot(2,2,4)

x = matchZSmooth;
y = RainGaugedata(matchTimeAdj); % only look at rain gauge data that match the times

% This portion will take rain rates and average them over their respective
    % reflectivity bin

for i = 1:(length(xavg)-1)
    
    % The averages the rain rate values within a reflectivity bin
    yavg(i) = mean(y(x<xavg(i+1) & x>=xavg(i)),'all','omitnan');
    n(i,4) = sum((x<xavg(i+1) & x>=xavg(i)),'all');
    
end
yavg(isnan(yavg)) = 0;

% Significance/Variance testing (itter=1000, w/ repeated)
yavgSig = zeros(length(xavg),2);
temp = zeros(itter,2);

for i = 1:(length(xavg)-1)
    
    % The averages the rain rate values within a reflectivity bin
    a = y(x<xavg(i+1) & x>=xavg(i));
    if isempty(a)
        continue
    else
        for j = 1:itter

            randIndex = randi(length(a),[length(a) 1]);
            temp(j) = mean(a(randIndex),'all','omitnan');
        end
    
    end
    
    temp = sort(temp);
    yavgSig(i,:) =  [temp(itter*bound) temp(itter*(1-bound))];
    
end

bar(xavg,yavg,.75)
hold on
er = errorbar(xavg(n(:,4)>0),yavg(n(:,4)>0),yavgSig((n(:,4)>0),1)-yavg((n(:,4)>0)),yavgSig((n(:,4)>0),2)-yavg((n(:,4)>0)));
er.Color = [0,0,0];
er.LineWidth = 2;
er.LineStyle = 'none';

% Plot the exponential fit
xtemp = xavg(n(:,4)>0);
ytemp = yavg(n(:,4)>0);
expfit = fit(xtemp',ytemp,'exp1');
p = plot(xavg,expfit(xavg));
p.LineWidth = 3;
p.Color = '#D95319';
temp = sqrt((1/length(ytemp))*sum(abs(ytemp-expfit(xtemp)).^2));

a = sprintf('Total RMSE: %1.4f',mean(temp)); % this is the RMSE from ALL the observations
text(-23,40,a,'fontsize',14,'HorizontalAlignment','left','VerticalAlignment','bottom');

legend(p,sprintf('%1.4f*exp(%1.4fx)',expfit.a,expfit.b),'Location','northwest')

ax = gca;
ax.FontSize = 16;
ax.Box = 'on';
ax.XGrid = 'on';
ax.YGrid = 'on';
xticks(-20:10:70)
xlabel('Reflectivity [dBZ]')
ylabel('Rain Gauge Rain Rate [mm/hr]')
ylim([0,50])
xlim([-25,70])
title('Rain Rate and Reflectivity Relationship')
subtitle('Time Adjusted, Smoothed')

% exportgraphics(f,'Rain_Rate-Reflectivity-Methods.pdf','ContentType','vector')
% clear a i j x y f ax temp itter er ans xavg yavg xavgSig yavgSig bound

%% Model testing with the best fit Z-R relationship 

% Time adjusted + Smoothed

% This portion will take rain rates and average them over their respective
    % reflectivity bin
xavg = -20:5:70;   % min Z is -17.5, max is 65.5 over the whole time
yavg = zeros(length(xavg),1);
n = zeros(length(xavg),4);    % This will provide the number of samples

% Create a random vector with a fraction of the rain gauge points to train the model
% This method assumes that each rain gauge provided equal amounts of data
    % for the model, which is a good assumption since a hurricane is a large
    % scale event
frac = 1;
ri = randperm(size(matchZSmooth,2),round(size(matchZSmooth,2)*frac));
rinot = 1:150;

% This look isolates the sites that are not in the training set. These
    % sites will be used to test the skill of the trained Z-R exp
    % relationship 
varb = 0;
temp = sort(ri);
for i = 1:150
    if rinot(i) == temp(i-varb)
        rinot(i) = 0;
%         varb = varb + 1;
    else
        varb = varb + 1;
    end
end
clear varb temp i

x = matchZSmooth(:,ri);
y = RainGaugedata(matchTimeAdj(:,ri)); % only look at rain gauge data that match the times

% This portion will take rain rates and average them over their respective
    % reflectivity bin

for i = 1:(length(xavg)-1)
    
    % The averages the rain rate values within a reflectivity bin
    yavg(i) = mean(y(x<xavg(i+1) & x>=xavg(i)),'all','omitnan');
    n(i,4) = sum((x<xavg(i+1) & x>=xavg(i)),'all');
    
end
yavg(isnan(yavg)) = 0;

% Significance/Variance testing (itter=1000, w/ repeated)
itter = 1000;
bound = .10/2;  % 90% significance
yavgSig = zeros(length(xavg),2);
temp = zeros(itter,2);

for i = 1:(length(xavg)-1)
    
    % The averages the rain rate values within a reflectivity bin
    a = y(x<xavg(i+1) & x>=xavg(i));
    if isempty(a)
        continue
    else
        for j = 1:itter

            randIndex = randi(length(a),[length(a) 1]);
            temp(j) = mean(a(randIndex),'all','omitnan');
        end
    
    end
    
    temp = sort(temp);
    yavgSig(i,:) =  [temp(itter*bound) temp(itter*(1-bound))];
    
end

f = figure;
bar(xavg,yavg,.75)
hold on
er = errorbar(xavg(n(:,4)>0),yavg(n(:,4)>0),yavgSig((n(:,4)>0),1)-yavg((n(:,4)>0)),yavgSig((n(:,4)>0),2)-yavg((n(:,4)>0)));
er.Color = [0,0,0];
er.LineWidth = 2;
er.LineStyle = 'none';


% % Significane testing the exp fit
% expfitSig = zeros(itter,2); % This will hold the coefficients (a*exp(b*x))
% xavg2 = -20:5:70;   % min Z is -17.5, max is 65.5 over the whole time
% yavg2 = zeros(length(xavg2),1);
% 
% ri = randperm(size(matchZSmooth,2),round(size(matchZSmooth,2)*frac));
% x = matchZSmooth(:,ri);
% y = RainGaugedata(matchTimeAdj(:,ri)); % only look at rain gauge data that match the times
% 
% 
% for i = 1:itter
%         
%     for j = 1:(length(xavg2)-1)
% 
%         % The averages the rain rate values within a reflectivity bin
%         yavg2(j) = mean(y(x<xavg2(j+1) & x>=xavg2(j)),'all','omitnan');
%         n(j,4) = sum((x<xavg2(j+1) & x>=xavg2(j)),'all');
% 
%     end
%     xtemp = xavg(n(:,4)>10);
%     ytemp = yavg2(n(:,4)>10);
%     ztemp = fit(xtemp',ytemp,'exp1');
%     expfitSig(i,:) = coeffvalues(ztemp);    % Gives you coefficients a & b
%     
% end
% yavg(isnan(yavg)) = 0;


% Plot the exponential fit
    %  n(:,4)>10 makes it so the bins with 0 values don't bias the exp fit 
xtemp = xavg(n(:,4)>0);
ytemp = yavg(n(:,4)>0);
expfit = fit(xtemp',ytemp,'exp1');
expcoefaSmooth = expfit.a;
expcoefbSmooth = expfit.b;
p = plot(xavg,expfit(xavg));

a = expfit(xtemp);
% Root-mean Square Error
expRMSE = sqrt((1/length(ytemp))*sum(abs(ytemp-a).^2));

b = sprintf('Total RMSE: %1.4f',mean(expRMSE)); % this is the RMSE from ALL the observations
text(-23,40,b,'fontsize',14,'HorizontalAlignment','left','VerticalAlignment','bottom');


% % Plotting the 95% significance (don't need to include)
% temp = confint(expfit);
% a = temp(:,1); b = temp(:,2);
% plow = plot(xavg,(a(1)*exp(b(1)*xavg)),':');
% phigh = plot(xavg,(a(2)*exp(b(2)*xavg)),'--');
p.LineWidth = 3;
p.Color = '#D95319';
legend(p,sprintf('%1.4f*exp(%1.4fx)',expfit.a,expfit.b),'Location','northwest')

ax = gca;
ax.FontSize = 16;
ax.Box = 'on';
ax.XGrid = 'on';
ax.YGrid = 'on';
xticks(-20:10:70)
xlabel('Reflectivity [dBZ]')
ylabel('Rain Gauge Rain Rate [mm/hr]')
ylim([0,50])
xlim([-25,70])
title('Rain Rate and Reflectivity Relationship')
subtitle('Time Adjusted, Smoothed')

% clear a c i j x y f p ax temp itter er ans xavg yavg xavgSig yavgSig bound 


% Model testing with the best fit Z-R relationship 

% This is a little complex, but essentially, we are going to select a
    % fraction of the sites to "train" the exp function, and then we will
    % apply that function to the fraction of sites that weren't trained on
    % to see how well they predict/quantify precipitation over those sites.
    % We also want to be able to have a diverse sample of sites that were and
    % weren't trained on so this will be looped over 'itterk' times to be
    % more robust.
% We will be utilizing the time adjusted + smoothed data set as that had
    % the best fit

itterk = 100;
expcoefa = zeros(itterk,1);
expcoefb = zeros(itterk,1);
expRMSE = zeros(itterk,1);

% This will keep track of the sites that were not used so we can
    % calculate how well the exp function preformed on them
rinot = repmat(1:150,itterk,1);
for k = 1:itterk
    
    % This portion will take rain rates and average them over their respective
        % reflectivity bin
    xavg = -20:5:70;   % min Z is -17.5, max is 65.5 over the whole time
    yavg = zeros(length(xavg),1);
%     n = zeros(length(xavg),1);    % This variable provides the number of samples

    % Below we will crete a random vector with a fraction of the rain gauge 
        % points to train the model.
    % This method assumes that each rain gauge provided equal amounts of data
        % for the model, which is a good assumption since a hurricane is a large
        % scale event
    frac = .75;     % The fraction of total sites we will use
    % This creates a random non-repeated vector
    ri = randperm(size(matchZSmooth,2),round(size(matchZSmooth,2)*frac));
   
    if frac~=1
        % This look isolates the sites that are not in the training set. These
            % sites will be used to test the skill of the trained Z-R exp
            % relationship 
        varb = 0;
        temp = sort(ri);
        for i = 1:max(temp)
    %         if rinot(i-1)==150
    %             continue
            if rinot(k,i) == temp(i-varb)
                rinot(k,i) = 0;
        %         varb = varb + 1;
            else
                varb = varb + 1;
            end
        end
        clear varb temp i
    end

    x = matchZSmooth(:,ri);
    y = RainGaugedata(matchTimeAdj(:,ri)); % only look at rain gauge data that match the times

    % This portion will take rain rates and average them over their respective
        % reflectivity bin

    for i = 1:(length(xavg)-1)

        % The averages the rain rate values within a reflectivity bin
        yavg(i) = mean(y(x<xavg(i+1) & x>=xavg(i)),'all','omitnan');
        n(i,4) = sum((x<xavg(i+1) & x>=xavg(i)),'all');

    end
    yavg(isnan(yavg)) = 0;

    
    xtemp = xavg(n(:,4)>0);
    ytemp = yavg(n(:,4)>0);
    expfit = fit(xtemp',ytemp,'exp1');
    % Save the coefficients of the 
    expcoefa(k) = expfit.a;
    expcoefb(k) = expfit.b;
    a = expfit(xtemp);
    % Root-mean Square Error
    expRMSE(k) = sqrt((1/length(ytemp))*sum(abs(ytemp-a).^2));
    
end

testRMSE = zeros(itterk,1);
temp = expcoefa.*exp(expcoefb*xavg);
for k = 1:itterk
    x = matchZSmooth(:,rinot(k,:)~=0);
    y = RainGaugedata(matchTimeAdj(:,rinot(k,:)~=0)); % only look at rain gauge data that match the times

    % This portion will take rain rates and average them over their respective
        % reflectivity bin

    for i = 1:(length(xavg)-1)

        % The averages the rain rate values within a reflectivity bin
        yavg(i) = mean(y(x<xavg(i+1) & x>=xavg(i)),'all','omitnan');
        n(i,4) = sum((x<xavg(i+1) & x>=xavg(i)),'all');

    end
    yavg(isnan(yavg)) = 0;

    xtemp = xavg(n(:,4)>0);
    ytemp = yavg(n(:,4)>0);
    a = temp(k,n(:,4)>0);

    testRMSE(k) = sqrt((1/length(ytemp))*sum(abs(ytemp-a').^2));
end

% b = sprintf('Total RMSE: %1.4f',mean(expRMSE)); % this is the RMSE from ALL the observations
a = sprintf('Test RMSE: %1.4f',mean(testRMSE)); % This is the RMSE when using the trained exp function on a frac of testing data sites

text(-23,42.5,a,'fontsize',14,'HorizontalAlignment','left','VerticalAlignment','bottom');
% text(-23,40,b,'fontsize',14,'HorizontalAlignment','left','VerticalAlignment','bottom');

% Save the figure
% exportgraphics(f,'Time-Adjusted_Smoothed.pdf','ContentType','vector')

%%
% %% Model testing with the best fit Z-R relationship 
% 
% % Non-Time adjusted + Non-Smoothed
% 
% % This portion will take rain rates and average them over their respective
%     % reflectivity bin
% xavg = -20:5:70;   % min Z is -17.5, max is 65.5 over the whole time
% yavg = zeros(length(xavg),1);
% n = zeros(length(xavg),4);    % This will provide the number of samples
% 
% % Create a random vector with a fraction of the rain gauge points to train the model
% % This method assumes that each rain gauge provided equal amounts of data
%     % for the model, which is a good assumption since a hurricane is a large
%     % scale event
% frac = 1;
% ri = randperm(size(matchZ,2),round(size(matchZ,2)*frac));
% rinot = 1:150;
% 
% % This look isolates the sites that are not in the training set. These
%     % sites will be used to test the skill of the trained Z-R exp
%     % relationship 
% varb = 0;
% temp = sort(ri);
% for i = 1:150
%     if rinot(i) == temp(i-varb)
%         rinot(i) = 0;
%     else
%         varb = varb + 1;
%     end
% end
% clear varb temp i
% 
% x = matchZ(:,ri);
% y = RainGaugedata(matchTime,ri); % only look at rain gauge data that match the times
% 
% % This portion will take rain rates and average them over their respective
%     % reflectivity bin
% 
% for i = 1:(length(xavg)-1)
%     
%     % The averages the rain rate values within a reflectivity bin
%     yavg(i) = mean(y(x<xavg(i+1) & x>=xavg(i)),'all','omitnan');
%     n(i,4) = sum((x<xavg(i+1) & x>=xavg(i)),'all');
%     
% end
% yavg(isnan(yavg)) = 0;
% 
% % Significance/Variance testing (itter=1000, w/ repeated)
% itter = 1000;
% bound = .10/2;  % 90% significance
% yavgSig = zeros(length(xavg),2);
% temp = zeros(itter,2);
% 
% for i = 1:(length(xavg)-1)
%     
%     % The averages the rain rate values within a reflectivity bin
%     a = y(x<xavg(i+1) & x>=xavg(i));
%     if isempty(a)
%         continue
%     else
%         for j = 1:itter
% 
%             randIndex = randi(length(a),[length(a) 1]);
%             temp(j) = mean(a(randIndex),'all','omitnan');
%         end
%     
%     end
%     
%     temp = sort(temp);
%     yavgSig(i,:) =  [temp(itter*bound) temp(itter*(1-bound))];
%     
% end
% 
% figure;
% bar(xavg,yavg,.75)
% hold on
% er = errorbar(xavg(n(:,4)>0),yavg(n(:,4)>0),yavgSig((n(:,4)>0),1)-yavg((n(:,4)>0)),yavgSig((n(:,4)>0),2)-yavg((n(:,4)>0)));
% er.Color = [0,0,0];
% er.LineWidth = 2;
% er.LineStyle = 'none';
% 
% 
% % % Significane testing the exp fit
% % expfitSig = zeros(itter,2); % This will hold the coefficients (a*exp(b*x))
% % xavg2 = -20:5:70;   % min Z is -17.5, max is 65.5 over the whole time
% % yavg2 = zeros(length(xavg2),1);
% % 
% % ri = randperm(size(matchZSmooth,2),round(size(matchZSmooth,2)*frac));
% % x = matchZSmooth(:,ri);
% % y = RainGaugedata(matchTimeAdj(:,ri)); % only look at rain gauge data that match the times
% % 
% % 
% % for i = 1:itter
% %         
% %     for j = 1:(length(xavg2)-1)
% % 
% %         % The averages the rain rate values within a reflectivity bin
% %         yavg2(j) = mean(y(x<xavg2(j+1) & x>=xavg2(j)),'all','omitnan');
% %         n(j,4) = sum((x<xavg2(j+1) & x>=xavg2(j)),'all');
% % 
% %     end
% %     xtemp = xavg(n(:,4)>0);
% %     ytemp = yavg2(n(:,4)>0);
% %     ztemp = fit(xtemp',ytemp,'exp1');
% %     expfitSig(i,:) = coeffvalues(ztemp);    % Gives you coefficients a & b
% %     
% % end
% % yavg(isnan(yavg)) = 0;
% 
% 
% % Plot the exponential fit
%     %  n(:,4)>0 makes it so the bins with 0 values don't bias the exp fit 
% xtemp = xavg(n(:,4)>0);
% ytemp = yavg(n(:,4)>0);
% expfit = fit(xtemp',ytemp,'exp1');
% p = plot(xavg,expfit(xavg));
% 
% % % Plotting the 95% significance (don't need to include)
% % temp = confint(expfit);
% % a = temp(:,1); b = temp(:,2);
% % plow = plot(xavg,(a(1)*exp(b(1)*xavg)),':');
% % phigh = plot(xavg,(a(2)*exp(b(2)*xavg)),'--');
% p.LineWidth = 3;
% p.Color = '#D95319';
% legend(p,sprintf('%1.4f*exp(%1.4fx)',expfit.a,expfit.b),'Location','northwest')
% 
% ax = gca;
% ax.FontSize = 16;
% ax.Box = 'on';
% ax.XGrid = 'on';
% ax.YGrid = 'on';
% xticks(-20:10:70)
% xlabel('Reflectivity [dBZ]')
% ylabel('Rain Gauge Rain Rate [mm/hr]')
% ylim([0,50])
% xlim([-25,70])
% title('Rain Rate and Reflectivity Relationship')
% subtitle('Non-Time Adjusted, Non-Smoothed')
% 
% % clear a c i j x y f p ax temp itter er ans xavg yavg xavgSig yavgSig bound 
% 
% 
% % Model testing with the best fit Z-R relationship 
% 
% % This is a little complex, but essentially, we are going to select a
%     % fraction of the sites to "train" the exp function, and then we will
%     % apply that function to the fraction of sites that weren't trained on
%     % to see how well they predict/quantify precipitation over those sites.
%     % We also want to be able to have a diverse sample of sites that were and
%     % weren't trained on so this will be looped over 'itterk' times to be
%     % more robust.
% % We will be utilizing the time adjusted + smoothed data set as that had
%     % the best fit
% 
% itterk = 100;
% expcoefa = zeros(itterk,1);
% expcoefb = zeros(itterk,1);
% expRMSE = zeros(itterk,1);
% 
% % This will keep track of the sites that were not used so we can
%     % calculate how well the exp function preformed on them
% rinot = repmat(1:150,itterk,1);
% for k = 1:itterk
%     
%     % This portion will take rain rates and average them over their respective
%         % reflectivity bin
%     xavg = -20:5:70;   % min Z is -17.5, max is 65.5 over the whole time
%     yavg = zeros(length(xavg),1);
% %     n = zeros(length(xavg),1);    % This variable provides the number of samples
% 
%     % Below we will crete a random vector with a fraction of the rain gauge 
%         % points to train the model.
%     % This method assumes that each rain gauge provided equal amounts of data
%         % for the model, which is a good assumption since a hurricane is a large
%         % scale event
%     frac = .75;     % The fraction of total sites we will use
%     % This creates a random non-repeated vector
%     ri = randperm(size(matchZ,2),round(size(matchZ,2)*frac));
%    
%     if frac~=1
%         % This look isolates the sites that are not in the training set. These
%             % sites will be used to test the skill of the trained Z-R exp
%             % relationship 
%         varb = 0;
%         temp = sort(ri);
%         for i = 1:max(temp)
%     %         if rinot(i-1)==150
%     %             continue
%             if rinot(k,i) == temp(i-varb)
%                 rinot(k,i) = 0;
%         %         varb = varb + 1;
%             else
%                 varb = varb + 1;
%             end
%         end
%         clear varb temp i
%     end
% 
%     x = matchZ(:,ri);
%     y = RainGaugedata(matchTime,:); % only look at rain gauge data that match the times
% 
%     % This portion will take rain rates and average them over their respective
%         % reflectivity bin
% 
%     for i = 1:(length(xavg)-1)
% 
%         % The averages the rain rate values within a reflectivity bin
%         yavg(i) = mean(y(x<xavg(i+1) & x>=xavg(i)),'all','omitnan');
%         n(i,4) = sum((x<xavg(i+1) & x>=xavg(i)),'all');
% 
%     end
%     yavg(isnan(yavg)) = 0;
% 
%     
%     xtemp = xavg(n(:,4)>0);
%     ytemp = yavg(n(:,4)>0);
%     expfit = fit(xtemp',ytemp,'exp1');
%     % Save the coefficients of the 
%     expcoefa(k) = expfit.a;
%     expcoefb(k) = expfit.b;
%     a = expfit(xtemp);
%     % Root-mean Square Error
%     expRMSE(k) = sqrt((1/length(ytemp))*sum(abs(ytemp-a).^2));
%     
% end
% 
% testRMSE = zeros(itterk,1);
% temp = expcoefa.*exp(expcoefb*xavg);
% for k = 1:itterk
%     x = matchZ(:);
%     y = RainGaugedata(matchTime,:); % only look at rain gauge data that match the times
% 
%     % This portion will take rain rates and average them over their respective
%         % reflectivity bin
% 
%     for i = 1:(length(xavg)-1)
% 
%         % The averages the rain rate values within a reflectivity bin
%         yavg(i) = mean(y(x<xavg(i+1) & x>=xavg(i)),'all','omitnan');
%         n(i,4) = sum((x<xavg(i+1) & x>=xavg(i)),'all');
% 
%     end
%     yavg(isnan(yavg)) = 0;
% 
%     xtemp = xavg(n(:,4)>0);
%     ytemp = yavg(n(:,4)>0);
%     a = temp(k,n(:,4)>0);
% 
%     testRMSE(k) = sqrt((1/length(ytemp))*sum(abs(ytemp-a').^2));
% end
% 
% b = sprintf('Total RMSE: %1.4f',mean(expRMSE)); % this is the RMSE from ALL the observations
% a = sprintf('Test RMSE: %1.4f',mean(testRMSE)); % This is the RMSE when using the trained exp function on a frac of testing data sites
% 
% text(-23,42.5,a,'fontsize',14,'HorizontalAlignment','left','VerticalAlignment','bottom');
% text(-23,39.2,b,'fontsize',14,'HorizontalAlignment','left','VerticalAlignment','bottom');


%% Model testing with the best fit Z-R relationship 

% Time adjusted + Non-Smoothed

% This portion will take rain rates and average them over their respective
    % reflectivity bin
xavg = -20:5:70;   % min Z is -17.5, max is 65.5 over the whole time
yavg = zeros(length(xavg),1);
n = zeros(length(xavg),4);    % This will provide the number of samples

% Create a random vector with a fraction of the rain gauge points to train the model
% This method assumes that each rain gauge provided equal amounts of data
    % for the model, which is a good assumption since a hurricane is a large
    % scale event
frac = 1;
ri = randperm(size(matchZ,2),round(size(matchZ,2)*frac));
rinot = 1:150;

% This look isolates the sites that are not in the training set. These
    % sites will be used to test the skill of the trained Z-R exp
    % relationship 
varb = 0;
temp = sort(ri);
for i = 1:150
    if rinot(i) == temp(i-varb)
        rinot(i) = 0;
    else
        varb = varb + 1;
    end
end
clear varb temp i

x = matchZ(:,ri);
y = RainGaugedata(matchTimeAdj(:,ri)); % only look at rain gauge data that match the times

% This portion will take rain rates and average them over their respective
    % reflectivity bin

for i = 1:(length(xavg)-1)
    
    % The averages the rain rate values within a reflectivity bin
    yavg(i) = mean(y(x<xavg(i+1) & x>=xavg(i)),'all','omitnan');
    n(i,4) = sum((x<xavg(i+1) & x>=xavg(i)),'all');
    
end
yavg(isnan(yavg)) = 0;

% Significance/Variance testing (itter=1000, w/ repeated)
itter = 1000;
bound = .10/2;  % 90% significance
yavgSig = zeros(length(xavg),2);
temp = zeros(itter,2);

for i = 1:(length(xavg)-1)
    
    % The averages the rain rate values within a reflectivity bin
    a = y(x<xavg(i+1) & x>=xavg(i));
    if isempty(a)
        continue
    else
        for j = 1:itter

            randIndex = randi(length(a),[length(a) 1]);
            temp(j) = mean(a(randIndex),'all','omitnan');
        end
    
    end
    
    temp = sort(temp);
    yavgSig(i,:) =  [temp(itter*bound) temp(itter*(1-bound))];
    
end

f = figure;
bar(xavg,yavg,.75)
hold on
er = errorbar(xavg(n(:,4)>0),yavg(n(:,4)>0),yavgSig((n(:,4)>0),1)-yavg((n(:,4)>0)),yavgSig((n(:,4)>0),2)-yavg((n(:,4)>0)));
er.Color = [0,0,0];
er.LineWidth = 2;
er.LineStyle = 'none';


% % Significane testing the exp fit
% expfitSig = zeros(itter,2); % This will hold the coefficients (a*exp(b*x))
% xavg2 = -20:5:70;   % min Z is -17.5, max is 65.5 over the whole time
% yavg2 = zeros(length(xavg2),1);
% 
% ri = randperm(size(matchZSmooth,2),round(size(matchZSmooth,2)*frac));
% x = matchZSmooth(:,ri);
% y = RainGaugedata(matchTimeAdj(:,ri)); % only look at rain gauge data that match the times
% 
% 
% for i = 1:itter
%         
%     for j = 1:(length(xavg2)-1)
% 
%         % The averages the rain rate values within a reflectivity bin
%         yavg2(j) = mean(y(x<xavg2(j+1) & x>=xavg2(j)),'all','omitnan');
%         n(j,4) = sum((x<xavg2(j+1) & x>=xavg2(j)),'all');
% 
%     end
%     xtemp = xavg(n(:,4)>0);
%     ytemp = yavg2(n(:,4)>0);
%     ztemp = fit(xtemp',ytemp,'exp1');
%     expfitSig(i,:) = coeffvalues(ztemp);    % Gives you coefficients a & b
%     
% end
% yavg(isnan(yavg)) = 0;


% Plot the exponential fit
    %  n(:,4)>0 makes it so the bins with 0 values don't bias the exp fit 
xtemp = xavg(n(:,4)>0);
ytemp = yavg(n(:,4)>0);
expfit = fit(xtemp',ytemp,'exp1');
expcoefaNonSmooth = expfit.a;
expcoefbNonSmooth = expfit.b;
p = plot(xavg,expfit(xavg));

a = expfit(xtemp);
% Root-mean Square Error
expRMSE = sqrt((1/length(ytemp))*sum(abs(ytemp-a).^2));

b = sprintf('Total RMSE: %1.4f',mean(expRMSE)); % this is the RMSE from ALL the observations
text(-23,40,b,'fontsize',14,'HorizontalAlignment','left','VerticalAlignment','bottom');

% % Plotting the 95% significance (don't need to include)
% temp = confint(expfit);
% a = temp(:,1); b = temp(:,2);
% plow = plot(xavg,(a(1)*exp(b(1)*xavg)),':');
% phigh = plot(xavg,(a(2)*exp(b(2)*xavg)),'--');
p.LineWidth = 3;
p.Color = '#D95319';
legend(p,sprintf('%1.4f*exp(%1.4fx)',expfit.a,expfit.b),'Location','northwest')

ax = gca;
ax.FontSize = 16;
ax.Box = 'on';
ax.XGrid = 'on';
ax.YGrid = 'on';
xticks(-20:10:70)
xlabel('Reflectivity [dBZ]')
ylabel('Rain Gauge Rain Rate [mm/hr]')
ylim([0,50])
xlim([-25,70])
title('Rain Rate and Reflectivity Relationship')
subtitle('Time Adjusted, Non-Smoothed')

% clear a c i j x y f p ax temp itter er ans xavg yavg xavgSig yavgSig bound 


% Model testing with the best fit Z-R relationship 

% This is a little complex, but essentially, we are going to select a
    % fraction of the sites to "train" the exp function, and then we will
    % apply that function to the fraction of sites that weren't trained on
    % to see how well they predict/quantify precipitation over those sites.
    % We also want to be able to have a diverse sample of sites that were and
    % weren't trained on so this will be looped over 'itterk' times to be
    % more robust.
% We will be utilizing the time adjusted + smoothed data set as that had
    % the best fit

itterk = 100;
expcoefa = zeros(itterk,1);
expcoefb = zeros(itterk,1);
expRMSE = zeros(itterk,1);

% This will keep track of the sites that were not used so we can
    % calculate how well the exp function preformed on them
rinot = repmat(1:150,itterk,1);
for k = 1:itterk
    
    % This portion will take rain rates and average them over their respective
        % reflectivity bin
    xavg = -20:5:70;   % min Z is -17.5, max is 65.5 over the whole time
    yavg = zeros(length(xavg),1);
%     n = zeros(length(xavg),1);    % This variable provides the number of samples

    % Below we will crete a random vector with a fraction of the rain gauge 
        % points to train the model.
    % This method assumes that each rain gauge provided equal amounts of data
        % for the model, which is a good assumption since a hurricane is a large
        % scale event
    frac = .75;     % The fraction of total sites we will use
    % This creates a random non-repeated vector
    ri = randperm(size(matchZ,2),round(size(matchZ,2)*frac));
   
    if frac~=1
        % This look isolates the sites that are not in the training set. These
            % sites will be used to test the skill of the trained Z-R exp
            % relationship 
        varb = 0;
        temp = sort(ri);
        for i = 1:max(temp)
    %         if rinot(i-1)==150
    %             continue
            if rinot(k,i) == temp(i-varb)
                rinot(k,i) = 0;
        %         varb = varb + 1;
            else
                varb = varb + 1;
            end
        end
        clear varb temp i
    end

    x = matchZ(:,ri);
    y = RainGaugedata(matchTimeAdj,:); % only look at rain gauge data that match the times

    % This portion will take rain rates and average them over their respective
        % reflectivity bin

    for i = 1:(length(xavg)-1)

        % The averages the rain rate values within a reflectivity bin
        yavg(i) = mean(y(x<xavg(i+1) & x>=xavg(i)),'all','omitnan');
        n(i,4) = sum((x<xavg(i+1) & x>=xavg(i)),'all');

    end
    yavg(isnan(yavg)) = 0;

    
    xtemp = xavg(n(:,4)>0);
    ytemp = yavg(n(:,4)>0);
    expfit = fit(xtemp',ytemp,'exp1');
    % Save the coefficients of the 
    expcoefa(k) = expfit.a;
    expcoefb(k) = expfit.b;
    a = expfit(xtemp);
    % Root-mean Square Error
    expRMSE(k) = sqrt((1/length(ytemp))*sum(abs(ytemp-a).^2));
    
end

testRMSE = zeros(itterk,1);
temp = expcoefa.*exp(expcoefb*xavg);
for k = 1:itterk
    x = matchZ(:);
    y = RainGaugedata(matchTimeAdj,:); % only look at rain gauge data that match the times

    % This portion will take rain rates and average them over their respective
        % reflectivity bin

    for i = 1:(length(xavg)-1)

        % The averages the rain rate values within a reflectivity bin
        yavg(i) = mean(y(x<xavg(i+1) & x>=xavg(i)),'all','omitnan');
        n(i,4) = sum((x<xavg(i+1) & x>=xavg(i)),'all');

    end
    yavg(isnan(yavg)) = 0;

    xtemp = xavg(n(:,4)>0);
    ytemp = yavg(n(:,4)>0);
    a = temp(k,n(:,4)>0);

    testRMSE(k) = sqrt((1/length(ytemp))*sum(abs(ytemp-a').^2));
end

% b = sprintf('Total RMSE: %1.4f',mean(expRMSE)); % this is the RMSE from ALL the observations
a = sprintf('Test RMSE: %1.4f',mean(testRMSE)); % This is the RMSE when using the trained exp function on a frac of testing data sites

text(-23,42.5,a,'fontsize',14,'HorizontalAlignment','left','VerticalAlignment','bottom');
% text(-23,40,b,'fontsize',14,'HorizontalAlignment','left','VerticalAlignment','bottom');

% Save figure
% exportgraphics(f,'Time-Adjusted_Non-Smoothed.pdf','ContentType','vector')

%% Qualitatively, view the estimated precipitation 

f = figure('Position',[-210 1340 716 562]);

% Plot reflectivity factor (Z) first
% subplot(1,3,1);
ax = gca;
i = 331;
contourf(lonZ{i},latZ{i},ZData{i}','Linestyle','None')
axis square
hold on;

% Plot radius in minutes it takes rain to fall to ground
% (1-to-5) minutes * 60 seconds * 8 m/s (average terminal velocity of rain drops)
    % The above gives you the meters above the ground the rain drop
    % needs to be for it to fall in 1-to-5 minutes.
% Take invtan(0.46)

% 1 deg = 110.574 km · Longitude: 1 deg = 111.320*cos(latitude)
minute = 1;
pos1 = [lonZ{i}(centerCoord(i,1)) - (minute*60*8/1000)/tand(0.46)/111 ...
    latZ{i}(centerCoord(i,2)) - (minute*60*8/1000)/tand(0.46)/111 ...
    ((minute*60*8/1000)/tand(0.46)/111)*2 ...
    (minute*60*8/1000)/tand(0.46)/111*2];  % [x y w v]
r1 = rectangle('Position',pos1,'Curvature',[1 1],'LineStyle','-','LineWidth',2);

minute = 2;
pos2 = [lonZ{i}(centerCoord(i,1)) - (minute*60*8/1000)/tand(0.46)/111 ...
    latZ{i}(centerCoord(i,2)) - (minute*60*8/1000)/tand(0.46)/111 ...
    ((minute*60*8/1000)/tand(0.46)/111)*2 ...
    (minute*60*8/1000)/tand(0.46)/111*2];  % [x y w v]
r2 = rectangle('Position',pos2,'Curvature',[1 1],'LineStyle','-','LineWidth',2);

minute = 3;
pos3 = [lonZ{i}(centerCoord(i,1)) - (minute*60*8/1000)/tand(0.46)/111 ...
    latZ{i}(centerCoord(i,2)) - (minute*60*8/1000)/tand(0.46)/111 ...
    ((minute*60*8/1000)/tand(0.46)/111)*2 ...
    (minute*60*8/1000)/tand(0.46)/111*2];  % [x y w v]
r3 = rectangle('Position',pos3,'Curvature',[1 1],'LineStyle','-','LineWidth',2);

%     % Plots the center of the radar!
%     scatter(lonZ{i}(centerCoord(i,1)),latZ{i}(centerCoord(i,2)),150,'.','MarkerEdgeColor','m','LineWidth',2);
% Plot the locations of the rain gauges
a = scatter(GaugeLocations.longitude,GaugeLocations.latitude,15,'*','r');

hold off
ax.FontSize = 16;
c = colorbar;
c.Label.String = 'Reflectivity [dBZ]';
c.Label.FontSize = 16;
c.Limits = [0,80];
caxis([0,80]);
colormap(boonlib('zmap',[0 80]))
ax.YTick = 27:1:32;
ax.XTick = -98:1:-92;
yticklabels({'27\circN','28\circN','29\circN','30\circN','31\circN','32\circN'})
xticklabels({'98\circW','97\circW','96\circW','95\circW','94\circW','93\circW','92\circW'})
xlabel('Longitude')
ylabel('Latitude')
title('KHGX - Reflectivity')
t = strcat(char(datetime(timeRadar(i)/86400,'ConvertFrom','datenum')),' UTC');
%     legend([a, r1, r2, r3], ['Rain Gauges', '1-minute radius','2-minute radius','3-minute radius'])
legend(a,'Rain Gauges')
subtitle(t)

% Save Figure
% exportgraphics(f,'Radar_Reflectivity-Quicklook.pdf','ContentType','vector')

%%
% Time Adjusted, non-smoothed
f = figure('Position',[-210 1322 1500 580]);
subplot(1,2,1);
ax = gca;
i = 331;
z = expcoefaNonSmooth*exp(expcoefbNonSmooth*ZData{i});
contourf(lonZ{i},latZ{i},z','Linestyle','None')

axis square
hold on;

% Plot radius in minutes it takes rain to fall to ground
% (1-to-5) minutes * 60 seconds * 8 m/s (average terminal velocity of rain drops)
    % The above gives you the meters above the ground the rain drop
    % needs to be for it to fall in 1-to-5 minutes.
% Take invtan(0.46)

% 1 deg = 110.574 km · Longitude: 1 deg = 111.320*cos(latitude)
minute = 1;
pos1 = [lonZ{i}(centerCoord(i,1)) - (minute*60*8/1000)/tand(0.46)/111 ...
    latZ{i}(centerCoord(i,2)) - (minute*60*8/1000)/tand(0.46)/111 ...
    ((minute*60*8/1000)/tand(0.46)/111)*2 ...
    (minute*60*8/1000)/tand(0.46)/111*2];  % [x y w v]
r1 = rectangle('Position',pos1,'Curvature',[1 1],'LineStyle','-','LineWidth',2);

minute = 2;
pos2 = [lonZ{i}(centerCoord(i,1)) - (minute*60*8/1000)/tand(0.46)/111 ...
    latZ{i}(centerCoord(i,2)) - (minute*60*8/1000)/tand(0.46)/111 ...
    ((minute*60*8/1000)/tand(0.46)/111)*2 ...
    (minute*60*8/1000)/tand(0.46)/111*2];  % [x y w v]
r2 = rectangle('Position',pos2,'Curvature',[1 1],'LineStyle','-','LineWidth',2);

minute = 3;
pos3 = [lonZ{i}(centerCoord(i,1)) - (minute*60*8/1000)/tand(0.46)/111 ...
    latZ{i}(centerCoord(i,2)) - (minute*60*8/1000)/tand(0.46)/111 ...
    ((minute*60*8/1000)/tand(0.46)/111)*2 ...
    (minute*60*8/1000)/tand(0.46)/111*2];  % [x y w v]
r3 = rectangle('Position',pos3,'Curvature',[1 1],'LineStyle','-','LineWidth',2);

%     % Plots the center of the radar!
%     scatter(lonZ{i}(centerCoord(i,1)),latZ{i}(centerCoord(i,2)),150,'.','MarkerEdgeColor','m','LineWidth',2);
% Plot the locations of the rain gauges
a = scatter(GaugeLocations.longitude,GaugeLocations.latitude,15,'*','r');

hold off
ax.FontSize = 16;
c = colorbar;
c.Label.String = 'Rain Rate [mm/hr]';
c.Label.FontSize = 16;
c.Limits = [0,40];
caxis([0,40]);
colormap(turbo(80))
ax.YTick = 27:1:32;
ax.XTick = -98:1:-92;
yticklabels({'27\circN','28\circN','29\circN','30\circN','31\circN','32\circN'})
xticklabels({'98\circW','97\circW','96\circW','95\circW','94\circW','93\circW','92\circW'})
xlabel('Longitude')
ylabel('Latitude')
title('KHGX - Precipiation Rate w/ Time Adjusted, Non-Smooth')
t = strcat(char(datetime(timeRadar(i)/86400,'ConvertFrom','datenum')),' UTC');
%     legend([a, r1, r2, r3], ['Rain Gauges', '1-minute radius','2-minute radius','3-minute radius'])
legend(a,'Rain Gauges')
subtitle(t)

subplot(1,2,2);
ax = gca;
i = 331;
z = expcoefaSmooth*exp(expcoefbSmooth*ZData{i});
contourf(lonZ{i},latZ{i},z','Linestyle','None')

axis square
hold on;

% Plot radius in minutes it takes rain to fall to ground
% (1-to-5) minutes * 60 seconds * 8 m/s (average terminal velocity of rain drops)
    % The above gives you the meters above the ground the rain drop
    % needs to be for it to fall in 1-to-5 minutes.
% Take invtan(0.46)

% 1 deg = 110.574 km · Longitude: 1 deg = 111.320*cos(latitude)
minute = 1;
pos1 = [lonZ{i}(centerCoord(i,1)) - (minute*60*8/1000)/tand(0.46)/111 ...
    latZ{i}(centerCoord(i,2)) - (minute*60*8/1000)/tand(0.46)/111 ...
    ((minute*60*8/1000)/tand(0.46)/111)*2 ...
    (minute*60*8/1000)/tand(0.46)/111*2];  % [x y w v]
r1 = rectangle('Position',pos1,'Curvature',[1 1],'LineStyle','-','LineWidth',2);

minute = 2;
pos2 = [lonZ{i}(centerCoord(i,1)) - (minute*60*8/1000)/tand(0.46)/111 ...
    latZ{i}(centerCoord(i,2)) - (minute*60*8/1000)/tand(0.46)/111 ...
    ((minute*60*8/1000)/tand(0.46)/111)*2 ...
    (minute*60*8/1000)/tand(0.46)/111*2];  % [x y w v]
r2 = rectangle('Position',pos2,'Curvature',[1 1],'LineStyle','-','LineWidth',2);

minute = 3;
pos3 = [lonZ{i}(centerCoord(i,1)) - (minute*60*8/1000)/tand(0.46)/111 ...
    latZ{i}(centerCoord(i,2)) - (minute*60*8/1000)/tand(0.46)/111 ...
    ((minute*60*8/1000)/tand(0.46)/111)*2 ...
    (minute*60*8/1000)/tand(0.46)/111*2];  % [x y w v]
r3 = rectangle('Position',pos3,'Curvature',[1 1],'LineStyle','-','LineWidth',2);

%     % Plots the center of the radar!
%     scatter(lonZ{i}(centerCoord(i,1)),latZ{i}(centerCoord(i,2)),150,'.','MarkerEdgeColor','m','LineWidth',2);
% Plot the locations of the rain gauges
a = scatter(GaugeLocations.longitude,GaugeLocations.latitude,15,'*','r');

hold off
ax.FontSize = 16;
c = colorbar;
c.Label.String = 'Rain Rate [mm/hr]';
c.Label.FontSize = 16;
c.Limits = [0,40];
caxis([0,40]);
colormap(turbo(80))
ax.YTick = 27:1:32;
ax.XTick = -98:1:-92;
yticklabels({'27\circN','28\circN','29\circN','30\circN','31\circN','32\circN'})
xticklabels({'98\circW','97\circW','96\circW','95\circW','94\circW','93\circW','92\circW'})
xlabel('Longitude')
ylabel('Latitude')
title('KHGX - Precipiation Rate w/ Time Adjusted, Smoothed')
t = strcat(char(datetime(timeRadar(i)/86400,'ConvertFrom','datenum')),' UTC');
%     legend([a, r1, r2, r3], ['Rain Gauges', '1-minute radius','2-minute radius','3-minute radius'])
legend(a,'Rain Gauges')
subtitle(t)

% Save figure
% exportgraphics(f,'Quantifying_Precip-QuickLook.pdf','ContentType','vector')

%% Compare and contrast between the different methods ZDR

f = figure('Position',[-210 1060 1047 842]);

% Non-time adjusted + non-smooth
subplot(2,2,1)

% x = matchZDR(:,GaugeLocations.radarDistance>20);
% y = RainGaugedata(matchTime,GaugeLocations.radarDistance>20); % only look at rain gauge data that match the times
% x = matchZDR(:,std(matchZDR,1,'omitnan')<3);
% y = RainGaugedata(matchTime,std(matchZDR,1,'omitnan')<3); % only look at rain gauge data that match the times
% x = matchZDR(:,mean(matchZDR,1,'omitnan')>0);
% y = RainGaugedata(matchTime,mean(matchZDR,1,'omitnan')>0); % only look at rain gauge data that match the times
x = matchZDR;
y = RainGaugedata(matchTime,:); % only look at rain gauge data that match the times

% This portion will take rain rates and average them over their respective
    % Differential Reflectivity bin
xavg = -10:1:10;   % min ZDR is -7, max is 7 over the whole time
yavg = zeros(length(xavg),1);
n = zeros(length(xavg),4);    % This will provide the number of samples

for i = 1:(length(xavg)-1)
    
    % The averages the rain rate values within a Differential Reflectivity bin
    yavg(i) = mean(y(x<xavg(i+1) & x>=xavg(i)),'all','omitnan');
    
    % Counts the number of samples to check the potential of outliers
    n(i,1) = sum((x<xavg(i+1) & x>=xavg(i)),'all');
    
end
yavg(isnan(yavg)) = 0;

% Significance/Variance testing through the Monte Carlo method (itter=1000, w/ repetition)
itter = 1000;
bound = .10/2;  % 90% significance
yavgSig = zeros(length(xavg),2);
temp = zeros(itter,2);

for i = 1:(length(xavg)-1)
    
    % The averages the rain rate values within a Differential Reflectivity bin
    a = y(x<xavg(i+1) & x>=xavg(i));
    if isempty(a)
        continue
    else
        for j = 1:itter

            randIndex = randi(length(a),[length(a) 1]);
            temp(j) = mean(a(randIndex),'all','omitnan');
        end
    end
    
    temp = sort(temp);
    yavgSig(i,:) =  [temp(itter*bound) temp(itter*(1-bound))];
    
end

bar(xavg,yavg,.75)
hold on
% Plot the confidence intervals at each bin
er = errorbar(xavg(n(:,1)~=0),yavg(n(:,1)~=0),yavgSig((n(:,1)~=0),1)-yavg((n(:,1)~=0)),yavgSig((n(:,1)~=0),2)-yavg((n(:,1)~=0)));
er.Color = [0,0,0];
er.LineWidth = 1.5;
er.LineStyle = 'none';

% % Plot the exponential fit
% xtemp = xavg(n(:,1)~=0);
% ytemp = yavg(n(:,1)~=0);
% expfit = fit(xtemp',ytemp,'exp1');
% p = plot(xavg,expfit(xavg));
% p.LineWidth = 3;
% p.Color = '#D95319';
% legend(p,sprintf('%1.4f*exp(%1.4fx)',expfit.a,expfit.b),'Location','northeast')

ax = gca;
ax.FontSize = 16;
ax.Box = 'on';
ax.XGrid = 'on';
ax.YGrid = 'on';
xticks(-10:2:10)
xlabel('Differential Reflectivity [dBZ]')
ylabel('Rain Gauge Rain Rate [mm/hr]')
ylim([0,20])
xlim([-10,10])
title('Rain Rate and Differential Reflectivity Relationship')
subtitle('Non-Time Adjusted, Non-Smooth')


% Time-adjusted + non-smooth
subplot(2,2,2)
x = matchZDR;
y = RainGaugedata(matchTimeAdj,:); % only look at rain gauge data that match the times

% This portion will take rain rates and average them over their respective
    % Differential Reflectivity bin
    
for i = 1:(length(xavg)-1)
    
    % The averages the rain rate values within a Differential Reflectivity bin
    yavg(i) = mean(y(x<xavg(i+1) & x>=xavg(i)),'all','omitnan');
    n(i,2) = sum((x<xavg(i+1) & x>=xavg(i)),'all');
    
end
yavg(isnan(yavg)) = 0;

% Significance/Variance testing (itter=1000, w/ repeated)
yavgSig = zeros(length(xavg),2);
temp = zeros(itter,2);

for i = 1:(length(xavg)-1)
    
    % The averages the rain rate values within a Differential Reflectivity bin
    a = y(x<xavg(i+1) & x>=xavg(i));
    if isempty(a)
        continue
    else
        for j = 1:itter

            randIndex = randi(length(a),[length(a) 1]);
            temp(j) = mean(a(randIndex),'all','omitnan'); 
        end
    end
    
    temp = sort(temp);
    yavgSig(i,:) =  [temp(itter*bound) temp(itter*(1-bound))];
    
end

bar(xavg,yavg,.75)
hold on

er = errorbar(xavg(n(:,2)~=0),yavg(n(:,2)~=0),yavgSig((n(:,2)~=0),1)-yavg((n(:,2)~=0)),yavgSig((n(:,2)~=0),2)-yavg((n(:,2)~=0)));
er.Color = [0,0,0];
er.LineWidth = 2;
er.LineStyle = 'none';

% % Plot the exponential fit
% xtemp = xavg(n(:,2)~=0);
% ytemp = yavg(n(:,2)~=0);
% expfit = fit(xtemp',ytemp,'exp1');
% p = plot(xavg,expfit(xavg));
% p.LineWidth = 3;
% p.Color = '#D95319';
% legend(p,sprintf('%1.4f*exp(%1.4fx)',expfit.a,expfit.b),'Location','northwest')

ax = gca;
ax.FontSize = 16;
ax.Box = 'on';
ax.XGrid = 'on';
ax.YGrid = 'on';
xticks(-10:2:10)
xlabel('Differential Reflectivity [dBZ]')
ylabel('Rain Gauge Rain Rate [mm/hr]')
ylim([0,20])
xlim([-10,10])
title('Rain Rate and Differential Reflectivity Relationship')
subtitle('Time Adjusted, Non-Smooth')


% Non-time adjusted + Smoothed
subplot(2,2,3)

x = matchZDRSmooth;
y = RainGaugedata(matchTime,:); % only look at rain gauge data that match the times

% This portion will take rain rates and average them over their respective
    % Differential Reflectivity bin

for i = 1:(length(xavg)-1)
    
    % The averages the rain rate values within a Differential Reflectivity bin
    yavg(i) = mean(y(x<xavg(i+1) & x>=xavg(i)),'all','omitnan');
    n(i,3) = sum((x<xavg(i+1) & x>=xavg(i)),'all');
    
end
yavg(isnan(yavg)) = 0;
% Significance/Variance testing (itter=1000, w/ repeated)
yavgSig = zeros(length(xavg),2);
temp = zeros(itter,2);

for i = 1:(length(xavg)-1)
    
    % The averages the rain rate values within a Differential Reflectivity bin
    a = y(x<xavg(i+1) & x>=xavg(i));
    
    if isempty(a)
        continue
    else
        for j = 1:itter

            randIndex = randi(length(a),[length(a) 1]);
            temp(j) = mean(a(randIndex),'all','omitnan');

        end
    end
    
    temp = sort(temp);
    yavgSig(i,:) =  [temp(itter*bound) temp(itter*(1-bound))];
    
end

bar(xavg,yavg,.75)
hold on
er = errorbar(xavg(n(:,3)~=0),yavg(n(:,3)~=0),yavgSig((n(:,3)~=0),1)-yavg((n(:,3)~=0)),yavgSig((n(:,3)~=0),2)-yavg((n(:,3)~=0)));
er.Color = [0,0,0];
er.LineWidth = 2;
er.LineStyle = 'none';

% % Plot the exponential fit
% xtemp = xavg(n(:,3)~=0);
% ytemp = yavg(n(:,3)~=0);
% expfit = fit(xtemp',ytemp,'exp1');
% p = plot(xavg,expfit(xavg));
% p.LineWidth = 3;
% p.Color = '#D95319';
% legend(p,sprintf('%1.4f*exp(%1.4fx)',expfit.a,expfit.b),'Location','northeast')

ax = gca;
ax.FontSize = 16;
ax.Box = 'on';
ax.XGrid = 'on';
ax.YGrid = 'on';
xticks(-10:2:10)
xlabel('Differential Reflectivity [dBZ]')
ylabel('Rain Gauge Rain Rate [mm/hr]')
ylim([0,20])
xlim([-10,10])
title('Rain Rate and Differential Reflectivity Relationship')
subtitle('Non-Time Adjusted, Smoothed')


% Time adjusted + Smoothed
subplot(2,2,4)

x = matchZDRSmooth;
y = RainGaugedata(matchTimeAdj,:); % only look at rain gauge data that match the times

% This portion will take rain rates and average them over their respective
    % Differential Reflectivity bin

for i = 1:(length(xavg)-1)
    
    % The averages the rain rate values within a Differential Reflectivity bin
    yavg(i) = mean(y(x<xavg(i+1) & x>=xavg(i)),'all','omitnan');
    n(i,4) = sum((x<xavg(i+1) & x>=xavg(i)),'all');
    
end
yavg(isnan(yavg)) = 0;

% Significance/Variance testing (itter=1000, w/ repeated)
yavgSig = zeros(length(xavg),2);
temp = zeros(itter,2);

for i = 1:(length(xavg)-1)
    
    % The averages the rain rate values within a Differential Reflectivity bin
    a = y(x<xavg(i+1) & x>=xavg(i));
    if isempty(a)
        continue
    else
        for j = 1:itter

            randIndex = randi(length(a),[length(a) 1]);
            temp(j) = mean(a(randIndex),'all','omitnan');
        end
    
    end
    
    temp = sort(temp);
    yavgSig(i,:) =  [temp(itter*bound) temp(itter*(1-bound))];
    
end

bar(xavg,yavg,.75)
hold on
er = errorbar(xavg(n(:,4)>0),yavg(n(:,4)>0),yavgSig((n(:,4)>0),1)-yavg((n(:,4)>0)),yavgSig((n(:,4)>0),2)-yavg((n(:,4)>0)));
er.Color = [0,0,0];
er.LineWidth = 2;
er.LineStyle = 'none';

% % Plot the exponential fit
% xtemp = xavg(n(:,4)>0);
% ytemp = yavg(n(:,4)>0);
% expfit = fit(xtemp',ytemp,'exp1');
% p = plot(xavg,expfit(xavg));
% p.LineWidth = 3;
% p.Color = '#D95319';
% legend(p,sprintf('%1.4f*exp(%1.4fx)',expfit.a,expfit.b),'Location','northwest')

ax = gca;
ax.FontSize = 16;
ax.Box = 'on';
ax.XGrid = 'on';
ax.YGrid = 'on';
xticks(-10:2:10)
xlabel('Differential Reflectivity [dBZ]')
ylabel('Rain Gauge Rain Rate [mm/hr]')
ylim([0,20])
xlim([-10,10])
title('Rain Rate and Differential Reflectivity Relationship')
subtitle('Time Adjusted, Smoothed')

% clear a i j x y f ax temp itter er ans xavg yavg xavgSig yavgSig bound