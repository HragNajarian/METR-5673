%% Load in the file names
clc
clear
% Make sure pathway is the folder where your folders for reflectivity are
    % located
pathwayZ = '/Users/hragnajarian/Desktop/KHGX (Gridded) - 08:25-29/Z';
pathwayZDR = '/Users/hragnajarian/Desktop/KHGX (Gridded) - 08:25-29/ZDR';
ZFile = dir(pathwayZ);
ZDRFile = dir(pathwayZDR);

ZFile = ZFile(4:end);     % Remove first 3 files
ZDRFile = ZDRFile(3:end);     % Remove first 2 files

clear pathwayZ pathwayZDR
%% Match Z and ZDR files
% This is simply to make sure the files for Z and ZDR
    % match
% Check which files are in common and remove the ones that are not common
for i = 1:length(ZFile)-1
    % See if the file names match, when it doesn't delete that file
    bool = ZFile(i).name == ZDRFile(i).name;    % Should give all 1's if they match
    if ~all(bool)
        ZFile(i) = [];
    end
end

clear bool

%%
tempZ = cell(length(ZFile),1); 
tempZDR = cell(length(ZDRFile),1); 
for i = 1:length(ZFile)
    tempZ{i} = strcat(ZFile(i).folder,'/',ZFile(i).name);
    tempZDR{i} = strcat(ZDRFile(i).folder,'/',ZDRFile(i).name);
end

ZFile = tempZ;
ZDRFile = tempZDR;

fileInfoZ = ncinfo(ZFile{1});
fileInfoZDR = ncinfo(ZDRFile{1});

clear tempZ tempZDR
%% Read in times
% Data starts on August 25th, 2017 and ends August 29th, 2017 and data is
    % collected every 5-6 minutes
% Data is from lowest EL = 0.46°
    
% I only looked at ZFile since the times will match between ZDR and Z files
    
% Function "datenum" gives the number of days since January 0, 0000
daysSince = datenum(1970,1,1)*86400;

% 'time' provides SECONDS since 1970-1-1
baseTime = zeros(length(ZFile),1);
for i = 1:length(ZFile)
    baseTime(i) = ncread(ZFile{i},'time');
end

% This will give us seconds since January 0, 0000
timeRadar = baseTime + daysSince;

% Run this for a quick check to make sure dates are correct
datetime(timeRadar(1)/86400,'ConvertFrom','datenum')

%% Read in basic variables

% First check if the latitudes and long are all the same. They're not but
% they are close enough

% Lat and Lon are also in different sizes depending on file type...wonderful
latZ = cell(length(ZFile),1);
lonZ = cell(length(ZFile),1);
latZDR = cell(length(ZDRFile),1);
lonZDR = cell(length(ZDRFile),1);

for i = 1:length(ZFile)
    latZ{i} = ncread(ZFile{i},'lat');
    lonZ{i} = ncread(ZFile{i},'lon');
    latZDR{i} = ncread(ZDRFile{i},'lat');
    lonZDR{i} = ncread(ZDRFile{i},'lon');
end
    
%% Concatinate all Z and ZDR data into one matrix
ZData = cell(length(ZFile),1);
ZDRData = cell(length(ZDRFile),1);

for i = 1:length(ZFile)
    ZData{i} = ncread(ZFile{i},'Reflectivity');
    ZDRData{i} = ncread(ZDRFile{i},'DifferentialReflectivity');
end

%% Where is the center of the radar?
% There is a radius around the center of the radar that are all NaNs, so
    % this is not the only point, but we have
% lon(1), lat(2) indices
centerCoord = zeros(length(ZData),2);
temp = zeros(length(ZData),1);

for i = 1:length(ZData)
    
    centerCoord(i,1) = size(ZData{i},1)/2 + 0;
    centerCoord(i,2) = size(ZData{i},2)/2 + 4;
    
    % The Z values at the center of the radar, should be all NaN's
    temp(i) = ZData{i}(centerCoord(i,1),centerCoord(i,2));
    
end

% Make sure the reflectivity values at the center make sense, they should 
    % all be NaN's since the radar is not able to see the reflectivity
    % above.
isnan(mean(temp))   % if 1, then good

clear ans temp i

save ('Z&ZDRData.mat','-v7.3')
%% Plot the evolution of reflectivity
set(groot, 'DefaultAxesFontName', 'Arial')
% 1 deg = 110.574 km · Longitude: 1 deg = 111.320*cos(latitude)

f = figure('Position',[-210 1340 716 562]);
apple = sum(isnan(matchZ),1);
% for i = 1:30:length(ZData)
for i = 331
    ax = gca;
    contourf(lonZ{i},latZ{i},ZData{i}','Linestyle','None')
    axis square
    hold on;
%     % Plots the center of the radar!
%     scatter(lonZ{i}(centerCoord(i,1)),latZ{i}(centerCoord(i,2)),150,'.','MarkerEdgeColor','m','LineWidth',2);
    % Plot the locations of the rain gauges
    a = scatter(GaugeLocations.longitude(apple<10),GaugeLocations.latitude(apple<10),75,'.','k');
    
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
    r1 = rectangle('Position',pos1,'Curvature',[1 1],'LineStyle','-');
    
    minute = 2;
    pos2 = [lonZ{i}(centerCoord(i,1)) - (minute*60*8/1000)/tand(0.46)/111 ...
        latZ{i}(centerCoord(i,2)) - (minute*60*8/1000)/tand(0.46)/111 ...
        ((minute*60*8/1000)/tand(0.46)/111)*2 ...
        (minute*60*8/1000)/tand(0.46)/111*2];  % [x y w v]
    r2 = rectangle('Position',pos2,'Curvature',[1 1],'LineStyle','--');
    
    minute = 3;
    pos3 = [lonZ{i}(centerCoord(i,1)) - (minute*60*8/1000)/tand(0.46)/111 ...
        latZ{i}(centerCoord(i,2)) - (minute*60*8/1000)/tand(0.46)/111 ...
        ((minute*60*8/1000)/tand(0.46)/111)*2 ...
        (minute*60*8/1000)/tand(0.46)/111*2];  % [x y w v]
    r3 = rectangle('Position',pos3,'Curvature',[1 1],'LineStyle',':');
    
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
    pause(0.5)
    
end
clear ax c t i ans pos1 pos2 pos3 r1 r2 r3 minute f a

%%
save('Z&ZDRData.mat','-v7.3')