%% Import RainGaugedata

% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 152);

% Specify sheet and range
opts.Sheet = "Sheet1 (2)";
opts.DataRange = "A9:EV1452";

% Specify column names and types
opts.VariableNames = ["startTime", "endTime", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15", "VarName16", "VarName17", "VarName18", "VarName19", "VarName20", "VarName21", "VarName22", "VarName23", "VarName24", "VarName25", "VarName26", "VarName27", "VarName28", "VarName29", "VarName30", "VarName31", "VarName32", "VarName33", "VarName34", "VarName35", "VarName36", "VarName37", "VarName38", "VarName39", "VarName40", "VarName41", "VarName42", "VarName43", "VarName44", "VarName45", "VarName46", "VarName47", "VarName48", "VarName49", "VarName50", "VarName51", "VarName52", "VarName53", "VarName54", "VarName55", "VarName56", "VarName57", "VarName58", "VarName59", "VarName60", "VarName61", "VarName62", "VarName63", "VarName64", "VarName65", "VarName66", "VarName67", "VarName68", "VarName69", "VarName70", "VarName71", "VarName72", "VarName73", "VarName74", "VarName75", "VarName76", "VarName77", "VarName78", "VarName79", "VarName80", "VarName81", "VarName82", "VarName83", "VarName84", "VarName85", "VarName86", "VarName87", "VarName88", "VarName89", "VarName90", "VarName91", "VarName92", "VarName93", "VarName94", "VarName95", "VarName96", "VarName97", "VarName98", "VarName99", "VarName100", "VarName101", "VarName102", "VarName103", "VarName104", "VarName105", "VarName106", "VarName107", "VarName108", "VarName109", "VarName110", "VarName111", "VarName112", "VarName113", "VarName114", "VarName115", "VarName116", "VarName117", "VarName118", "VarName119", "VarName120", "VarName121", "VarName122", "VarName123", "VarName124", "VarName125", "VarName126", "VarName127", "VarName128", "VarName129", "VarName130", "VarName131", "VarName132", "VarName133", "VarName134", "VarName135", "VarName136", "VarName137", "VarName138", "VarName139", "VarName140", "VarName141", "VarName142", "VarName143", "VarName144", "VarName145", "VarName146", "VarName147", "VarName148", "VarName149", "VarName150", "VarName151", "VarName152"];
opts.VariableTypes = ["datetime", "datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "startTime", "InputFormat", "");
opts = setvaropts(opts, "endTime", "InputFormat", "");

% Import the data
RainGaugedata = readtable("/Users/hragnajarian/Library/Mobile Documents/com~apple~CloudDocs/OU/Classes/METR 5673-Weather Radar Theory and Practice/Final Project/Rain Gauge data.xlsx", opts, "UseExcel", false);


% Clear temporary variables
clear opts

%% Quality of life changes
% Remove the first, thrid, and fourth rows

RainGaugedata(4,:) = [];
RainGaugedata(3,:) = [];
RainGaugedata(1,:) = [];

% Site numbers are row 1
siteNum = table2array(RainGaugedata(1,3:end))';

% Gather end times, we don't need start times
endTime = table2array(RainGaugedata(2:end,2));
% Convert to datenum to make comparison between radar times and rain gauge
    % times easier
endTime = datenum(endTime)*86400;

% % Remove the start and end times
RainGaugedata(:,2) = [];
RainGaugedata(:,1) = [];

% Remove the site numbers
RainGaugedata(1,:) = [];

% Make it into an array
RainGaugedata = table2array(RainGaugedata);

% Sort it to be accending
[~,b] = sort(siteNum);
temp = zeros(size(RainGaugedata));
for i = 1:size(RainGaugedata,2)
    temp(:,i) =  RainGaugedata(:,b(i));
end

RainGaugedata = temp;

% Convert to mm/5min from in/5min
RainGaugedata = RainGaugedata*25.4; % 1 in = 25.4 mm
% Convert to mm/hr
RainGaugedata = RainGaugedata*12;
% Minimum reading will now be 1 mm with a delta of 1 mm

[~,y] = max(RainGaugedata(:,43));
% Remove the outlier
RainGaugedata(y,43) = 0;


clear b i temp y
%% Import Site Locations
% There will be extra site locations FYI

opts = spreadsheetImportOptions("NumVariables", 3);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:C187";

% Specify column names and types
opts.VariableNames = ["siteNum", "latitude", "longitude"];
opts.VariableTypes = ["double", "double", "double"];

% Import the data
GaugeLocations = readtable("/Users/hragnajarian/Library/Mobile Documents/com~apple~CloudDocs/OU/Classes/METR 5673-Weather Radar Theory and Practice/Final Project/Site Locations.xlsx", opts, "UseExcel", false);
GaugeLocations.longitude(180) = -95.532730;
GaugeLocations.longitude(31) = -95.339294;

% Clear temporary variables
clear opts

%% Quality of life changes for site locations

% Only keep the site locations that are in 'siteNum'

for i = length(GaugeLocations.siteNum):-1:1

    if sum(GaugeLocations.siteNum(i)==siteNum) == 1
        continue
    else
        GaugeLocations(i,:) = [];
    end
end

clear siteNum i

save('RainGaugeData.mat')