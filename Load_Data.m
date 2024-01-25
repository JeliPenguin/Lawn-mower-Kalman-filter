% Load in data
pseudoRangeFile = readtable("Pseudo_ranges.csv");
pseudoRangeRateFile = readtable("Pseudo_range_rates.csv");
% pseudoRangeFile = readtable("Workshop1_Pseudo_ranges.csv");
% pseudoRangeRateFile = readtable("Workshop1_Pseudo_range_rates.csv");

satelliteNumRow = table2array(pseudoRangeFile(1,:));
satelliteNumbers = satelliteNumRow(2:end);
numberOfSatellites = length(satelliteNumbers);

pseudoRanges = table2array(pseudoRangeFile(2:end,2:end));
pseudoRangeRates = table2array(pseudoRangeRateFile(2:end,2:end));
times = table2array(pseudoRangeFile(2:end,1));
epochNum = height(times);