% Load in data
pseudo_range_file = readtable("Pseudo_ranges.csv");
pseudo_range_rate_file = readtable("Pseudo_range_rates.csv");

% pseudoRangeFile = readtable("Workshop1_Pseudo_ranges.csv");
% pseudoRangeRateFile = readtable("Workshop1_Pseudo_range_rates.csv");

% pseudoRangeFile = readtable("Workshop2_Pseudo_ranges.csv");
% pseudoRangeRateFile = readtable("Workshop2_Pseudo_range_rates.csv");

satellite_num_row = table2array(pseudo_range_file(1,:));
satellite_numbers = satellite_num_row(2:end);
satellite_count = length(satellite_numbers);

pseudo_ranges = table2array(pseudo_range_file(2:end,2:end));
pseudo_range_rates = table2array(pseudo_range_rate_file(2:end,2:end));
times = table2array(pseudo_range_file(2:end,1));
epoch_num = height(times);