% Load in data from given csv files
pseudo_range_file = readtable("Pseudo_ranges.csv");
pseudo_range_data = table2array(pseudo_range_file);

pseudo_range_rate_file = readtable("Pseudo_range_rates.csv");
pseudo_range_rate_data = table2array(pseudo_range_rate_file);

deadReckoningFile = readtable("Dead_Reckoning.csv");
dr_measurement_data = table2array(deadReckoningFile);

% Array of time value for each epoch, being the first column of the files
times = pseudo_range_data(2:end,1); 




