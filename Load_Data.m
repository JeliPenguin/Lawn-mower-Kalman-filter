% Load in data
pseudo_range_file = readtable("Pseudo_ranges.csv");
pseudo_range_rate_file = readtable("Pseudo_range_rates.csv");
pseudo_range_data = table2array(pseudo_range_file);
pseudo_range_rate_data = table2array(pseudo_range_rate_file);

times = pseudo_range_data(2:end,1); % Array of time value for each epoch


deadReckoningFile = readtable("Dead_Reckoning.csv");
dr_measurement_data = table2array(deadReckoningFile);

