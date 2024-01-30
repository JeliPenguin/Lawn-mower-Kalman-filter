clear;

Load_Data;

%% GNSS Solution
gnss_solution = GNSS_KF(pseudo_range_data,pseudo_range_rate_data,times);
gnss_latitude = gnss_solution(:,2);
gnss_longitude = gnss_solution(:,3);
gnss_velocity_n = gnss_solution(:,5);
gnss_velocity_e = gnss_solution(:,6);

Plot_Graph(gnss_longitude,gnss_latitude,gnss_velocity_n,gnss_velocity_e,times,"GNSS",false)


%% DR Solution
dr_solution = Dead_Reckoning(dr_measurement_data,gnss_solution,times);
dr_latitude = dr_solution(:,2);
dr_longitude = dr_solution(:,3);
dr_velocity_n = dr_solution(:,4);
dr_velocity_e = dr_solution(:,5);

Plot_Graph(dr_longitude,dr_latitude,dr_velocity_n,dr_velocity_e,times,"DR",false)

%% DR/GNSS Solution
dr_gnss_solution = DR_GNSS(gnss_solution,dr_solution,times);
dr_gnss_latitude = dr_gnss_solution(:,2);
dr_gnss_longitude = dr_gnss_solution(:,3);
dr_gnss_velocity_n = dr_gnss_solution(:,4);
dr_gnss_velocity_e = dr_gnss_solution(:,5);

Plot_Graph(dr_gnss_longitude,dr_gnss_latitude,dr_gnss_velocity_n,dr_gnss_velocity_e,times,"DR-GNSS",false)

%% INS/GNSS Solution
% ins_gnss_solution = INS_GNSS();

% Plot_Graph(dr_longitude,dr_latitude,dr_velocity_n,dr_velocity_e,times,"DR",false)