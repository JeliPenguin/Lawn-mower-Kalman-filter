clear;

Load_Data;

%% GNSS Solution
gnss_solution = GNSS_KF(pseudo_range_data,pseudo_range_rate_data,times);
gnss_latitude = gnss_solution(:,2);
gnss_longitude = gnss_solution(:,3);
gnss_velocity_n = gnss_solution(:,5);
gnss_velocity_e = gnss_solution(:,6);

%% Gyro-smoothed Magnetic Heading solution
gyro_mag_heading_solution = Gyro_Magnetometer_KF(dr_measurement_data,times);
% gyro_mag_heading_solution = dr_measurement_data(:,7);

%% DR Solution
dr_solution = Dead_Reckoning(dr_measurement_data,gnss_solution,gyro_mag_heading_solution,times);
dr_latitude = dr_solution(:,2);
dr_longitude = dr_solution(:,3);
dr_velocity_n = dr_solution(:,4);
dr_velocity_e = dr_solution(:,5);

%% DR-GNSS Solution
dr_gnss_solution = DR_GNSS(gnss_solution,dr_solution,times);
dr_gnss_latitude = dr_gnss_solution(:,2);
dr_gnss_longitude = dr_gnss_solution(:,3);
dr_gnss_velocity_n = dr_gnss_solution(:,4);
dr_gnss_velocity_e = dr_gnss_solution(:,5);


%% Plotting
% figure;
% plot(gnss_longitude,gnss_latitude);
% hold on
% plot(dr_longitude,dr_latitude);
% hold on
% plot(dr_gnss_longitude,dr_gnss_latitude);
% legend(["GNSS","DR","GR-DNSS"]);

Plot_Graph(gnss_longitude,gnss_latitude,gnss_velocity_n,gnss_velocity_e,times,"GNSS",false)

Plot_Graph(dr_longitude,dr_latitude,dr_velocity_n,dr_velocity_e,times,"DR",false)

Plot_Graph(dr_gnss_longitude,dr_gnss_latitude,dr_gnss_velocity_n,dr_gnss_velocity_e,times,"DR-GNSS",false)

Plot_heading(dr_gnss_longitude,dr_gnss_latitude,gyro_mag_heading_solution,false);