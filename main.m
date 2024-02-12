clear;

Load_Data;

%% GNSS With Outlier Solution
gnss_outlier_solution = GNSS_KF(pseudo_range_data,pseudo_range_rate_data,times,false);
gnss_outlier_latitude = gnss_outlier_solution(:,2);
gnss_outlier_longitude = gnss_outlier_solution(:,3);
gnss_outlier_velocity_n = gnss_outlier_solution(:,5);
gnss_outlier_velocity_e = gnss_outlier_solution(:,6);

%% GNSS Outlier Removed Solution
gnss_solution = GNSS_KF(pseudo_range_data,pseudo_range_rate_data,times,true);
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

%% Final Output
final_output = [dr_gnss_solution,gyro_mag_heading_solution];
outputTable = array2table(final_output);
writetable(outputTable,"Solutions/final_solution.csv",'WriteVariableNames',0)


%% Plotting
% figure;
% plot(times, gnss_velocity_e,"xr-");
% hold on
% plot(times, dr_velocity_e,"*g-");
% hold on
% plot(times,dr_gnss_velocity_e,"b");
% grid on
% title("East Velocities")
% xlabel("Time (s)");
% ylabel("Velocity (m/s)");
% legend(["GNSS","DR","GR-DNSS"]);

% figure;
% plot(gnss_longitude,gnss_latitude,"b-o");
% hold on
% plot(gnss_outlier_longitude,gnss_outlier_latitude);
% title("Position");
% xlabel("Longitude (deg)");
% ylabel("Latitude (deg)");
% grid on
% legend(["GNSS Outlier Removed","GNSS With Outlier"]);

Plot_Graph(gnss_longitude,gnss_latitude,gnss_velocity_n,gnss_velocity_e,times,"GNSS",false)

Plot_Graph(gnss_outlier_longitude,gnss_outlier_latitude,gnss_outlier_velocity_n,gnss_outlier_velocity_e,times,"GNSS-Outlier",false)

Plot_Graph(dr_longitude,dr_latitude,dr_velocity_n,dr_velocity_e,times,"DR",false)

Plot_Graph(dr_gnss_longitude,dr_gnss_latitude,dr_gnss_velocity_n,dr_gnss_velocity_e,times,"DR-GNSS",false)

Plot_heading(dr_gnss_longitude,dr_gnss_latitude,gyro_mag_heading_solution,false);