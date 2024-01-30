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


%% INS/GNSS Solution

% Configure final outputs
% outputFile = "Output_Profile";
% % 1: Time
% % 2: Geodetic latitude in degrees
% % 3: Geodetic longitude in degrees 
% % 4: North velocity in m/s
% % 5: East velocity in m/s
% % 6: Heading in degrees.
% outputs = zeros(epochNum,6);

% Write to table and export
% outputTable = table(outputs);
% writetable(outputTable,"Solutions/"+outputFile+".csv",'WriteVariableNames',0)
