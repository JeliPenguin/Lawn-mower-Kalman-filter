clear;

Load_Data;

% Configure outputs
outputFile = "Output_Profile";
% 1: Time
% 2: Geodetic latitude in degrees
% 3: Geodetic longitude in degrees 
% 4: North velocity in m/s
% 5: East velocity in m/s
% 6: Heading in degrees.
% outputs = zeros(epochNum,6);


gnss_solution = GNSS_KF(pseudo_ranges,pseudo_range_rates,satellite_numbers,satellite_count,times,epoch_num);
gnss_latitude = gnss_solution(:,2);
gnss_longitude = gnss_solution(:,3);
gnss_velocity_n = gnss_solution(:,5);
gnss_velocity_e = gnss_solution(:,6);

Plot_Graph(gnss_longitude,gnss_latitude,"GNSS Position","Longitude","Latitude")
Plot_Graph(times,gnss_velocity_n,"GNSS North Velocity","Time","Velocity")
Plot_Graph(times,gnss_velocity_e,"GNSS East Velocity","Time","Velocity")

% Write to table and export
% outputTable = table(outputs);
% writetable(outputTable,outputFile+".csv",'WriteVariableNames',0)
