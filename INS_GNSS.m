function ins_gnss_solution = INS_GNSS()
% Calculates the integrated solution of GNSS and Inertial Navigation System
% Inputs:
%   dr_measurement_data                 Array of raw dead reckoning measurement data
%   times                               Array of time value for each epoch
%
% Outputs:
%   ins_gnss_solution                         Calculated INS/GNSS solution
%                                               Col 1: Time (s)
%                                               Col 2: Latitude (deg)
%                                               Col 3: Longitude (deg)
%                                               Col 4: North Velocity (m/s)
%                                               Col 5: East Velocity (m/s)

Define_Constants;




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
end
