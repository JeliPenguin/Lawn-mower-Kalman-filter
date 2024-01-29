% Configure outputs
outputFile = "Output_Profile";
% 1: Time
% 2: Geodetic latitude in degrees
% 3: Geodetic longitude in degrees 
% 4: North velocity in m/s
% 5: East velocity in m/s
% 6: Heading in degrees.
% outputs = zeros(epochNum,6);

% GNSS Kalman Filter estimation
GNSS_KF;

% Write to table and export
% outputTable = table(outputs);
% writetable(outputTable,outputFile+".csv",'WriteVariableNames',0)
