function dr_solution = Dead_Reckoning_modified(dr_measurement_data,gnss_solution,heading_solution,times)
% Calculates Dead Reckoning solution
% Inputs:
%   dr_measurement_data                 Array of raw dead reckoning measurement data
%   gnss_solution                       Array gnss kalman filter solutions
%   heading_solution                    Array of gyro-smooth magnetic heading solutions
%   times                               Array of time value for each epoch
%
% Outputs:
%   dr_solution                         Calculated DR solution, for each row
%                                           Col 1: Time (s)
%                                           Col 2: Latitude (deg)
%                                           Col 3: Longitude (deg)
%                                           Col 4: Damped North Velocity (m/s)
%                                           Col 5: Damped East Velocity (m/s)

Define_Constants;

dr_measurement = dr_measurement_data(:,2:end);

% Total count of epochs
epoch_num = height(times); 

lat = gnss_solution(1,2) * deg_to_rad;
long = gnss_solution(1,3)*deg_to_rad;
h = gnss_solution(1,4);

% Array for storing calculated dead reckoning solutions
dr_solution = zeros(epoch_num,5);
delta_r_N_prev = zeros(epoch_num, 1);
delta_r_E_prev = zeros(epoch_num, 1);


for i=1:epoch_num
    time = times(i);
    
    % Caulates average speed of vehicle given wheel measurements
    wheel_speed_measurements = dr_measurement(i,1:4);
    gyro_measurement = dr_measurement(i,5);
    compass_measurement = dr_measurement(i, 6);
    

    heading_k = heading_solution(i) * deg_to_rad;

    if i==1
        % Assume the instantaneous velocity at first epoch is given by the
        % speed and heading measurements at that time
        v_k = Calc_Avg_Speed(wheel_speed_measurements);
        v_N_k = v_k*cos(heading_k);
        v_E_k = v_k*sin(heading_k);
    else
        [R_N,R_E]= Radii_of_curvature(lat);

        % Heading from previous epoch
        heading_k_m_1 = dr_measurement(i-1,6) * deg_to_rad;

        delta_r = Calc_Change_Position(wheel_speed_measurements,gyro_measurement, compass_measurement);
        delta_r_N = delta_r(1);
        delta_r_E = delta_r(2);

        % Store the delta_r values
        delta_r_N_prev(i) = delta_r_N;
        delta_r_E_prev(i) = delta_r_E;

        % Longitude and latitude calculation
        lat = lat + delta_r_N/(R_N+h);
        long = long + delta_r_E/((R_E+h)*cos(lat));
        
        % Calculate to save velocities
        v_N_k = delta_r_N / tau;
        v_E_k = delta_r_E / tau;
    end
    
    dr_solution(i,:)=[time,lat*rad_to_deg,long*rad_to_deg,v_N_k,v_E_k];
end

solutionTable = array2table(dr_solution);
writetable(solutionTable,"Solutions/DR_Solution.csv",'WriteVariableNames',0)

end



