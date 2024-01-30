function dr_solution = Dead_Reckoning(dr_measurement_data,gnss_solution,times)
% Calculates Dead Reckoning solution
% Inputs:
%   dr_measurement_data                 Array of raw dead reckoning measurement data
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

dr_solution = [];

for i=1:epoch_num
    time = times(i);
    wheel_speed_measurements = dr_measurement(i,1:4);
    v_k = Calc_Avg_Speed(wheel_speed_measurements);
    heading_k = dr_measurement(i,6) * deg_to_rad;
    if i==1
        v_N_k = v_k*cos(heading_k);
        v_E_k = v_k*sin(heading_k);
    else
        [R_N,R_E]= Radii_of_curvature(lat);
        heading_k_m_1 = dr_measurement(i-1,6) * deg_to_rad;
        v_vec = 0.5*v_k*[cos(heading_k)+cos(heading_k_m_1);sin(heading_k)+sin(heading_k_m_1)];
        v_cap_N_k = v_vec(1);
        v_cap_E_k = v_vec(2);
        v_N_k=1.7*v_cap_N_k-0.7*v_N_k;
        v_E_k=1.7*v_cap_E_k-0.7*v_E_k;

        % Longitude and latitude calculation
        prevTime = times(i-1);
        lat = lat + v_cap_N_k*(time - prevTime)/(R_N+h);
        long = long + v_cap_E_k*(time - prevTime)/((R_E+h)*cos(lat));
    end
    dr_solution=[dr_solution;[time,lat*rad_to_deg,long*rad_to_deg,v_N_k,v_E_k]];
end

solutionTable = array2table(dr_solution);
writetable(solutionTable,"Solutions/DR_Solution.csv",'WriteVariableNames',0)


end



