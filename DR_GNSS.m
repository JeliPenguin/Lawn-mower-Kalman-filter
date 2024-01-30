function dr_gnss_solution = DR_GNSS(gnss_solution,dr_solution,times)
% Calculates the GNSS corrected Dead Reckoning solution
% Inputs:
%   gnss_solution                       
%   times                               Array of time value for each epoch
%
% Outputs:
%   dr_gnss_solution                         Calculated DR/GNSS solution
%                                               Col 1: Time (s)
%                                               Col 2: Latitude (deg)
%                                               Col 3: Longitude (deg)
%                                               Col 4: North Velocity (m/s)
%                                               Col 5: East Velocity (m/s)

Define_Constants;

epoch_num = height(times);

% Initialization

gnss_h = gnss_solution(1,4);
gnss_lat = gnss_solution(1,2) * deg_to_rad;
[R_N,R_E]= Radii_of_curvature(gnss_lat);
P_k_m_1_plus = [sigma_v^2,0,0,0
                0,sigma_v^2,0,0
                0,0,sigma_r^2/(R_N+gnss_h)^2,0
                0,0,0,sigma_r^2/((R_E+gnss_h)^2*(cos(gnss_lat))^2)];

x_caret_k_m_1_plus = zeros(4,1);

Phi_k_m_1 = [1,0,0,0
             0,1,0,0
             tau/(R_N+gnss_h),0,1,0
             0,tau/((R_E+gnss_h)*(cos(gnss_lat))),0,1];


Q_k_m_1 = [S_DR*tau,0,0.5*S_DR*tau^2/(R_N+gnss_h),0
           0,S_DR*tau,0,0.5*S_DR*tau^2/((R_E+gnss_h)*cos(gnss_lat))
           0.5*S_DR*tau^2/(R_N+gnss_h),0,S_DR*tau^3/(3*(R_N+gnss_h)^2),0
           0,0.5*S_DR*tau^2/((R_E+gnss_h)*cos(gnss_lat)),0,S_DR*tau^3/(3*(R_E+gnss_h)^2*(cos(gnss_lat)^2))];

dr_gnss_solution = [gnss_solution(1,1),gnss_solution(1,2),gnss_solution(1,3),gnss_solution(1,5),gnss_solution(1,6)];

for i = 2:epoch_num

    time = times(i);
    gnss_lat = gnss_solution(i,2) * deg_to_rad;
    gnss_long = gnss_solution(i,3) * deg_to_rad;
    gnss_h = gnss_solution(i,4);
    gnss_v_N = gnss_solution(i,5);
    gnss_v_E = gnss_solution(i,6);
    % v_D = gnss_solution(i,7);
    [R_N,R_E]= Radii_of_curvature(gnss_lat);

    x_caret_k_minus = Phi_k_m_1*x_caret_k_m_1_plus;
    
    P_k_minus = Phi_k_m_1*P_k_m_1_plus*Phi_k_m_1'+Q_k_m_1;
    
    H_k = [0,0,-1,0
           0,0,0,-1
           -1,0,0,0
           0,-1,0,0];

    R_k = [sigma_Gr^2/(R_N+gnss_h)^2,0,0,0
           0,sigma_Gr^2/((R_E+gnss_h)^2*cos(gnss_lat)^2),0,0
           0,0,sigma_Gv^2,0
           0,0,0,sigma_Gv^2];

    K_k = P_k_minus*H_k'*inv(H_k*P_k_minus*H_k'+R_k);
    
    dr_lat = dr_solution(i,2);
    dr_long = dr_solution(i,3);
    dr_v_N = dr_solution(i,4);
    dr_v_E = dr_solution(i,5);
    delta_z_m_k = [gnss_lat-dr_lat;
                   gnss_long - dr_long;
                   gnss_v_N-dr_v_N;
                   gnss_v_E-dr_v_E
                   ]-H_k*x_caret_k_minus;
    x_caret_k_m_1_plus = x_caret_k_minus + K_k*delta_z_m_k;
    P_k_m_1_plus = (eye(4) - K_k*H_k)*P_k_minus;

    correction = [dr_v_N;dr_v_E;dr_lat;dr_long]-x_caret_k_m_1_plus;
    solution = [time,correction(3:4)'*rad_to_deg,correction(1:2)'];
    dr_gnss_solution = [dr_gnss_solution;solution];

    Phi_k_m_1 = [1,0,0,0
                 0,1,0,0
                 tau/(R_N+gnss_h),0,1,0
                 0,tau/((R_E+gnss_h)*(cos(gnss_lat))),0,1];


    Q_k_m_1 = [S_DR*tau,0,0.5*S_DR*tau^2/(R_N+gnss_h),0
               0,S_DR*tau,0,0.5*S_DR*tau^2/((R_E+gnss_h)*cos(gnss_lat))
               0.5*S_DR*tau^2/(R_N+gnss_h),0,S_DR*tau^3/(3*(R_N+gnss_h)^2),0
               0,0.5*S_DR*tau^2/((R_E+gnss_h)*cos(gnss_lat)),0,S_DR*tau^3/(3*(R_E+gnss_h)^2*(cos(gnss_lat)^2))];

end

% Writing solution to csv file
outputTable = array2table(dr_gnss_solution);
writetable(outputTable,"Solutions/DR_GNSS_Solution.csv",'WriteVariableNames',0)


end