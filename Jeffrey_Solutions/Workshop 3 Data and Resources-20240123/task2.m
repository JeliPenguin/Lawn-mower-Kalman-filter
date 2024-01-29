clear;
Define_Constants;

gnssPosVelTable = readtable("Workshop3_GNSS_Pos_Vel_NED.csv");
task1SolutionTable = readtable("task1Solution.csv");
task1Solution = table2array(task1SolutionTable);
gnssPosVel = table2array(gnssPosVelTable);
epochNum = height(gnssPosVel);

% Initialization
tau_s = 0.5;
sigma_v = 0.1;
sigma_r = 10;
sigma_Gr = 5;
sigma_Gv = 0.02;
S_DR = 0.2;
h = gnssPosVel(1,4);
lat = gnssPosVel(1,2) * deg_to_rad;
[R_N,R_E]= Radii_of_curvature(lat);
P_plus_k_m_1 = [sigma_v^2,0,0,0
                0,sigma_v^2,0,0
                0,0,sigma_r^2/(R_N+h)^2,0
                0,0,0,sigma_r^2/((R_E+h)^2*(cos(lat))^2)];

x_caret_plus_k_m_1 = zeros(4,1);

Phi_k_m_1 = [1,0,0,0
             0,1,0,0
             tau_s/(R_N+h),0,1,0
             0,tau_s/((R_E+h)*(cos(lat))),0,1];


Q_k_m_1 = [S_DR*tau_s,0,0.5*S_DR*tau_s^2/(R_N+h),0
           0,S_DR*tau_s,0,0.5*S_DR*tau_s^2/((R_E+h)*cos(lat))
           0.5*S_DR*tau_s^2/(R_N+h),0,S_DR*tau_s^3/(3*(R_N+h)^2),0
           0,0.5*S_DR*tau_s^2/((R_E+h)*cos(lat)),0,S_DR*tau_s^3/(3*(R_E+h)^2*(cos(lat)^2))];

solutionTable = [gnssPosVel(1,1),gnssPosVel(1,2)*deg_to_rad,gnssPosVel(1,3)*deg_to_rad,gnssPosVel(1,4),gnssPosVel(1,5),gnssPosVel(1,6),gnssPosVel(1,7)];

for i = 2:epochNum

    time = gnssPosVel(i,1);
    lat = gnssPosVel(i,2) * deg_to_rad;
    long = gnssPosVel(i,3) * deg_to_rad;
    h = gnssPosVel(i,4);
    v_N = gnssPosVel(i,5);
    v_E = gnssPosVel(i,6);
    v_D = gnssPosVel(i,7);
    [R_N,R_E]= Radii_of_curvature(lat);

    x_caret_m_k = Phi_k_m_1*x_caret_plus_k_m_1;
    
    P_minus_k = Phi_k_m_1*P_plus_k_m_1*Phi_k_m_1'+Q_k_m_1;
    
    H_k = [0,0,-1,0
           0,0,0,-1
           -1,0,0,0
           0,-1,0,0];

    R_k = [sigma_Gr^2/(R_N+h)^2,0,0,0
           0,sigma_Gr^2/((R_E+h)^2*cos(lat)^2),0,0
           0,0,sigma_Gv^2,0
           0,0,0,sigma_Gv^2];

    K_k = P_minus_k*H_k'*inv(H_k*P_minus_k*H_k'+R_k);
    
    latDR = task1Solution(i,2);
    longDR = task1Solution(i,3);
    v_N_DR = task1Solution(i,4);
    v_E_DR = task1Solution(i,5);
    delta_z_m_k = [lat-latDR;
                   long - longDR;
                   v_N-v_N_DR;
                   v_E-v_E_DR
                   ]-H_k*x_caret_m_k;
    x_caret_plus_k_m_1 = x_caret_m_k + K_k*delta_z_m_k;
    P_plus_k_m_1 = (eye(4) - K_k*H_k)*P_minus_k;

    solution = [v_N_DR;v_E_DR;latDR;longDR]-x_caret_plus_k_m_1;
    solution = [time,solution(3:4)',gnssPosVel(1,4),solution(1:2)',gnssPosVel(1,7)];
    solutionTable = [solutionTable;solution];

    Phi_k_m_1 = [1,0,0,0
                 0,1,0,0
                 tau_s/(R_N+h),0,1,0
                 0,tau_s/((R_E+h)*(cos(lat))),0,1];


    Q_k_m_1 = [S_DR*tau_s,0,0.5*S_DR*tau_s^2/(R_N+h),0
               0,S_DR*tau_s,0,0.5*S_DR*tau_s^2/((R_E+h)*cos(lat))
               0.5*S_DR*tau_s^2/(R_N+h),0,S_DR*tau_s^3/(3*(R_N+h)^2),0
               0,0.5*S_DR*tau_s^2/((R_E+h)*cos(lat)),0,S_DR*tau_s^3/(3*(R_E+h)^2*(cos(lat)^2))];

end
% disp(solutionTable)
solutionTable = array2table(solutionTable);
writetable(solutionTable,"task2Solution.csv",'WriteVariableNames',0)
