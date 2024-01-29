Define_Constants;

speedHeadingTable = readtable("Workshop3_Speed_Heading.csv");
speedHeading = table2array(speedHeadingTable);

epochNum =height(speedHeading);

h = 37.4;
lat = 50.424958 * deg_to_rad;
long = -3.5957974*deg_to_rad;
v_N_k = 0;
v_E_k = 0;
v_cap_N_k = 0;
v_cap_E_k = 0;
v_cap_k = 0;

solution = [];

for i=1:epochNum
    time = speedHeading(i,1);
    v_k = speedHeading(i,2);
    heading_k = speedHeading(i,3) * deg_to_rad;
    if i==1
        v_N_k = v_k*cos(heading_k);
        v_E_k = v_k*sin(heading_k);
    else
        [R_N,R_E]= Radii_of_curvature(lat);
        heading_k_m_1 = speedHeading(i-1,3) * deg_to_rad;
        v_vec = 0.5*v_k*[cos(heading_k)+cos(heading_k_m_1);sin(heading_k)+sin(heading_k_m_1)];
        v_cap_N_k = v_vec(1);
        v_cap_E_k = v_vec(2);
        v_N_k=1.7*v_cap_N_k-0.7*v_N_k;
        v_E_k=1.7*v_cap_E_k-0.7*v_E_k;
        % Longitude and latitude calculation
        prevTime = speedHeading(i-1,1);
        lat = lat + v_cap_N_k*(time - prevTime)/(R_N+h);
        long = long + v_cap_E_k*(time - prevTime)/((R_E+h)*cos(lat));
    end
    solution=[solution;[time,lat,long,v_N_k,v_E_k]];
end

solutionTable = array2table(solution);
writetable(solutionTable,"task1Solution.csv",'WriteVariableNames',0)




