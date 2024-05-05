function Plot_heading(dr_gnss_longitude,dr_gnss_latitude,heading_solution,visible)
% Function for plotting heading of the vehicle using quiver plots and
% saving them into Figures folder
% Inputs:
%   dr_gnss_longitude             Array of integrated longitude solutions of the vehicle
%   dr_gnss_latitude              Array of integrated longitude solutions of the vehicle
%   heading_solution              Array of heading solutions
%   show                          Boolean indicating to show the plots or not

Define_Constants;
f = figure("Visible",visible);
epoch_num = height(heading_solution);
x = zeros(epoch_num,1);
y = zeros(epoch_num,1);

% Magnitude of each quiver vector
mag = 0.0000000001;

[y(1),x(1)] = pol2cart(heading_solution(1)*deg_to_rad,mag);

for i = 2: epoch_num
    % mag = sqrt((dr_gnss_latitude(i)-dr_gnss_latitude(i-1))^2+(dr_gnss_longitude(i)-dr_gnss_longitude(i-1))^2);

    % Compute vector from given heading 
    [yVec,xVec] = pol2cart(heading_solution(i)*deg_to_rad,mag);
    x(i) = xVec;
    y(i) = yVec;
end

quiver(dr_gnss_longitude,dr_gnss_latitude,x,y);
title("Heading");
xlabel("Longitude");
ylabel("Latitude");
grid on;
saveas(f,"Figures/Heading/heading.png")
end