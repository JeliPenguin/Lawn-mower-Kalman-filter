function Plot_Graph(longitude,latitude,v_N,v_E,times,plotType,show)
% Function for plotting figures of vehicle position and velocities and
% saving them into Figures folder
% Inputs:
%   longitude             Array of longitude solutions of the vehicle
%   latitude              Array of longitude solutions of the vehicle
%   v_N                   Array of north velocity solutions of the vehicle
%   v_E                   Array of east velocity solutions of the vehicle
%   times                 Array of time value for each epoch
%   plotType              String indicating the type of the solution (GNSS/DR)
%   show                  Boolean indicating to show the plots or not

% Setting visbility of figures
visible = "off";
if show
    visible = "on";
end

f = figure("Visible",visible);
plot(longitude, latitude);
title("Position");
xlabel("Longitude (deg)");
ylabel("Latitude (deg)");
grid on;
saveas(f,"Figures/"+plotType+"/"+plotType+"_position.png")

f = figure("Visible",visible);
plot(times, v_N);
title("North Velocity");
xlabel("Time (s)");
ylabel("Velocity (m/s)");
grid on;
saveas(f,"Figures/"+plotType+"/"+plotType+"_velocity_north.png")

f = figure("Visible",visible);
plot(times, v_E);
title("East Velocity");
xlabel("Time (s)");
ylabel("Velocity (m/s)");
grid on;
saveas(f,"Figures/"+plotType+"/"+plotType+"_velocity_east.png")

f = figure("Visible",visible);
geoplot(latitude,longitude);
title("Geography Position");
grid on;
geobasemap streets
saveas(f,"Figures/"+plotType+"/"+plotType+"_geoplot.png")

end