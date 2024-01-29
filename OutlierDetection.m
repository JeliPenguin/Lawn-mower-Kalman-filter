function hasOutlier = OutlierDetection(numberOfSatellites,H_G_e,delta_z_min,sigma_p,T)
    % Outlier detection
    hasOutlier = false;
    v = (H_G_e*inv(H_G_e.'*H_G_e)*H_G_e.'-eye(numberOfSatellites))*delta_z_min; % residuals vector
    C_v = (eye(numberOfSatellites)-H_G_e*inv(H_G_e.'*H_G_e)*H_G_e.')*sigma_p^2;
    maxJ = -1;
    maxResidual = -inf;
    for j=1:numberOfSatellites
        normalized_residual = norm(v(j))/sqrt(C_v(j,j));
        if normalized_residual > T
            % disp(["Outlier detected at: ",j," Residual: ",normalized_residual])
            % sprintf("Outlier detected with Satellite %d, Residual: %f",satelliteNumbers(j),normalized_residual)
            if normalized_residual > maxResidual
                maxJ = j;
                maxResidual = normalized_residual;
            end
            hasOutlier = true;
        end
    end
end