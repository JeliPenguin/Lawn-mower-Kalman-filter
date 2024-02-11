function [valid_indicies, maxJ] = Outlier_Detection(number_of_satellites,H_G_e,delta_z_min,sigma)
    % Outlier detection
    Define_Constants;

    % Compute residuals vector
    v = (H_G_e * inv(H_G_e.' * H_G_e) * H_G_e.' - eye(number_of_satellites)) * delta_z_min; 

    % Compute residuals covariance matrix
    C_v = (eye(number_of_satellites) - H_G_e * inv(H_G_e.' * H_G_e) * H_G_e.') * sigma^2;

    % Setup
    valid_indicies = zeros(1,number_of_satellites);
    maxResidual = -inf;
    maxJ = 0;

    % Check for outliers
    for j=1:number_of_satellites
        normalized_residual = norm(v(j))/sqrt(C_v(j,j));
        if normalized_residual > T
            % sprintf("Outlier detected with Satellite %d, Residual: %f",j,normalized_residual)
            if normalized_residual > maxResidual
                maxJ = j;
                maxResidual = normalized_residual;
            end
        else
            valid_indicies(j) = 1;
        end
    end
end

