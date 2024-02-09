function [corrected_data] = Outlier_Data_Correction(valid_indicies,data)
% Corrects inputting data according to an array of valid indicies, where
% index with 1 means the data at that index is valid and vice versa
% Inputs:
%   valid_indicies          Valid indicies array
%   data                    Data to be corrected
%
% Outputs:
%   corrected_data         Data corrected using the valid indicies array
corrected_data = (nonzeros(valid_indicies.*data))';
end