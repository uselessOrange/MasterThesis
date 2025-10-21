function [dataWithoutNans,nan_indices]  = GetRidOfNans(data)
%Function finds NaN value entries in data and removes them
%data should be a nxm matrix with n samples of m features 
%( i.e column vercors )
% Transpose data if there are more columns than rows
[n, m] = size(data);
if n < m
    data = data';
    [n, m] = size(data);
end

% Find indices of NaN values
nan_indices = any(isnan(data), 2);

% Select rows without NaN values using logical indexing
dataWithoutNans = data(~nan_indices, :);
end
