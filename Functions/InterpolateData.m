function data_interp = InterpolateData(data,threshold,x)
% data_interp = InterpolateData(data,threshold,x)
%threshlod = [lower,upper];
%Interpolates data by remooving samples below a threshold and replacing
%them with a linear function between two 
error_indices = find(data >= threshold(1) & data <= threshold(2));
non_error_indices = setdiff(1:numel(data), error_indices);

% Linear interpolation needs two points to interpolate between. If the last
% sample gets under the threshold the code wouldn't work
if data(end)<= threshold
    data(end)=data(non_error_indices(end));
end

data_interp  = interp1(x(non_error_indices), data(non_error_indices), x, 'linear');
end