function [ features ] = featfun( x_tensor , c2f)
% x hast Elem2param format

% Each feature must be returned in vector format NOT in Elem2param format. 
% The corresponding conversion for an arbitrary X_tensor in Elem2param 
% format is done as follows: X = reshape(X_tensor,[numel(X_tensor), 1]);
% Every column of the features matrix contains one such feature.

N_features = 1;

N_coarse_elements = numel(c2f);
x_tensor_size = size(x_tensor);
features_size = [x_tensor_size, N_features];
features = zeros(features_size);

ndimdata = ndims(data);
firstdims = repmat({':'},1,ndims(data)-1);

% Feature 1: local mean
for ce = 1:N_coarse_elements
    features(firstdims{:},ce,1) = mean(x_tensor(firstdims{:},c2f{ce}),ndimdata);
end


% Convert to vector format
features = reshape(features,[numel(x_tensor(firstdims{:},1))*N_coarse_elements, N_features]);


end






