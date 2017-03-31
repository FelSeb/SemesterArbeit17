function [ features ] = featfun2( x_tensor , c2f)
% x has Elem2param format

% Each feature must be returned in vector format -NOT in Elem2param format! 
% The corresponding conversion for an arbitrary X_tensor in Elem2param 
% format is done as follows: X = reshape(X_tensor,[numel(X_tensor), 1]);
% Every column of the features matrix contains one feature.

% How to use Constantin's feature functions:
% phi = someFeatFun(lambda)
% phi the value of the effective parameter (or feature) of one coarse
% element. 
% lambda is a vector containing the parameter values of each fine grained
% element within one coarse element.

% The following steps are taken to ensure the generality of the featfun2
% function. It may be thus applied to parameter arrays of the following 
% formats:
% X_ = [X_11]
% X_ = [X_11; X_22]
% X_ = [X_11, X_12; X_21, X_22] 


ndimdata = ndims(x_tensor);
firstdims = repmat({':'},1,ndims(x_tensor)-1);

N_features = 1;

N_coarse_elements = numel(c2f(:));
x_tensor_size = size(x_tensor(firstdims{:},1));
features_size = [x_tensor_size, N_features];
features = zeros(features_size);

% Feature 1: local mean
for ce = 1:N_coarse_elements
    features(firstdims{:},ce,1) = mean(x_tensor(firstdims{:},c2f{ce}),ndimdata);
end

% Feature 2: lineal path (there are multiple lineal paths for multiple directions)
for ce = 1:N_coarse_elements
    features(firstdims{:},ce,2) = mean(x_tensor(firstdims{:},c2f{ce}),ndimdata);
end



% Feature 3: mean image properties


% Convert to vector format
features = reshape(features,[numel(x_tensor(firstdims{:},1))*N_coarse_elements, N_features]);


end






