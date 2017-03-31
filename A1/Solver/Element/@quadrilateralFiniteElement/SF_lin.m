function [ val, dval ] = SF_lin( xi, eta , nf)
% this function computes the value of the shape functions and its derivatve
% on the domain [0,1].
%       Usage:
%       [ val, dval ] = shape_lin( xi , nf)
% The first parameter is the transformed spatial parameter xi, the second
% is the number of the testfunction (in this case ranging from 1 to 2).
%   xi  (to the right)
%   eta (upwards)
%   
%       (-1,1)  4---3  (1,1)
%               |   |
%       (-1,-1) 1---2  (1,-1)
%
%
%


if(nf == 1)
    val = 0.25 * (1-xi) * (1-eta);
    dval(1) = -0.25 * (1-eta);
    dval(2) = -0.25 * (1-xi);
elseif(nf == 2)
    val = 0.25 * (1+xi) * (1-eta);
    dval(1) = 0.25 * (1-eta);
    dval(2) = -0.25 * (1+xi);
elseif(nf == 3)
    val = 0.25 * (1+xi) * (1+eta);
    dval(1) = 0.25 * (1+eta);
    dval(2) = 0.25 * (1+xi);
elseif(nf == 4)
    val = 0.25 * (1-xi) * (1+eta);
    dval(1) = -0.25 * (1+eta);
    dval(2) = 0.25 * (1-xi);
else
    error('Invalid second argument')
end


end

