function [ val, dval ] = SF_lin1D( xi , nf)
% this function computes the value of the shape functions and its derivatve
% on the domain [0,1].
%       Usage:
%       [ val, dval ] = shape_lin( xi , nf)
% The first parameter is the transformed spatial parameter xi, the second
% is the number of the testfunction (in this case ranging from 1 to 2).

if(nf == 1)
    val = -0.5*xi + 0.5;
    dval = -0.5;
elseif(nf == 2)
    val = 0.5*xi + 0.5;
    dval = 0.5;
else
    error('Invalid second argument')
end


end