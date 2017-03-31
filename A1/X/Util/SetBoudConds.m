% Temperature field and gradient generating the boundary conditions
boundaryCoeffs = [0 0 0 1];

BoundVals = @(x,y,boundaryCoeffs) boundaryCoeffs(1) + boundaryCoeffs(2)*x + boundaryCoeffs(3)*y + boundaryCoeffs(4)*x*y;

BoundGrads = @(x,y,boundaryCoeffs) [boundaryCoeffs(2)+ boundaryCoeffs(4)*y; boundaryCoeffs(3) + boundaryCoeffs(4)*x];