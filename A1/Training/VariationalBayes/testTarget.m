function [val,info] = testTarget(X, mu, Sigma)
% Return value of logQopt_i = log_pc_i + log_pcf_i

    

    R = chol(Sigma);
    log_det_Sigma = 2 * sum(log(diag(R)));
    v = R'\(X'-mu');
    
    val = -0.5 * log_det_Sigma - 0.5 * (v'*v);
    
    dval = - (R\v)';
    
    info{1,1} = [];
    info{1,2} = dval;

%     F = @(X) log(mvnpdf(X,mu,Sigma));
%     Check whether gradient is computed correctly
%     dFdX = computeNumericalGradient(F, X)
%     plot1DDistribution(F,X)




end



% Compute the gradient numerically and compare it to the analytical result:
function dFdX = computeNumericalGradient(F, X)
% F is a function handle. X is a input vector;

dFdX = zeros(1,length(X));
h = 1e-8;
dX = zeros(1,length(X));
for d = 1:length(X)
   dX(d) = h;
   dFdX(d) = (F(X+dX) - F(X)) /h;
   %set back
   dX(d) = 0;
end

end


function plot1DDistribution(F,Xloc)
n = 200;
X = linspace(-5,5,n);
Xfix = repmat(Xloc(1:end-1)',1,n);
X = [Xfix; X];

X_ref = X(:,1);
X_ref(end) = -5;

p=zeros(n,1);
normalization = -F(X_ref);
for i = 1:n
   p(i) = exp(F(X(:,i))+normalization);
end

figure;
plot(X(end,:),p,'x')

end










