function [val,info] = logQopt_i(this,X,x_i,y_i, varargin)
% Return value of logQopt_i = log_pc_i + log_pcf_i

if(numel(varargin)>0)
    compute_grad = varargin{1};
else
    compute_grad = 0;
end

if(~compute_grad)
    val_pc = this.PC.evalLogPDF(X, x_i);
    [val_pcf, y] = this.PCF.evalLogPDF(y_i, X);
    val = val_pc + val_pcf;
    
    info{1} = y;
else
    [val_pc, dval_pc] = this.PC.evalLogPDF(X, x_i,1);
    [val_pcf, y, dval_pcf] = this.PCF.evalLogPDF(y_i, X,1);

%     % Compare the fine solution to the interpolated coarse solution
%     PostSol.T = y_i;
%     this.PC.Fine_grid.plotSol(PostSol, 0);
%     
%     PostSol.T = this.PCF.Interpolation_matrix*y;
%     this.PC.Fine_grid.plotSol(PostSol, 0);
    
    
    val = val_pc + val_pcf;
    
    dval = dval_pc + dval_pcf;
    
%     dval_norm = norm(dval)
    
    info{1,1} = y;
    info{1,2} = dval;
    
%     Illustrate Pc and Pcf:
%        F = @(XX) this.PC.evalLogPDF(XX, x_i,1) + this.PCF.evalLogPDF(y_i, XX,1);
%      plot1DDistribution(F,X)
%     F1 = @(XX) this.PC.evalLogPDF(XX, x_i,1);
%     F2 = @(XX) this.PCF.evalLogPDF(y_i, XX,1);
%     plot1DDistribution(F1,X)
%     plot1DDistribution(F2,X)
%     Check whether gradient is computed correctly
%      dFdX = computeNumericalGradient(F, X)
%      plot1DDistribution(F,X)

    
end


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










