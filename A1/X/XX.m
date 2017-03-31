
clear;
clc;
close all;

% logpdf = @(x) log(normpdf(x,1,1));
logpdf = @(x) alternativelogpdf(x);



% 1 Chain
% logproppdf = @(x,y) log(normpdf(x,y,1));
% proprnd = @(y) normrnd(y,1);

% 2 Chains
propvar = 2;
logproppdf = @(x,y) [log(normpdf(x(1),y(1),propvar)); log(normpdf(x(2),y(2),propvar));];
proprnd = @(y) [normrnd(y(1),propvar); normrnd(y(2),propvar)];


N_samples = 10;
% [smpl, accept] = MYmhsample([0;0],N_samples,'logpdf',logpdf,'logproppdf',logproppdf,'proprnd',proprnd,'returninfo',0,'nchain',2);
[smpl, accept] = MYmhsample([0;0],N_samples,'logpdf',logpdf,'logproppdf',logproppdf,'proprnd',proprnd,'returninfo',1,'nchain',2);


accept