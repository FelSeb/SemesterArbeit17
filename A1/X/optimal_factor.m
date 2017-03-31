lo = 0.5;
hi = 2;
opt = 0.4;

c = lo;
b = ((1-lo) -(hi-lo)*opt^2)/(opt-opt^2);
a = hi-lo-b;


x = linspace(-0.1,1.1,30);

fac = a*x.^2 + b*x + c;

plot(x,fac)