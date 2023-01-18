function [f,df] = f(x)

f = -x.^4 + x + 3;
df = 4*x.^3 + 1;