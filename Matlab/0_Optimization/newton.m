clear all;

x0 = 10;
diff = 1e+4;
iter = 1;

while (diff>1e-10)
    
    disp([iter x0]);

    [f0,df0] = f(x0);
    x1 = x0 - f0/df0;
    diff = abs(x1-x0);
    iter = iter+1;
    x0 = x1;
    
end