clear all;

xL = 0.1;
xH = 2.0;
fL = f(xL);
fH = f(xH);
diff = 1e+4;
iter = 1;

while (diff>1e-10)
    
    disp([iter xL xH fL fH]);
    
    x0 = (xL+xH)/2;
    f0 = f(x0);
    
    if (fL*f0<0)
        xH = x0;
        fH = f0;
    else % fH*f0<0
        xL = x0;
        fL = f0;
    end
    
    diff = abs(xH-xL);
    iter = iter+1;
    
end