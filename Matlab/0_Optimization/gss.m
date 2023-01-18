clear all;

%x = [-10:0.01:10];
%plot(x,f(x));
phi = (-1+sqrt(5))/2; % 1/phi = golden ratio

xL = 0.0;
xH = 6.0;
fL = f_gss(xL);
fH = f_gss(xH);
x1 = xL + phi*(xH-xL);
f1 = f_gss(x1);
x2 = xH - phi*(xH-xL);
f2 = f_gss(x2);
disp([xL x2 x1 xH])
disp([fL f2 f1 fH])
pause;

diff = 1e+4;
iter = 1;

while (diff>1e-10)
    
    disp([iter xL x2 x1 xH fL fH]);
    
    if (f2>f1)
        xH = x1;
        fH = f1;
        x1 = x2;
        f1 = f2;
        x2 = xH - phi*(xH-xL);
        f2 = f_gss(x2);
    else
        xL = x2;
        fL = f2;
        x2 = x1;
        f2 = f1;
        x1 = xL + phi*(xH-xL);
        f1 = f_gss(x1);
    end
    
    diff = abs(xH-xL);
    iter = iter+1;
    
end