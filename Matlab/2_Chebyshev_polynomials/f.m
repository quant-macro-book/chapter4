function f0 = f(x0)

    nx = size(x0,1);
    f0 = ones(nx,1)./(ones(nx,1)+25*x0.^2);