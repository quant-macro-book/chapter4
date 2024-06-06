function err = shooting_err(c0,k0,T,alpha,beta,delta,Abar,kss)

[k,c] = shooting(c0,k0,T,alpha,beta,delta,Abar);

err = k(T)-kss;