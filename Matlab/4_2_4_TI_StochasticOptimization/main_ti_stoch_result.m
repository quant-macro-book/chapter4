clear all;

load stoch_result.mat;

pfcn = zeros(nk,nz);

for j = 1:11
    
    pfcn(:,j) = exp(zgrid(j))*kgrid.^alpha + (1-delta)*kgrid - cfcn0(:,j);

end

figure;
plot(kgrid,pfcn(:,1));
hold on;
plot(kgrid,pfcn(:,6));
plot(kgrid,pfcn(:,11));
xlim([kgrid(1) kgrid(11)]);

figure;
subplot(311);
plot(cvec);
ylabel('c');
subplot(312);
plot(kvec(1:T));
ylabel('k');
subplot(313);
plot(exp(zgrid(ivec(1:T))));
ylabel('A');
xlabel('time');

print -depsc2 TI_stochsim.eps
