function k = polygrid(kmin,kmax,N)

x = -cos(pi*[0:N-1]'/(N-1)); % チェビシェフ極値点
%x = [0;-cos((2*[1:N-1]'-1)*pi/2/(N-1))]; % チェビシェフゼロ点

% x\in[-1,1]からkに変換
k = (x+1)*(kmax-kmin)/2 + kmin.*ones(N,1);