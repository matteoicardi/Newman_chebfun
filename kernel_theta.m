K=600;
theta=1;
x=chebfun('x',[pi,(K+1)*pi]);
f=theta*sin(x)-x.*cos(x);
alpha=roots(f)/pi;

L=100;
fs=1/1000;
N=L/fs+1;
x=linspace(1,N,N).*fs*10;
%close all
f=exp(x*0)-1;
ft=f;
for k=1:K
    f=f+exp(-alpha(k)^2*x);
%    ft=ft+2*alpha(k)^2./(x.^2+alpha(k)^4);   % fourier
    ft=ft+1./(x+alpha(k)^2);   % laplace
    
%     figure(1)
%     hold on
%     loglog(x,f)
%     figure(2)
%     hold on
%     %ff=fft(f);
%     %loglog(x(1:(N-1)/2),abs(ff(1:(N-1)/2)))
%     loglog(ft)
%     figure(3)
%     hold on
%     loglog(ft.*sqrt(1+x))
%     %loglog(x(1:(N-1)/2),abs(ff(1:(N-1)/2)).*sqrt(1+1e1*x(1:(N-1)/2)))
end

%  figure()
%  hold on
%  loglog(x,f)
 figure(6)
 hold on
% ff=fft(f);
% loglog(x(1:(N-1)/2),abs(ff(1:(N-1)/2)))
 loglog(x,ft)
 figure(7)
 hold on
 loglog(x,ft.*sqrt(1+x))
% loglog(x(1:(N-1)/2),abs(ff(1:(N-1)/2)).*sqrt(1+1e1*x(1:(N-1)/2)))