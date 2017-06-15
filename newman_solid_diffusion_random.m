function out=newman_solid_diffusion_random(tt,mu,l,v)
%load randomcube.mat
%sim=mu03;
%D = 1e-4;
amp = 1;
omega=2*pi;
%R=1; %0.398288;
V=4*pi/3; %1.03543;

tolerance=1e-4;

fact= 3;%/R;

t=chebfun('t',[0,tt(end)]);
t0=-log(tolerance)/mu*max(l)
tm=chebfun('t',[-t0,tt(end)]);
tp=chebfun(@(t) exp(-mu*t) ,[0,t0]);
%mu=D/(R)^2;

splitting off
KK=tp*0+1;
for k=1:length(l)
    KK=KK+V*v(k)*tp.^(1/l(k));
end

resp=conv(cos(omega*tm)*amp,fact*KK,'same');

out=resp(tt);