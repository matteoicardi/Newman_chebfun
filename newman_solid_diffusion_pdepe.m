function res=newman_solid_diffusion_pdepe(t,mu)
m = 2;  % m=2 spherical
omega=2*pi;

x = linspace(0,1,1000);
%t = linspace(0,T*omega,T)*2*pi;

sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);
% Extract the first solution component as u.
uu = sol(:,:,1);

res = uu(:,end);
% % A surface plot is often a good way to study a solution.
% figure(1)
% hold on
% mesh(x,t,uu) 
% xlabel('Distance x')
% ylabel('Time t')
% 
% figure(2)
% plot(x,uu(end,:))
% xlabel('Distance x')
% ylabel('u(x,T)')
% 
% figure(3)
% plot(t,uu(:,end))
% xlabel('Time t')
% ylabel('u(R,t)')



% --------------------------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,DuDx)
Ds = mu;  %diffusion coefficient
c = 1;
f = DuDx*Ds;
s = 0;
end
% --------------------------------------------------------------
function u0 = pdex1ic(x)
u0 = 0;
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = 0;
ql = 0;
pr = 1;%-sinc(omega*t/pi);
qr = 1;
%pl = -cos(omega*t);
%ql = 1;
end

end