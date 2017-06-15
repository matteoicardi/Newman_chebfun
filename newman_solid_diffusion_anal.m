%% Analytical solution of the solid diffusion for periodic forcing
function f=newman_solid_diffusion_anal(t,mu,alpha)

%K = 100;  % number of modes

%x=chebfun('x',[pi,(K+1)*pi]);
%f=sin(x)-x.*cos(x);
%alpha=roots(f);

%% parameters
omega = 2*pi;  % frequency I=cos(omega*t)
%mu = 1e-5/4;  % Ds/R^2
amp = 1;

% beta=(alpha.^2*mu/omega)./sqrt(1+((mu*alpha.^2)./omega).^2);
% amp=amp*beta./(alpha.^4)*omega/mu^2*2;
% 
% %% function
% f = 0*t;
% f = 3/omega*sin(omega*t) + 1/(5*mu)*cos(omega*t);
% 
% for i=1:length(alpha)
%     f = f - amp(i)*cos(omega*t+asin(beta(i)));
% end

%% New formula

beta=1./sqrt(1+((mu*alpha.^2)./omega).^2);
amp=2./sqrt(omega^2+alpha.^4*mu^2);

%% function
f = 0*t;
f = 3/omega*sin(omega*t);

for i=1:length(alpha)
    f = f + amp(i)*sin(omega*t+asin(beta(i)));
end



end
