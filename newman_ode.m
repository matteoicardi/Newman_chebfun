function res=newman_ode(mu,param)
% Function to be used for solving a ODE using chebfun
% not used now as I am still not sure which ODE I should solve
% mu is the relaxation time of the particle (Radius^2/Diffusivity)

T=10/mu;  % Total time
tf=chebfun('t',[0,T]);
xf=sinc(pi*tf*10*mu);  % input function and its derivatives
xf1=diff(xf);
xf2=diff(xf,2);
N = chebop(0,T);
N.lbc = @(x,y) [x; y];
figure()
plot(xf)
hold on

% param = [mu*p0 p1 p2/mu mu^2*q0]

res=solve([mu*param(1) param(2) param(3)/mu mu^2*param(4)]);

    function y=solve(p)
    N.op = @(t,y,z) [diff(y) - z + p(2)*p(4)*y/p(1);
                     diff(z) + p(4)*y - p(1)*xf - p(2)*xf1 - p(3)*xf2 ];
    [y,z] = N\0;


    end

plot(res)
end