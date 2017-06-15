function out=newman_solid_diffusion_random_quad(t,mu,l,v)
%load randomcube.mat
%sim=mu03;
%D = 1e-4;
amp = 1;
omega=2*pi;
%R=1; %0.398288;
V=4*pi/3; %1.03543;

tolerance=1e-2;

fact= -3;%/R;

t0=-log(tolerance)/mu*max(l)

    function kk=kernel(tau)
        kk=ones(size(tau));
        for k=1:length(l)
            kk=kk-V*v(k)*exp(-mu*tau/l(k));
        end
        %kk(isinf(kk)) = 0;
        kk=kk.*(tau>0);
    end

myfun = @(tt) fact*kernel(t-tt)*cos(omega*tt)*amp;

out=integral(myfun,-t0,max(t), 'RelTol',tolerance, 'ArrayValued', true);

end

