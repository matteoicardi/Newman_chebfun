function y=f(t,y0,param)

m=param{1};
c0=param{2};
D=param{3};
j=param{4};
dj=param{5};
Ds=param{6};

% hold off
% plot(sum(-Ds*D*y0),'r')
% hold on
% plot(sum(-chebmatrix(c0')*chebmatrix(dj)),'b')
% plot(sum(Ds*chebmatrix(m')*chebmatrix(j)),'k')
y = -Ds*D*y0 - chebmatrix(c0')*chebmatrix(dj)/Ds ...
    + chebmatrix(m')*chebmatrix(j);

