function y=f(t,y0,param)

M=param(1);
D=param(2);
R=param(3);
j=param(4);

y = ones(M,1);

yy=y0;
% yy(3:end) = y0(3:end)-2*R*y(2:end-1)+R^2*y(1:end-2);
% yy(2) = y0(2)-R*y(1);
% yy(1) = y0(1);

global abl abml
[ab,abm]=chebyshev(floor(M/2),yy);
ab=real(ab);
abm=abs(abm);
ml=floor(M/2);
while ml>1
    %q=real(radau(ml-1,ab,R));
    %q=real(lobatto(ml-2,ab,0,R));
    q=real(gauss(ml,ab));
    if min(q(:,1))<0 || max(q(:,1))<R
        ml=ml-1;
    else
        break
    end
end
disp(ml)


coeff = (polyeval(abl,abml,q(:,1))*q(:,2))/R;
Cs = sum(polyeval(abl,abml,R).*coeff);


ym2=[0;0;y0];

for i = 0:M-1
    y(i+1) = D*j*R^i - i*D*Cs*R^(i-1) + (i)*(i-1)*D*ym2(i+1);
end

