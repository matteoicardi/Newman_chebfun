newman_param;

syms xce xpe xcs xps
xj =  3*ce_sq(xce).*cs_sq(xcs).*(n1*K1*sinh(alpha*F*(xps-U01(xcs)-xpe)/RT));
h=hessian(xj,[xce,xcs,xpe,xps]);
g=gradient(xj,[xce,xcs,xpe,xps]);

cc = chebfun(@(c) c, [.1,.9]);
Ice = linspace(min(cc),max(cc),5);
Ics = Ice;
Ips = linspace(min(U01(cc)),max(U01(cc)),5);
Ipe = Ips;
Ivol = [ Ice(end)-Ice(1) Ics(end)-Ics(1) Ipe(end)-Ipe(1) Ips(end)-Ips(1)];
Ivol2 = Ivol'*Ivol;
[X0,Y0,Xce,Xcs,Xpe,Xps]=ndgrid([0],[0],Ice,Ics,Ipe,Ips);
%a=single(subs(subs(h,xpe,0),{xce,xcs,xps},{Xce,Xcs,Xps}));
gg=single(subs(g,{xce,xcs,xpe,xps},{Xce,Xcs,Xpe,Xps}));
%hh=single(subs(h,{xce,xcs,xpe,xps},{Xce,Xcs,Xpe,Xps}));
multimean(abs(gg),3:6).*Ivol'
