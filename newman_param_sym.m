%% coefficients (from Uddin, Perera, Widanalage)
RT = 8.3144598 * 300; % gas constant times temperature
R1 = 5e-6; % radius active particles
R2 = 5e-6; % radius active particles
F  = 96487;   % faraday constant
A  = 0.04;   % electrode surface/cell width*depth
A1  = 0.04;   % electrode surface/cell width*depth
A2  = 0.04;   % electrode surface/cell width*depth
t0 = .363;  % constant
n1  = .3;  % porosity
n2  = .3;  % porosity
ne  = .5;
brugg = 3./2;
De = 2e-10;  % bulk, overwritten
Dee = De * ne^(brugg);
De1 = De * n1^(brugg);
De2 = De * n2^(brugg);
Ds1 = 1e-14; % internal diffusion
Ds2 = 1e-14; % internal diffusion
Ds1x = Ds1 * n1^(brugg);  % inter-particle diffusion
Ds2x = Ds2 * n2^(brugg);
sigma1 = 100; % electrical conductivity
sigma2 = 100; % electrical conductivity
k  = 1;   % bulk overwritten, Fraunhofer thesis, Dao et al 2012
k1  = k * n1^(brugg);   
ke  = k * ne^(brugg);   
k2  = k * n2^(brugg);
kD = k * 2* RT*(1-t0)/F;
kD1 = k1 * 2* RT*(1-t0)/F;
kD2 = k2 * 2* RT*(1-t0)/F;
kDe = ke * 2* RT*(1-t0)/F;

alpha  = .5; % constant
K1  = 1e-11;   % reaction kinetic constant
K2  = 1e-11;   % reaction kinetic constant
L  = 3e-4; % cell length
lambda = .45; % relative length of catode/anode
L1 = lambda*L;
L2 = (1-lambda)*L;
cmax = 50000;  % concentration max
hour = 3600;
I_ref  = A*L1*n1*F*cmax/hour;   % applied current  % COMPUTE based on 1C

% % overpotential (Doyle Foller Newman 1993)
% cs = chebfun(@(cs) cs, [0,1]);
% U01  = (-4.656 + 88.669*cs.^2 - 401.119*cs.^4 + 342.909*cs.^6 - 462.471*cs.^8 + 433.434*cs.^10) ...
%        ./(-1 +18.933*cs.^2 - 79.532*cs.^4 + 37.311*cs.^6 - 73.083*cs.^8 + 95.96*cs.^10);
% U02  = @(cs) cs;

% overpotential (Ramos)
% U01  = chebfun(@(cs) RT*log(K1*(1-(limit01(cs)))./(limit01(cs)))/F,[-1,2]);
% U02  = chebfun(@(cs) RT*log(K2*(1-(limit01(cs)))./(limit01(cs)))/F,[-1,2]);
 U01  = @(y) (RT*log(K1*(y)./(1-y))/F)+5;  % cathode
 U02  = @(y) (RT*log(K2*(y)./(1-y))/F);  % anode
