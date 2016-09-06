%% coefficients (from Uddin, Perera, Widanalage)
RT = 8.3144598 * 300; % gas constant times temperature
R1 = 10.7e-6; % radius active particles
R2 = 5.7e-6; % radius active particles
F  = 96487;   % faraday constant
A  = 0.04;   % electrode surface/cell width*depth
A1  = 0.04;   % electrode surface/cell width*depth
A2  = 0.04;   % electrode surface/cell width*depth
t0 = .363;  % constant
n1  = .3;  % porosity
n2  = .4;  % porosity
ne  = .5;
brugg = 3./2;
De = 2.6e-10;  % bulk, overwritten
Dee = De * ne^(brugg);
De1 = De * n1^(brugg);
De2 = De * n2^(brugg);
Ds1 = 1e-14; % internal diffusion
Ds2 = 3.9e-14; % internal diffusion
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
K1  = 2.344e-11;   % reaction kinetic constant
K2  = 5.0307e-11;   % reaction kinetic constant
L  = 3e-4; % cell length
lambda = .45; % relative length of catode/anode
L1 = lambda*L;
L2 = (1-lambda)*L;
cmax = 50000;  % concentration max
hour = 3600;
I_ref  = A*L1*n1*F*cmax/hour;   % applied current  % COMPUTE based on 1C