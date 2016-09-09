pe=res_pe{i};
ps1=res_ps1{i};
ps2=res_ps2{i};
ce=res_ce{i};
cs1=res_cs1{i};
cs2=res_cs2{i};
ps = chebfun(@(x) ps1(x).*IOm1(x)+ps2(x).*IOm3(x),Om);
cs = chebfun(min(max(chebfun(@(x) cs1(x).*IOm1(x)+cs2(x).*IOm3(x),Om),1e-3),1-1e-3));
res_pe=res_pe{1:i};
res_ce=res_ce{1:i};
res_ps1=res_ps1{1:i};
res_ps2=res_ps2{1:i};
res_cs1=res_cs1{1:i};
res_cs2=res_cs2{1:i};

dj1=0;
dj2=0;