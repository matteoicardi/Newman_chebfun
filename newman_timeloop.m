    t = tt(i);
    disp([num2str(i), ' of ',num2str(nt) , ' ',...
        num2str(t), ' of ', num2str(ttot)])
%     figure(1)
%     plot(ps1)
%     hold on
    if mod(i,10)==0 && draw
        figure(1)
        plot(ps)
        hold on
        plot(pe)
        figure(2)
        plot(ce)
        hold on
        plot(cs)
        figure(3)
        plot(j(ce,cs,pe,ps)./F./cmax)
        hold on
        res_ce(i,:)=(ce);
        res_cs1(i,:)=(cs1);
        res_cs2(i,:)=(cs2);
        res_pe(i,:)=(pe);
        res_ps1(i,:)=(ps1);
        res_ps2(i,:)=(ps2);
    end
    if draw
        figure(4)
        hold off
        plot(SoC2(1:i-1),VV(1:i-1))
        drawnow
    end
    %I=0;
    %if max(cs1)>.95 || max(cs2)>.95 || max(1-cs1)>.95 || max(1-cs2)>.95
    %    I = -I;
    %    init = 0;
    %    disp('Charge/Discharge finished, invert current')
    %end
    pde_ps1.op = @(x,ps) -sigma1*diff(ps,2) + j1(ce,cs1,pe,ps);
    pde_ps1.lbc = @(ps) sigma1*diff(ps)-I(t)/A1 ;
    if init
        pde_ps1.init = ps1;
    end
    
%     plot(ps2)
    pde_ps2.op = @(x,ps) -sigma2*diff(ps,2) + j2(ce,cs2,pe,ps);
    pde_ps2.rbc = @(ps) sigma2*diff(ps) - I(t)/A2;
    if init
        pde_ps2.init = ps2;
    end
    
%     pde_e.op= @(ce,pe) [De*diff(ce,2) + (1-t0)*j(ce,cs,pe,ps)/F/n/cmax - ce/dt;...
%         -k*diff(pe,2)-diff(kD*diff(ce)./ce) - j(ce,cs,pe,ps)];
%     ce = pde_e \ [(ce/dt); 0];
%     plot(pe)
    pde_pe.op= @(x,pe) -k.*diff(pe,2)-kD.*diff(diff(ce)./ce) - j(ce,cs,pe,ps);
    if init
        pde_pe.init = pe;
    end
    
%     figure(2)
%     plot(ce)
%     hold on
    pde_ce.op= @(x,ce) De.*diff(ce,2) + ...
        (1-t0).*j(ce,cs,pe,ps)/F/cmax./(1-nn) - ce/dt;
       % (- 15.99*ce.^6 + 53.59*ce.^5 - 70.62*ce.^4 + 46.64*ce.^3 - 16.62*ce.^2 + 3.945*ce + 0.05083)*...
    %if init
        pde_ce.init = ce;
    %end

%     pde_cs1.op= @(cs) Ds1*diff(cs,2) - ...
%         (1-t0)*j1(ce,cs1,pe,ps1)/F/cmax/n1 - cs/dt;
%        % (- 310*cs.^8 + 1.24e3*cs.^7 - 2.06e3*cs.^6 + 1.83e3*cs.^5 - 943*cs.^4 + 286*cs.^3 - 51*cs.^2 + 5.8*cs + 0.0266)*...
%     pde_cs1.init = cs1;
%     cs1 = pde_cs1 \(-cs1/dt);
% 
%     pde_cs2.op= @(cs) Ds2*diff(cs,2) - ...
%         (1-t0)*j2(ce,cs2,pe,ps2)/F/cmax/n2 - cs/dt;
%        % (- 310*cs.^8 + 1.24e3*cs.^7 - 2.06e3*cs.^6 + 1.83e3*cs.^5 - 943*cs.^4 + 286*cs.^3 - 51*cs.^2 + 5.8*cs + 0.0266)*...
%     pde_cs2.init = cs2;
%     cs2 = pde_cs2 \(-cs2/dt);
%     
%     cs = chebfun(@(x) cs1(x).*IOm1(x)+cs2(x).*IOm3(x)+.5*IOm2(x),Om);

%     pde_cs1.op = @(x,y,Cs) diff(Cs,2,1)  + diff(Cs,2,2) ... + 2*diff(Cs,1,2)./y...
%         - Cs/dt/Ds1 + Cs1/dt/Ds1;
%     pde_cs1.ubc = @(x,Cs) diff(Cs)+ R1*j1(ce,cs1,pe,ps1)/(3*n1*F)/cmax;
%     Cs1 = pde_cs1 \ 0;
%     
%     pde_cs2.op = @(x,y,Cs) diff(Cs,2,1)  + diff(Cs,2,2) ... + 2*diff(Cs,1,2)./y ...
%         - Cs/dt/Ds2 + Cs2/dt/Ds2;
%     pde_cs2.ubc = @(x,Cs) diff(Cs)+R2*j2(ce,cs2,pe,ps2)/(3*n2*F)/cmax;
%     Cs2 = pde_cs2 \ 0;
%     cs1 = chebfun(@(t) feval(Cs1,t,R1)',Om1);
%     cs2 = chebfun(@(t) feval(Cs1,t,R2)',Om2);


    jj1 = j1(ce,cs1,pe,ps1)/F/cmax/n1/(R1)*4*pi;
    jj2 = j2(ce,cs2,pe,ps2)/F/cmax/n2/(R2)*4*pi;
    m1 = m1 - eig_sph_diff(t,m1,{ncoef1,vcoef1,dcoef1,jj1,dj1,Ds1})*dt;
    m2 = m2 - eig_sph_diff(t,m2,{ncoef2,vcoef2,dcoef2,jj2,dj2,Ds2})*dt;

    
    % compute all
    ps1_n = pde_ps1 \ 0;
    ps2_n = pde_ps2 \ 0;
    ps_n = chebfun(@(x) ps1(x).*IOm1(x)+ps2(x).*IOm3(x),Om);
    pe_n = pde_pe \0;
    pe_n = pe_n - mean(pe_n);
    ce_n = pde_ce \(-ce/dt);
    cs1_n=feval(ef1,R1)'*m1-jj1.*v1(R1);
    cs2_n=feval(ef2,R2)'*m2-jj2.*v2(R2);
    cs1_n=chebfun(cs1_n{1});
    cs2_n=chebfun(cs2_n{1});
    cs_n = chebfun(min(max(chebfun(@(x) cs1_n(x).*IOm1(x)+cs2_n(x).*IOm3(x),Om),1e-3),1-1e-3));
    
    % update all
    ce=ce_n;
    cs1=cs1_n;
    cs2=cs2_n;
    cs=cs_n;
    ps1=ps1_n;
    ps2=ps2_n;
    ps=ps_n;
    pe=pe_n;

    dj1 = (j1(ce,cs1,pe,ps1)/F/cmax/n1/(R1)*4*pi-jj1)/dt;
    dj2 = (j2(ce,cs2,pe,ps2)/F/cmax/n2/(R2)*4*pi-jj2)/dt;
    

    % QoI
    VV(i) = ps1(0)-ps2(L);
    II(i) = I(t);
    SoC(i) = mean(cs1);
    SoC2(i) = mean(cs2);
    
    disp([VV(i), II(i), SoC(i), SoC2(i)])

    if I(t)*I(t+dt)<0
        init = 0;
    else
        init = 1;
    end
        
    i=i+1;