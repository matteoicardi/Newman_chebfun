    %% smoothed ramp to initialise simulation
    I_imposed = I(t);%/(1+exp((-t+dt)));

    %% solid potentials
    % solve for variable (phi_s-phi_e)
    pde_ps1.op = @(x,ps) sigma1*diff(ps,2) - j1(ce,cs1,pe,ps);
    pde_ps1.rbc = @(ps) diff(ps) + feval(diff(pe),L1) ;
    pde_ps1.lbc = @(ps) sigma1*diff(ps)-I_imposed/A1 ;
    if init
        pde_ps1.init = ps1-restrict(pe,Om1);
    end
    
%     plot(ps2)
    pde_ps2.op = @(x,ps) sigma2*diff(ps,2) - j2(ce,cs2,pe,ps);
    pde_ps2.lbc = @(ps) diff(ps) + feval(diff(pe),L2) ;
    pde_ps2.rbc = @(ps) sigma2*diff(ps) - I_imposed/A2;
    if init
        pde_ps2.init = ps2-restrict(pe,Om2);
    end
    
    %% electrolyte potential
    
%     pde_e.op= @(ce,pe) [De*diff(ce,2) + (1-t0)*j(ce,cs,pe,ps)/F/n/cmax - ce/dt;...
%         -k*diff(pe,2)-diff(kD*diff(ce)./ce) - j(ce,cs,pe,ps)];
%     ce = pde_e \ [(ce/dt); 0];
%     plot(pe)
    pde_pe.op= @(x,pe) k.*diff(pe,2) + kD.*diff(diff(ce)./ce);
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

