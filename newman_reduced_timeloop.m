    t = tt(i);
    disp([num2str(i), ' of ',num2str(nt) , ' ',...
        num2str(t), ' of ', num2str(ttot)])    

    %if max(cs1)>.95 || max(cs2)>.95 || max(1-cs1)>.95 || max(1-cs2)>.95
    %    I = -I;
    %    init = 0;
    %    disp('Charge/Discharge finished, invert current')
    %end

    %% update operators (implicitly in time)
    newman_update_operators


    %% compute solid potentials
    ps1_n = pde_ps1 \ 0;
    ps2_n = pde_ps2 \ 0;
    ps_n = chebfun(@(x) ps1_n(x).*IOm1(x)+ps2_n(x).*IOm3(x),Om);
    if ~update_after
        ps1=ps1_n;
        ps2=ps2_n;
        ps=ps_n;
        newman_update_operators
    end
    
    
    %% update solid diffusion coefficients explicitly
    m10 = m1;
    m20 = m2;
    dj1 = 0;
    dj2 = 0;
    jj10 = -j1(ce,cs1,pe,ps1)/3/F/cmax/n1*(R1);
    jj20 = -j2(ce,cs2,pe,ps2)/3/F/cmax/n2*(R2);
    jj1 = jj10;
    jj2 = jj20;
    ps1_n = ps1;
    ps2_n = ps2;
    for iter=1:3
        m1 = m10 + eig_sph_diff(t,m1,{ncoef1,vcoef1,dcoef1,jj1,dj1,Ds1})*dt;
        m2 = m20 + eig_sph_diff(t,m2,{ncoef2,vcoef2,dcoef2,jj2,dj2,Ds2})*dt;
        for jter=1:3
            cs1_n=chebfun(feval(ef1,R1)'*m1+jj1.*v1(R1)/Ds1);
            cs2_n=chebfun(feval(ef2,R2)'*m2+jj2.*v2(R2)/Ds2);
            jj1 = -j1(ce,cs1_n,pe,ps1_n)/3/F/cmax/n1*(R1);
            jj2 = -j2(ce,cs2_n,pe,ps2_n)/3/F/cmax/n2*(R2);
            %disp(norm(jj1))
        end
        dj1 = (jj1-jj10)/dt;
        dj2 = (jj2-jj20)/dt;
        pde_ps1.op = @(x,ps) sigma1*diff(ps,2) - j1(ce,cs1_n,pe,ps);
        pde_ps2.op = @(x,ps) sigma2*diff(ps,2) - j2(ce,cs2_n,pe,ps);
        % compute solid potentials
        pde_ps1.init=ps1_n;
        pde_ps2.init=ps2_n;
        ps1_n = pde_ps1 \ 0;
        ps2_n = pde_ps2 \ 0;
        ps_n = chebfun(@(x) ps1_n(x).*IOm1(x)+ps2_n(x).*IOm3(x),Om);

        disp([sum(jj1) sum(dj1) mean(cs1_n) mean(ps1_n)])
        
    end
    
    %% evaluate solid concentration at surface
    cs1_n=feval(ef1,R1)'*m1+jj1.*v1(R1)/Ds1;
    cs2_n=feval(ef2,R2)'*m2+jj2.*v2(R2)/Ds2;
    cs1_n=chebfun(cs1_n{1});
    cs2_n=chebfun(cs2_n{1});
    cs_n = chebfun(min(max(chebfun(@(x) cs1_n(x).*IOm1(x)+cs2_n(x).*IOm3(x),Om),1e-3),1-1e-3));
    if ~update_after
        cs1=cs1_n;
        cs2=cs2_n;
        cs=cs_n;
        newman_update_operators
    end
    
    
    %% update all fields
    if update_after
        cs1=cs1_n;
        cs2=cs2_n;
        cs=cs_n;
        ps1=ps1_n;
        ps2=ps2_n;
        ps=ps_n;
    end
        
    % store and plot
    if store
        res_ce(i+1,:)=(ce);
        res_cs1(i+1,:)=(cs1);
        res_cs2(i+1,:)=(cs2);
        res_pe(i+1,:)=(pe);
        res_ps1(i+1,:)=(ps1);
        res_ps2(i+1,:)=(ps2);
    end
    if mod(i,10)==0 && draw
        figure(1)
        plot(ps)
        hold on
        figure(2)
        hold on
        plot(cs)
        figure(3)
        plot(j(ce,cs,pe,ps)./F./cmax)
        hold on
        figure(4)
        hold off
        plot(SoC2(1:i-1),VV(1:i-1))
        figure(5)
        plot(feval(m1,0).*ef1)
        drawnow
    end

    % QoI
    VV(i) = ps1(0)-ps2(L);
    II(i) = I_imposed;
    SoC(i) = mean(ef1{end})*mean(chebfun(m1(end,:)));
    SoC2(i) = mean(ef2{end})*mean(chebfun(m2(end,:)));
    
    disp([VV(i), II(i), SoC(i), SoC2(i)])

    if I(t)*I(t+dt)<0
        init = 0;
    else
        init = 1;
    end
        
    i=i+1;