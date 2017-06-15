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
    ps1_n = ps1; %pde_ps1 \ 0;
%    ps2_n = ps2; %pde_ps2 \ 0;
%    ps_n = ps; %chebfun(@(x) ps1_n(x).*IOm1(x)+ps2_n(x).*IOm3(x),Om);
%     if ~update_after
%         ps1=ps1_n;
%         ps2=ps2_n;
%         ps=ps_n;
%         newman_update_operators
%     end
    
relax = 1;
    
    %% update solid diffusion coefficients explicitly
    m10 = m1;
%    m20 = m2;
    dj1 = 0;
%    dj2 = 0;
    jj10 = jj1;
%    jj20 = jj2;
%     ps1_n = ps1;
%     ps2_n = ps2;
    cs1_n = cs1;
%    cs2_n = cs2;
    for iter=1:100
        %% %% compute solid potentials
        % non linear versions
        pde_ps1.op = @(x,ps) sigma1*diff(ps,2) - j1(ce,cs1_n,pe,ps);
%        pde_ps2.op = @(x,ps) sigma2*diff(ps,2) - j2(ce,cs2_n,pe,ps);
        %% linearised versions
%         pde_ps1.op = @(x,ps) sigma1*diff(ps,2) - j1(ce,cs1_n,pe,ps1) ...
%             + j1dps(ce,cs1_n,pe,ps1).*(ps1_n-ps);
%         pde_ps2.op = @(x,ps) sigma2*diff(ps,2) - j2(ce,cs2_n,pe,ps2)...
%             + j2dps(ce,cs2_n,pe,ps2).*(ps2_n-ps);
        if init
            pde_ps1.init=ps1_n;
%            pde_ps2.init=ps2_n;
        end
        ps1_n = (relax*(pde_ps1 \ 0) + (1-relax)*ps1_n);
%        ps2_n = (relax*(pde_ps2 \ 0) + (1-relax)*ps2_n);
      
%         % fake solver for potential
%         ps1_n = chebfun(@(x) U01(mean(cs1)), Om1)+I_imposed/A1/sigma1/2*(-x1.^2 + 2*L1*x1);
%         ps2_n = chebfun(@(x) U02(mean(cs2)), Om2)+I_imposed/A2/sigma2/2*(-(L-x2).^2 + 2*L1*(L-x2));
        
        %% %% compute solid concentration
        jj1 = -j1(ce,cs1_n,pe,ps1_n)/3/F/cmax/n1*(R1);
%        jj2 = -j2(ce,cs2_n,pe,ps2_n)/3/F/cmax/n2*(R2);
        jjcheck = jj1;
        pscheck = ps1_n;
        m1 = m10 + eig_sph_diff(t,m1,{ncoef1,vcoef1,dcoef1,jj1,dj1,Ds1})*dt;
%        m2 = m20 + eig_sph_diff(t,m2,{ncoef2,vcoef2,dcoef2,jj2,dj2,Ds2})*dt;
        %for jter=1:3
            cs1_n=chebfun(feval(ef1,R1)'*m1+jj1.*v1(R1)/Ds1);
%           cs2_n=chebfun(feval(ef2,R2)'*m2+jj2.*v2(R2)/Ds2);
            jj1 = -j1(ce,cs1_n,pe,ps1_n)/3/F/cmax/n1*(R1);
%            jj2 = -j2(ce,cs2_n,pe,ps2_n)/3/F/cmax/n2*(R2);
            %disp(norm(jj1))
        %end
        dj1 = (jj1-jj10)/dt;
%        dj2 = (jj2-jj20)/dt;

        psnormcheck = norm(ps1_n-pscheck)/norm(pscheck);
        normcheck = norm(jj1-jjcheck)/norm(jjcheck);
        %disp([sum(jj1) sum(dj1) mean(cs1_n) mean(ps1_n) normcheck psnormcheck])
        if ((iter>2 && max(normcheck,psnormcheck)<releps*relax) || isinf(normcheck))
            break
        end
        
    end
    
    jj(i)=mean(jj1);
    %% evaluate solid concentration at surface
    cs1_n=feval(ef1,R1)'*m1+jj1.*v1(R1)/Ds1;
%    cs2_n=feval(ef2,R2)'*m2+jj2.*v2(R2)/Ds2;
    cs1_n=chebfun(cs1_n{1});
%    cs2_n=chebfun(cs2_n{1});
    mcs1 = mean(ef1{end})*(chebfun(m1(end,:)));
%    mcs2 = mean(ef2{end})*(chebfun(m2(end,:)));
    
    % merge domains
%    ps_n = chebfun(@(x) ps1_n(x).*IOm1(x)+ps2_n(x).*IOm3(x),Om,'splitting','on');
%    cs_n = chebfun(@(x) cs1_n(x).*IOm1(x)+.5.*IOm2(x)+cs2_n(x).*IOm3(x),Om,'splitting','on');
%    Cs1  = chebfun2(@(x,y) (feval(ef1,y))'*feval(m1,x)+jj1(x)*(v1(y))/Ds1, [Om1,Oms1],'vectorize');
%    Cs2  = chebfun2(@(x,y) (feval(ef2,y))'*feval(m2,x)+jj2(x)*(v2(y))/Ds2, [Om2,Oms2],'vectorize');
    P1(i)   = mean(cs1_n-mcs1);
%    P2(i)   = mean(cs2_n-mcs2);
    
    cs1=cs1_n;
%    cs2=cs2_n;
%    cs=cs_n;
    ps1=ps1_n;
%    ps2=ps2_n;
%    ps=ps_n;
    newman_update_operators
    
    
    % QoI
    VV(i) = ps1(0)-ps1(L1);
    II(i) = I_imposed;
    SoC(i) = mean(mcs1);
%    SoC2(i) = mean(mcs2);
    
    disp([VV(i), II(i), SoC(i), P1(i), iter])

    % store and plot
    if store
%        res_ce(i+1,:)=(ce);
        res_cs1(i+1,:)=(cs1);
%        res_cs2(i+1,:)=(cs2);
%        res_pe(i+1,:)=(pe);
        res_ps1(i+1,:)=(ps1);
%        res_ps2(i+1,:)=(ps2);
    end
    if mod(i,drawstep)==0 && draw
        figure(1)
        plot(ps1)
        hold on
        figure(2)
        hold on
        plot(cs1)
        figure(3)
        plot(j1(ce,cs1,pe,ps1)./F./cmax)
        hold on
        figure(4)
        hold off
        plot(SoC(1:i-1),VV(1:i-1))
        Cs1=chebfun2(@(x,y) feval(ef1,y)'*feval(m1,x)+jj1(x)*v1(y)/Ds1,[Om1,Oms1],'vectorize');
        figure(5)
        hold on
        plot(Cs1)
        drawnow
    end


%     if I(t)*I(t+dt)<0
%         init = 0;
%     else
%         init = 1;
%     end
        
    i=i+1;