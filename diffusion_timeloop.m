    t = tt(i);
    disp([num2str(i), ' of ',num2str(nt) , ' ',...
        num2str(t), ' of ', num2str(ttot)])    

    %if max(cs1)>.95 || max(cs2)>.95 || max(1-cs1)>.95 || max(1-cs2)>.95
    %    I = -I;
    %    init = 0;
    %    disp('Charge/Discharge finished, invert current')
    %end

    %% update solid diffusion coefficients explicitly
    jc=1e-5*t;
    m10 = m1;
    m20 = m2;
    jj10 = cs1.*(1-cs1).*x1*jc;
    jj20 = -cs2.*(1-cs2).*(L-x2)*jc;
    if i==0
        dj1 = jj10; % when starting from constant solution
        dj2 = jj20; % the derivative is non zero
    end
    jj1 = jj10;
    jj2 = jj20;
    %m1corr=-chebmatrix(vcoef1')*chebmatrix(jj1)/Ds1;
    %m2corr=-chebmatrix(vcoef2')*chebmatrix(jj2)/Ds2;        
    for iter=1:5
        m1 = m10 + eig_sph_diff(t,m1,{ncoef1,vcoef1,dcoef1,jj1,dj1,Ds1})*dt;
        m2 = m20 + eig_sph_diff(t,m2,{ncoef2,vcoef2,dcoef2,jj2,dj2,Ds2})*dt;
        for jter=1:3
            cs1_n=chebfun(feval(ef1,R1)'*m1+jj1.*v1(R1)/Ds1);
            cs2_n=chebfun(feval(ef2,R2)'*m2+jj2.*v2(R2)/Ds2);
            jj1 = cs1_n.*(1-cs1_n).*x1*jc;
            jj2 = -cs2_n.*(1-cs2_n).*(L-x2)*jc;
            %disp(norm(jj1))
        end
        dj1 = (jj1-jj10)/dt;
        dj2 = (jj2-jj20)/dt;

        disp([sum(jj1) sum(jj2) sum(dj1) mean(cs1_n) mean(cs2_n) mean(ef1{end})*mean(chebfun(m1(end,:))) mean(ef2{end})*mean(chebfun(m2(end,:)))])
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
    end
    
    
    %% update all fields
    if update_after
        cs1=cs1_n;
        cs2=cs2_n;
        cs=cs_n;
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
        figure(2)
        hold on
        plot(cs1)
        figure(3)
        plot(cs2)
        hold on
        figure(5)
        hold off
        plot(feval(m1,0).*ef1)
        hold on
        plot(feval(m1,0.00014).*ef1)
        drawnow
    end


    if I(t)*I(t+dt)<0
        init = 0;
    else
        init = 1;
    end
        
    i=i+1;