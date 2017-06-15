    t = tt(i);
    disp([num2str(i), ' of ',num2str(nt) , ' ',...
        num2str(t), ' of ', num2str(ttot)])    

    %if max(cs1)>.95 || max(cs2)>.95 || max(1-cs1)>.95 || max(1-cs2)>.95
    %    I = -I;
    %    init = 0;
    %    disp('Charge/Discharge finished, invert current')
    %end

    %% update solid diffusion coefficients explicitly
    jc=I(t)/F/A1;
    m10 = m1;
    jj10 = cs1.*(1-cs1).*x1*jc;
    if i==0
        dj1 = jj10; % when starting from constant solution
    end
    jj1 = jj10;
    %m1corr=-chebmatrix(vcoef1')*chebmatrix(jj1)/Ds1;
    for iter=1:5
        m1 = m10 + eig_sph_diff(t,m1,{ncoef1,vcoef1,dcoef1,jj1,dj1,Ds1})*dt;
        for jter=1:3
            cs1_n=chebfun(feval(ef1,R1)'*m1+jj1.*v1(R1)/Ds1);
            jj1 = cs1_n.*(1-cs1_n).*x1*jc;
            %disp(norm(jj1))
        end
        dj1 = (jj1-jj10)/dt;

        disp([sum(jj1) sum(dj1) mean(cs1_n)  mean(ef1{end})*mean(chebfun(m1(end,:))) ])
    end
    
    %% evaluate solid concentration at surface
    cs1_n=feval(ef1,R1)'*m1+jj1.*v1(R1)/Ds1;
    cs1_n=chebfun(cs1_n{1});
    if ~update_after
        cs1=cs1_n;
    end
    
    
    %% update all fields
    if update_after
        cs1=cs1_n;
    end
        
    % store and plot
    if store
        res_cs1(i+1,:)=(cs1);
    end
    if mod(i,10)==0 && draw
        figure(2)
        hold on
        plot(cs1)
        figure(5)
        hold off
        plot(feval(m1,0).*ef1)
        hold on
        plot(feval(m1,0.00014).*ef1)
        drawnow
    end

        
    i=i+1;