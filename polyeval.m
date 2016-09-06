function p=polyeval(ab,abm,x)
x=reshape(x,[1,length(x)]);
M = size(ab,1);
%f = ones(M,1);
f(1,:)=zeros(length(x),1);
f(2,:)=ones(length(x),1);


for i=1:M-1
    f(i+2,:) = (x-ab(i,1)).*f(i+1,:)-ab(i,2)*f(i,:);
end

p=sqrt(abm(1))*f(2:end,:)./(sqrt([abm])*ones(1,length(x)));