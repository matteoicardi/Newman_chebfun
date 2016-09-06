function m=multimean(a,d)
% compute mean over multiple dimensions
m=a;
for i=1:length(d)
    m=mean(m,d(i));
end
