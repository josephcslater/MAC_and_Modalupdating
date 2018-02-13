function a=hist2cdf(b)
  
lb=length(b);
a=zeros(size(b));
for i=1:lb
  a(i)=sum(b(1:i));
end
a=a/max(a);