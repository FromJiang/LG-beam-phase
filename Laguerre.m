%缔合拉盖尔高斯多项式
function y=Laguerre(p,m,x)
y=zeros(size(x));
for k=0:p
y=y+(-1)^k*factorial(m+p)*x.^k/factorial(p-k)/factorial(m+k)/factorial(k);
end
end

