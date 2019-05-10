%%%LG光束的光强
clear
N = 512;
delta =0.003;
% c=-1:0.01:1;
% r=-1:0.01:1;
x=(-N/2:N/2-1)*delta;
y=(-N/2:N/2-1)*delta;
[x,y]=meshgrid(x,y);
% [x,y]=meshgrid(c,r);
[theta,rho]=cart2pol(x,y);
w0=0.3;                      %束腰半径
m=3;p=3;z=0;                 %确定光束的阶数m,模数p，传输距离z
lambda=0.632e-3;k=2*pi/lambda;
zr=pi*w0^2/lambda;
wz=w0*sqrt(1+(z/zr)^2);
A=2*factorial(p)/(pi*factorial(m+p));    %归一化因子
u=(sqrt(2)*rho/w0).^m.*exp(-rho.^2/w0^2).*exp(1i*m*theta)*sqrt(2/factorial(m)/pi).*Laguerre(p,m,2*rho.^2/w0^2);
%aaa= A/wz^2.*(2*rho.^2/wz^2).^m.*(Laguerre(p,m,2*rho.^2/wz^2)).^2.*exp(-2*rho.^2/wz^2).*exp(1i*m*theta);    %光场表达式
beta=10;
I=u.*conj(u);
Iu=I/max(max(I));
%画光强相位
figure
surf(x,y,Iu);
shading interp                  %色彩插值处理
zlabel('intensity');

figure
pha=rem(angle(u)+pi*2,2*pi);
surf(x,y,pha);
shading interp                
%colormap gray;

%画光强相位像素图
figure
imshow(Iu)
figure
vortex_phase=rem(angle(u)+2*pi,2*pi);
imshow(vortex_phase,[])
