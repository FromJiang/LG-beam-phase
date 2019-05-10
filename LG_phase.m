 %模拟LG光束在湍流传播
clc                                             
clear                                      

Cn2=2.0e-16;                                    %   湍流结构常数
lamda=0.6328e-6;                                %   波长
k=2*pi/lamda;                                   %   波数
w0=20.0e-2;                                     %   束腰半径
z=2000;                                         %   设置传输距离
deltz=200;                                      %   设置相位屏间距
L=1.7;                                          %   屏的大小
caiyang=512; 
delta=L/caiyang;
N=512; 
dfx=1/(N*delta);
dfy=1/(N*delta);
buchang=L/caiyang;
C=2*pi/L;
%power=1;
%A0=(power/(pi*w0.^2))^(1/2); 
f=pi*w0^2/lamda;
k=2*pi/lamda;
wz=w0*sqrt(1+(z/f)^2);
R=z+f^2/z;

[x,y]=meshgrid(-L/2:buchang:L/2-buchang,-L/2:buchang:L/2-buchang);
%[x,y]=meshgrid(0:buchang:L-buchang,0:buchang:L-buchang);
%Eu2=A0*(exp(-(x.^2+y.^2)/w0.^2)).^(1/2);
%r2=x.^2+y.^2;
[theta,r]=cart2pol(x,y);

% 贝塞尔涡旋光束
% m=3;p=8;    % 拓扑荷数、径向指数
% Ou=p/w0;
%  F=inf;
%  w=1/(k*w0^2)+i/F;
%  u0=besselj(m,Ou*r).*exp(-1*k*r.^2*w).*exp(-i*m*theta);

%  拉盖尔涡旋光束
m=3;p=0;                     %确定光束的阶数m，模数p
u0=(sqrt(2)*r/w0).^m.*exp(-r.^2/w0^2).*exp(1i*m*theta)*sqrt(2/factorial(m)/pi).*(Laguerre(p,m,2*r.^2/w0^2));         %光场表达式
  I=u0.*conj(u0);
  Iu=I/max(max(I));
%Iu=I;
% 
% figure;
% mesh(x,y,Iu);
%Eu2=w0./wz.*exp(-r2./(wz.^2)).*exp(i*atan(z./f)-i*k.*(z+r2./(2.*R)));
%Eu2=w0./wz.*exp(-(r2./(wz.^2)+i*k*r2/(2*R)));         %基模高斯光束的场
[x1,y1]=meshgrid(-caiyang/2:1:caiyang/2-1,-caiyang/2:1:caiyang/2-1);       
%[x1,y1]=meshgrid(1:1:caiyang,1:1:caiyang); 
l0=0.01;
L0=1;
km = 5.92/l0;
k0 = 2*pi/L0;
kr=sqrt((2*pi*x1/L).^2+(2*pi*y1/L).^2);      %k为三维空间的波数，kr=(kx2+ky2)的1/2次方
%pusai=2*pi*k.^2*0.033*Cn2*(kr).^(-11/3)*deltz;       %大气相位功率谱函数表达式
pusai=2*pi*k.^2*0.033*Cn2*deltz * exp(-(kr/km).^2) ./ (kr.^2 + k0^2).^(11/6);
[m n]=find(pusai==inf);
pusai(m,n)=pusai(m-1,n);
pusai=fftshift(pusai);     %傅里叶反变换

bushu=z/deltz;
h=waitbar(0,'计算中，请等待...');
for l=1:bushu;
   
ra=randn(caiyang,caiyang);         %零均值，单位方差的高斯随机数   
rb=randn(caiyang,caiyang);
rr=ra+i.*rb;

ping=sqrt(C)*caiyang^2*ifft2(rr.*sqrt(pusai));           
ping=real(ping); 
figure(1);
mesh(ping);
figure(2);
imshow(mat2gray(ping));


    Eu1=fft2(u0.*exp(i.*ping));
    %% 自由空间衍射
fx=(-N/2:N/2-1)*dfx;
fy=(-N/2:N/2-1)*dfy; 
[Fx,Fy]=meshgrid(fx,fy);
H=exp(1i*k*sqrt(ones(N,N)-(lamda*Fx).^2-(lamda*Fy).^2));  % 传递函数(没有乘距离)
pfac = fftshift(H);  
    u1=Eu1.*pfac.^z;         %光场表达式
    
    Eu2=abs(ifft2(u1));
    %aa1=exp((-i.*deltz./(2*k)).*((2*pi*x1./L).^2+(2*pi*y1./L).^2));
%     aa2=aa1.*Eu1;
%     Eu2=abs(ifft2(aa2));

    shijian=num2str(l/bushu*100);
    shijian=num2str(fix(l/bushu*100));
    waitbar(l/bushu,h,['请等待，已完成',shijian,'%']);
end
close(h);
    
%I=Eu2.^2;
I=Eu2.*conj(Eu2);
Iu=I/max(max(I));
%Iu=I;

figure(3);
mesh(x,y,Iu);
%figure(2);
%mesh(ping);
xlabel('平面X坐标/m');
ylabel('平面Y坐标/m');
zlabel('光强');
figure(4)
imshow(mat2gray(Iu))
xlabel('平面X坐标/m');
ylabel('平面Y坐标/m');
