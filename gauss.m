 %模拟高斯光束在湍流传播
clc                                             
clear                                      

Cn2=2.0e-20;                                    %  湍流结构常数
lamda=0.6328e-6;                                %   波长
k=2*pi/lamda;                                   %   波数
w0=20.0e-3;                                     %   束腰半径
z=1000;                                        %   设置传输距离
deltz=200;                                      %   设置相位屏间距
L=0.6;                                          %   屏的大小
caiyang=256; 
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
r2=x.^2+y.^2;
Eu2=w0./wz.*exp(-r2./(wz.^2)).*exp(i*atan(z./f)-i*k.*(z+r2./(2.*R)));
%Eu2=w0./wz.*exp(-(r2./(wz.^2)+i*k*r2/(2*R)));         %基模高斯光束的场



[x1,y1]=meshgrid(-caiyang/2:1:caiyang/2-1,-caiyang/2:1:caiyang/2-1);       
%[x1,y1]=meshgrid(1:1:caiyang,1:1:caiyang); 
kr=sqrt((2*pi*x1/L).^2+(2*pi*y1/L).^2);      %k为三维空间的波数，kr=(kx2+ky2)的1/2次方
pusai=2*pi*k.^2*0.033*Cn2*(kr).^(-11/3)*deltz;       %大气相位功率谱函数表达式
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


    Eu1=fft2(Eu2.*exp(i.*ping));
    aa1=exp((-i.*deltz./(2*k)).*((2*pi*x1./L).^2+(2*pi*y1./L).^2));
    aa2=aa1.*Eu1;
    Eu2=abs(ifft2(aa2));

    shijian=num2str(l/bushu*100);
    shijian=num2str(fix(l/bushu*100));
    waitbar(l/bushu,h,['请等待，已完成',shijian,'%']);
end
close(h);
    
I=Eu2.^2;
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