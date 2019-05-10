 %ģ���˹��������������
clc                                             
clear                                      

Cn2=2.0e-20;                                    %  �����ṹ����
lamda=0.6328e-6;                                %   ����
k=2*pi/lamda;                                   %   ����
w0=20.0e-3;                                     %   �����뾶
z=1000;                                        %   ���ô������
deltz=200;                                      %   ������λ�����
L=0.6;                                          %   ���Ĵ�С
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
%Eu2=w0./wz.*exp(-(r2./(wz.^2)+i*k*r2/(2*R)));         %��ģ��˹�����ĳ�



[x1,y1]=meshgrid(-caiyang/2:1:caiyang/2-1,-caiyang/2:1:caiyang/2-1);       
%[x1,y1]=meshgrid(1:1:caiyang,1:1:caiyang); 
kr=sqrt((2*pi*x1/L).^2+(2*pi*y1/L).^2);      %kΪ��ά�ռ�Ĳ�����kr=(kx2+ky2)��1/2�η�
pusai=2*pi*k.^2*0.033*Cn2*(kr).^(-11/3)*deltz;       %������λ�����׺������ʽ
[m n]=find(pusai==inf);
pusai(m,n)=pusai(m-1,n);
pusai=fftshift(pusai);     %����Ҷ���任

bushu=z/deltz;
h=waitbar(0,'�����У���ȴ�...');
for l=1:bushu;
   
ra=randn(caiyang,caiyang);         %���ֵ����λ����ĸ�˹�����   
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
    waitbar(l/bushu,h,['��ȴ��������',shijian,'%']);
end
close(h);
    
I=Eu2.^2;
Iu=I/max(max(I));
%Iu=I;

figure(3);
mesh(x,y,Iu);
%figure(2);
%mesh(ping);
xlabel('ƽ��X����/m');
ylabel('ƽ��Y����/m');
zlabel('��ǿ');
figure(4)
imshow(mat2gray(Iu))
xlabel('ƽ��X����/m');
ylabel('ƽ��Y����/m');