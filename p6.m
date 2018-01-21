I=load('TrainingSamplesDCT_8_new.mat');
I1=imread('cheetah.bmp');
I2=imread('cheetah_mask.bmp');
[m n]=size(I1);
[m1 n1]=size(I.TrainsampleDCT_FG);
[m2 n2]=size(I.TrainsampleDCT_BG);
A=I.TrainsampleDCT_FG;
B=I.TrainsampleDCT_BG;
p1=m1/(m1+m2);
p2=m2/(m1+m2);
uC=zeros(64,1);
vC=zeros(64,64);
uG=zeros(64,1);
vG=zeros(64,64);
max=zeros(64,1);
for i=1:64
    sumF=0;
    sumB=0;
    vaF=0;
    vaB=0;
    for j=1:m1
        sumF=sumF+A(j,i);
    end
    u1=sumF/m1;
    uC(i)=u1;
    for j=1:m2
        sumB=sumB+B(j,i);
    end
    u2=sumB/m2;
    uG(i)=u2;
    for j=1:m1
        vaF=vaF+(A(j,i)-u1)^2;
    end
    v1=vaF/m1;
    vC(i,i)=v1;
    for j=1:m2
        vaB=vaB+(B(j,i)-u2)^2;
    end
    v2=vaB/m2;
    vG(i,i)=v2;
    u=u1-u2;
    v=v1+v2;
    max(i)=u^2/v;
    x=-1:0.01:1;
    y=(1/sqrt(2*pi*v1))*exp(-(x-u1).^2/2/v1);
    subplot(8,8,i);
    plot(x,y);
    hold on;
       if(i==1)
           x=-5:0.01:5
       else
           x=-1:0.01:1
       end
    y=(1/sqrt(2*pi*v2))*exp(-(x-u2).^2/2/v2);
    subplot(8,8,i);
    plot(x,y);
        str=['x=',num2str(i)];
    title(str);
end
max=abs(max);
[t index]=sort(max);
I8=[];
Iw=[];
uc8=zeros(8,1);
ug8=zeros(8,1);
ucw=zeros(8,1);
ugw=zeros(8,1);
for i=1:8
    I8(i)=index(65-i);
    Iw(i)=index(i);
    uc8(i)=uC(I8(i));
    ug8(i)=uG(I8(i));
    vc8(i,i)=vC(I8(i),I8(i));
    vg8(i,i)=vG(I8(i),I8(i));
    ucw(i)=uC(Iw(i));
    ugw(i)=uG(Iw(i));
    vcw(i,i)=vC(Iw(i),Iw(i));
    vgw(i,i)=vG(Iw(i),Iw(i));
end
for i=1:8
    u1=uc8(i);
    u2=ug8(i);
    v1=vc8(i,i);
    v2=vg8(i,i);
    x=-1:0.01:5;
    y=(1/sqrt(2*pi*v1))*exp(-(x-u1).^2/2/v1);
    subplot(4,2,i);
    plot(x,y);
    hold on;
    y=(1/sqrt(2*pi*v2))*exp(-(x-u2).^2/2/v2);
    plot(x,y);
    str=['x=',num2str(I8(i))];
    title(str);
end
figure;
for i=1:8
    u1=ucw(i);
    u2=ugw(i);
    v1=vcw(i,i);
    v2=vgw(i,i);
    x=-1:0.01:5;
    y=(1/sqrt(2*pi*v1))*exp(-(x-u1).^2/2/v1);
    subplot(4,2,i);
    plot(x,y);
    hold on;
    y=(1/sqrt(2*pi*v2))*exp(-(x-u2).^2/2/v2);
    plot(x,y);
    str=['x=',num2str(Iw(i))];
    title(str);
end
A=[];
X8=zeros(8,1);
for i=1:m-7
    for j=1:n-7
        A1=double(I1(i:i+7,j:j+7))/255.0;
        A1=dct2(A1)';
        H4=reshape(A1,[64,1]);
        H=[1 2 6 7 15 16 28 29 3 5 8 14 17 27 30 43 4 9 13 18 26 31 42 44 10 12 19 25 32 41 45 54 11 20 24 33 40 46 53 55 21 23 34 39 47 52 56 61 22 35 38 48 51 57 60 62 36 37 49 50 58 59 63 64];
        X=zeros(64,1);
        X(H)=H4;
        for k1=1:8
            X8(k1)=X(I8(k1));
        end
        i1=-0.5*(X-uC)'*inv(vC)*(X-uC)-0.5*log((2*pi)^64*det(vC))+log(p1);
        i2=-0.5*(X-uG)'*inv(vG)*(X-uG)-0.5*log((2*pi)^64*det(vG))+log(p2);
        if(i1>i2)
            A(i,j)=255;
        else
            A(i,j)=0;
        end
        i11=-0.5*(X8-uc8)'*inv(vc8)*(X8-uc8)-0.5*log((2*pi)^8*det(vc8))+log(p1);
        i22=-0.5*(X8-ug8)'*inv(vg8)*(X8-ug8)-0.5*log((2*pi)^8*det(vg8))+log(p2);
        if(i11>i22)
            A8(i,j)=255;
        else
            A8(i,j)=0;
        end
    end
end
%imshow(A);
%figure;
%imshow(A8);
sum1=0;
sum2=0;
for i=1:m
    for j=1:n
        if I2(i,j)==255
            sum1=sum1+1;
        else sum2=sum2+1;
        end
    end
end
k1=0;
k2=0;
Ae=[];
 for i=1:m-7
    for j=1:n-7
        if(I2(i,j)~=A(i,j)&&I2(i,j)==255)
            k1=k1+1;
            Ae(i,j)=0;
        else if(I2(i,j)~=A(i,j)&&I2(i,j)==0)
                k2=k2+1;
                Ae(i,j)=0;
            else
                Ae(i,j)=1;
            end
        end
    end
end
noise1=k1/sum1*p1+k2/sum2*p2;   
k1=0;
k2=0;
Ae1=[];
for i=1:m-7
    for j=1:n-7
        if(I2(i,j)~=A8(i,j)&&I2(i,j)==255)
            k1=k1+1;
            Ae1(i,j)=0;
        else if(I2(i,j)~=A8(i,j)&&I2(i,j)==0)
                k2=k2+1;
                Ae1(i,j)=0;
            else
                Ae1(i,j)=255;
            end
        end
    end
end
noise2=k1/sum1*p1+k2/sum2*p2; 
figure; 
imshow(A);
figure;
imshow(A8);
figure;
imshow(Ae);
figure;
imshow(Ae1);