I=load("TrainingSamplesDCT_8.mat");
[m n]=size(I.TrainsampleDCT_FG);
[x y]=size(I.TrainsampleDCT_BG);
M1=I.TrainsampleDCT_FG;
M2=I.TrainsampleDCT_BG;
Pf=m/(m+x);
Pb=x/(x+m);
H1=zeros(1,64);
H2=zeros(1,64);
for i=1:m
    [t index]=sort(M1(i,:));
    k=index(63);
    H1(k)=H1(k)+1;
end  
for i=1:x
    [t index]=sort(M2(i,:));
    k=index(63);
    H2(k)=H2(k)+1;
end
H3=zeros(1,64);
for i=1:64
    x1=H1(i)*Pf/double(m);
    x2=H2(i)*Pb/double(x);
    if(x1>x2)
        H3(i)=255;
    else
        H3(i)=0;
    end
end
bar(H3);
I1=imread('cheetah.bmp');
I1=double(I1)/255;
[m n]=size(I1);
A=[];
for i=1:m-7
    for j=1:n-7
        A1=I1(i:i+7,j:j+7);
        A1=abs(dct2(A1));
        H4=reshape(A1,[1,64]);
        H=[1 2 6 7 15 16 28 29 3 5 8 14 17 27 30 43 4 9 13 18 26 31 42 44 10 12 19 25 32 41 45 54 11 20 24 33 40 46 53 55 21 23 34 39 47 52 56 61 22 35 38 48 51 57 60 62 36 37 49 50 58 59 63 64];
       [t, index]=sort(H4);
       k=index(63);
       k=H(k);
       k=H3(k);
       A(i,j)=k;
    end
end
figure;
A1=imagesc(A);

A2=colormap(gray(255));
I2=imread('cheetah_mask.bmp');
k1=0;
k2=0;
sum1=0;
sum2=0;
for i=1:m-7
    for j=1:n-7
        if(I2(i,j)~=A(i,j)&&I2(i,j)==255)
            k1=k1+1;
        else if(I2(i,j)~=A(i,j)&&I2(i,j)==0)
                k2=k2+1;
            end
        end
        if(I2(i,j)==255)
            sum1=sum1+1;
        else
            sum2=sum2+1;
         end
    end
end
p1=sum1/(sum1+sum2);
p2=sum2/(sum1+sum2);
noise=k1/sum1*p1+k2/sum2*p2;

        
                        
                        
        