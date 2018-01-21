I=load('TrainingSamplesDCT_8_new.mat');
I1=imread('cheetah.bmp');
[m,n]=size(I1);
mask=imread('cheetah_mask.bmp');
sum1=0;
sum2=0;
[m1,n1]=size(I.TrainsampleDCT_FG);
[m2,n2]=size(I.TrainsampleDCT_BG);
p1=m1/(m1+m2);
p2=m2/(m1+m2);
for i=1:m
    for j=1:n
        if mask(i,j)==255
            sum1=sum1+1;
        else sum2=sum2+1;
        end
    end
end
bg=I.TrainsampleDCT_BG;
fg=I.TrainsampleDCT_FG;
[ub,sigmab,pib]=EMN(bg,32);
[uf,sigmaf,pif]=EMN(fg,32);
D=[1,2,4,8,16,24,32,40,48,56,64];
i1=6;
        A=[];
        for r=1:11
            d=D(r);
        for i=1:m-7
            for j=1:n-7
                A1=double(I1(i:i+7,j:j+7))/255.0;
                A1=dct2(A1)';
                H4=reshape(A1,[64,1]);
                H=[1 2 6 7 15 16 28 29 3 5 8 14 17 27 30 43 4 9 13 18 26 31 42 44 10 12 19 25 32 41 45 54 11 20 24 33 40 46 53 55 21 23 34 39 47 52 56 61 22 35 38 48 51 57 60 62 36 37 49 50 58 59 63 64];
                X=zeros(64,1);
                X(H)=H4;
                X1=X(1:d);
                u1=0;
                u2=0;
                for k=1:2^(i1-1)
                    uC=uf(k,1:d);
                    uG=ub(k,1:d);
                    vC=zeros(d);
                    vG=zeros(d);
                    piC=pif(k);
                    piG=pib(k);
                    for k1=1:d
                        vC(k1,k1)=sigmaf(k,k1);
                        vG(k1,k1)=sigmab(k,k1);
                    end
                    uC=uC';
                    uG=uG';
                    u1=u1+1/sqrt((2*pi)^d*det(vC))*exp(-0.5*(X1-uC)'*inv(vC)*(X1-uC))*pif(k);
                    u2=u2+1/sqrt((2*pi)^d*det(vG))*exp(-0.5*(X1-uG)'*inv(vG)*(X1-uG))*pib(k); 
                end
                u1=u1*p1;
                u2=u2*p2;
                if(u1>u2)
                    A(i,j)=255;
                else
                    A(i,j)=0;
                end
            end
        end
        k1=0;
        k2=0;
        for i=1:m-7
            for j=1:n-7
                if(mask(i,j)~=A(i,j)&&mask(i,j)==255)
                    k1=k1+1;
                elseif(mask(i,j)~=A(i,j)&&mask(i,j)==0)
                    k2=k2+1;
                end
            end
        end
        noise(r)=k1/sum1*p1+k2/sum2*p2; 
        end
        r1=1:11;
        plot(r1,noise);
        ylabel('PoE');
        title('D=32');
