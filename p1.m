D=load('TrainingSamplesDCT_subsets_8.mat');
I1=imread('cheetah.bmp');
mask=imread('cheetah_mask.bmp');
[m,n]=size(I1);
I0=load('Prior_2.mat');
alpha=load('Alpha.mat');
alpha=alpha.alpha;
mu0B=I0.mu0_BG;
mu0F=I0.mu0_FG;
w0=I0.W0;
D1B=D.D3_BG;
D1F=D.D3_FG;
covB=zeros(64);
covF=zeros(64);
[m1,n1]=size(D1B);
[m2,n2]=size(D1F);
pb=m1/(m1+m2);
pf=m2/(m1+m2);
B1=zeros(1,64);
F1=zeros(1,64);
for i=1:m1
    B1=B1+D1B(i,1:64);
end
B1=B1/m1;
for i=1:m1
    covB=covB+(D1B(i,1:64)-B1)'*(D1B(i,1:64)-B1);
end
covB=covB/m1;
for i=1:m2
    F1=F1+D1F(i,1:64);
end
F1=F1/m2;
for i=1:m2
    covF=covF+(D1F(i,1:64)-F1)'*(D1F(i,1:64)-F1);
end
covF=covF/m2;
cov0=zeros(64);
for i=1:64
    cov0(i,i)=w0(i);
end
R1=[];
k=1;
for i=1:m-7
    for j=1:n-7
        A1=double(I1(i:i+7,j:j+7))/255.0;
        A1=dct2(A1)';
        H4=reshape(A1,[64,1]);
        H=[1 2 6 7 15 16 28 29 3 5 8 14 17 27 30 43 4 9 13 18 26 31 42 44 10 12 19 25 32 41 45 54 11 20 24 33 40 46 53 55 21 23 34 39 47 52 56 61 22 35 38 48 51 57 60 62 36 37 49 50 58 59 63 64];
        X=zeros(64,1);
        X(H)=H4;
        R1(1:64,k)=X;
        k=k+1;
    end
end
sum1=0;
sum2=0;
for i=1:m
    for j=1:n
        if mask(i,j)==255
            sum1=sum1+1;
        else
            sum2=sum2+1;
        end
    end
end
MAP=[];
ML=[];
PD=[];
r=[];
[M,N]=size(R1);
for i=1:9
    alpha0=alpha(i);
    cov00=alpha0.*cov0;
    u1F=cov00*(cov00+1/m2.*covF)^(-1)*F1'+1/m2.*covF*(cov00+1/m2.*covF)^(-1)*mu0F';
    u1B=cov00*(cov00+1/m1.*covB)^(-1)*B1'+1/m1.*covB*(cov00+1/m1.*covB)^(-1)*mu0B';
    cov1F=cov00*(cov00+1/m2.*covF)^(-1)*covF./m2;
    cov1B=cov00*(cov00+1/m1.*covB)^(-1)*covB./m1;
    unF=u1F;
    covnF=cov1F+covF;
    unB=u1B;
    covnB=cov1B+covB;
    k=1;
    e1=0;
    e2=0;
    for j=1:m-7
        for h=1:n-7
            X=R1(1:64,k);
            i1=-0.5*(X-B1')'*inv(covB)*(X-B1')-0.5*log((2*pi)^64*det(covB))+log(pb);
            i2=-0.5*(X-F1')'*inv(covF)*(X-F1')-0.5*log((2*pi)^64*det(covF))+log(pf);
            if(i1>i2)
                A(j ,h)=0;
            else
                A(j,h)=255;
            end
            if(uint8(A(j,h))==255&&mask(j,h)==0)
                e1=e1+1;
            end
            if(uint8(A(j,h))==0&&mask(j,h)==255)
                e2=e2+1;
            end
            k=k+1;
        end
    end
    ML(i)=e1/sum2*pb+e2/sum1*pf;
    k=1;
    e1=0;
    e2=0;
    for j=1:m-7
        for h=1:n-7
            X=R1(1:64,k);
            i1=-0.5*(X-u1B)'*inv(covB)*(X-u1B)-0.5*log((2*pi)^64*det(covB))+log(pb);
            i2=-0.5*(X-u1F)'*inv(covF)*(X-u1F)-0.5*log((2*pi)^64*det(covF))+log(pf);
            if(i1>i2)
                A(j ,h)=0;
            else
                A(j,h)=255;
            end
            if(uint8(A(j,h))==255&&mask(j,h)==0)
                e1=e1+1;
            end
            if(uint8(A(j,h))==0&&mask(j,h)==255)
                e2=e2+1;
            end
            k=k+1;
        end
    end
    MAP(i)=e1/sum2*pb+e2/sum1*pf;
    k=1;
    e1=0;
    e2=0;
    for j=1:m-7
        for h=1:n-7
            X=R1(1:64,k);
            i1=-0.5*(X-u1B)'*inv(covnB)*(X-u1B)-0.5*log((2*pi)^64*det(covnB))+log(pb);
            i2=-0.5*(X-u1F)'*inv(covnF)*(X-u1F)-0.5*log((2*pi)^64*det(covnF))+log(pf);
            if(i1>i2)
                A(j,h)=0;
            else
                A(j,h)=255;
            end
            if(uint8(A(j,h))==255&&mask(j,h)==0)
                e1=e1+1;
            end
            if(uint8(A(j,h))==0&&mask(j,h)==255)
                e2=e2+1;
            end
            k=k+1;
        end
    end
    PD(i)=e1/sum2*pb+e2/sum1*pf;
end
 set(gca,'XScale', 'log');
 semilogx(alpha,ML);
 hold on;
 semilogx(alpha,MAP);
 hold on;
 semilogx(alpha,PD);
 xlabel('\alpha');
 ylabel('PoE');
 title('Strategty2, D3');