function [ub1,sigmab1,pib1] = EM(bg)
[m1,n1]=size(bg);
u=zeros(8,64);
sigma=zeros(8,64);
for i=1:8
    for j=1:64
        u(i,j)=5*rand-5*rand;
    end
end
for i=1:8
    sigma(i,1:64)=10+(5*rand-5*rand);
end
pi0=zeros(8,1);
o=0;
for i=1:7
    pi0(i)=1/7*rand;
    o=o+pi0(i,1);
end
pi0(8)=1-o;
p=0;
while(1)
    un=zeros(8,64);
    sigman=zeros(8,64);
    pin=zeros(8,1);
    for i=1:8
        H=0;
        H1=0;
        H2=0;
        v1=zeros(64);
        for k=1:64
            v1(k,k)=sigma(i,k);
        end
        u1=u(i,1:64);
        for k=1:m1
            x=bg(k,1:64);
            s=1/sqrt(det(v1))*exp(-1/2*(x-u1)*inv(v1)*(x-u1)')*pi0(i);
            s0=0;
            for n=1:8
                v0=zeros(64);
                for y=1:64
                    v0(y,y)=sigma(n,y);
                end
                u0=u(n,1:64);
                s0=s0+1/sqrt(det(v0))*exp(-1/2*(x-u0)*inv(v0)*(x-u0)')*pi0(n);
            end
            h1=s/s0;
            H=H+h1;
            H1=H1+h1.*x;
            H2=H2+h1.*(x-u1).*(x-u1);
        end
        un(i,1:64)=H1./H;
        pin(i)=H./m1;
        sigman(i,1:64)=H2./H;
        for k=1:64
            if(sigman(i,k)<0.0001)
                sigman(i,k)=0.0001;
            end
        end
    end
    p=0;
    for i=1:64
        p=p+(u(i)-un(i))^2;
    end
    p=p/64;
    if(p<0.00005)
        break;
    end
    u=un;
    pi0=pin;
    sigma=sigman;
end
ub1=u;
sigmab1=sigma;
pib1=pi0;
end

