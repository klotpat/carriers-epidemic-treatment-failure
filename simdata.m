function  [T_allend]= simdata(beta,f,Time,P_ext_ana)
% Parameters 
N= 1000;
alpha= 0.85; 
gamma= 0.07; 
xi= 0.19089; 
thc= 0.012;  
thi= 0.75; 
eta= 0.1; 
sig= 0.5; 
psi=0.01; 
var=0.15; 
rng(42); % random seed
%basic reproduction number
R0= beta*N*(alpha*eta*f*thi+var*eta*f*thi-alpha*eta*gamma-alpha*eta*thi-alpha*f*thc+alpha*gamma*xi-var*f*thc-eta*f*thi+alpha*thc+eta*gamma+eta*thi+f*thc+psi)/(var*f*gamma*xi*thi-var*f*gamma*thc-f*gamma*xi*thi-f*psi*thi-f*thi*thc+gamma*gamma*xi+gamma*xi*thi+gamma*psi+gamma*thc+psi*thi+thc*thi);

count=0;
time_cut_off = 0;
numpath=10000;
T_allend = [];

for k=1:numpath
    clear t s i c a
    t(1)=0; 
    s(1)=999; 
    i(1)=0; 
    c(1)=1; 
    a(1)=0;
    j=1; 
    
    force = beta*(i(j)+eta*c(j));
    while (i(j)>0 || c(j)>0 || a(j)>0) && t(j)<= Time 
        u1=rand; 
        force = beta*(i(j)+eta*c(j));
        lambda = force*s(j)+sig*a(j)+(gamma+thi)*i(j)+(gamma*xi+thc+psi)*c(j);
        t(j+1)  = -(log(u1)/lambda)+t(j);
        p1 = (1-alpha)*force*s(j)/lambda;
        p2 = alpha*force*s(j)/lambda;
        p3 = thc*c(j)/lambda;
        p4 = psi*c(j)/lambda;
        p5 = thi*i(j)/lambda;
        p6 = (sig*(1-f)*a(j))/lambda;
        p7 = (var*sig*f*a(j))/lambda;
        p8 = (1-var)*sig*f*a(j)/lambda;
        p9 = (gamma*i(j))/lambda;
        p10 = (gamma*xi*c(j))/lambda;
        
% updating event
        u2 = rand;
        if u2 <= p1
            s(j+1)= s(j)-1;
            i(j+1)= i(j);
            c(j+1)= c(j)+1;
            a(j+1)= a(j);
        elseif u2 > p1 && u2 <= p1+p2
            s(j+1)= s(j)-1;
            i(j+1)= i(j)+1;
            c(j+1)= c(j);
            a(j+1)= a(j);
        elseif u2 > p1+p2 && u2<=p1+p2+p3
            s(j+1)= s(j);
            i(j+1)= i(j);
            c(j+1)= c(j)-1;
            a(j+1)= a(j)+1;
        elseif u2 > p1+p2+p3 && u2<=p1+p2+p3+p4 %difference from model1
            s(j+1)= s(j);
            i(j+1)= i(j)+1;
            c(j+1)= c(j)-1;
            a(j+1)= a(j);
        elseif u2> p1+p2+p3+p4 && u2<=p1+p2+p3+p4+p5
            s(j+1)= s(j);
            i(j+1)= i(j)-1;
            c(j+1)= c(j);
            a(j+1)= a(j)+1;
        elseif u2>p1+p2+p3+p4+p5 && u2<=p1+p2+p3+p4+p5+p6
            s(j+1)= s(j)+1;
            i(j+1)= i(j);
            c(j+1)= c(j);
            a(j+1)= a(j)-1;
        elseif (u2>p1+p2+p3+p4+p5+p6) && (u2<=p1+p2+p3+p4+p5+p6+p7) % difference from model 1
            s(j+1)= s(j);
            i(j+1)= i(j);
            c(j+1)= c(j)+1;
            a(j+1)= a(j)-1;
        elseif (u2>p1+p2+p3+p4+p5+p6+p7) && (u2<=p1+p2+p3+p4+p5+p6+p7+p8)
            s(j+1)= s(j);
            i(j+1)= i(j)+1;
            c(j+1)= c(j);
            a(j+1)= a(j)-1;
        elseif (u2>p1+p2+p3+p4+p5+p6+p7+p8) && (u2<=p1+p2+p3+p4+p5+p6+p7+p8+p9)
            s(j+1)= s(j)+1;
            i(j+1)= i(j)-1;
            c(j+1)= c(j);
            a(j+1)= a(j);
        elseif (u2>p1+p2+p3+p4+p5+p6+p7+p8+p9) && (u2<=p1+p2+p3+p4+p5+p6+p7+p8+p9+p10)
            s(j+1)= s(j)+1;
            i(j+1)= i(j);
            c(j+1)= c(j)-1;
            a(j+1)= a(j);
        else
            s(j+1)= s(j);
            i(j+1)= i(j);
            c(j+1)= c(j);
            a(j+1)= a(j);
        end
        if i(j)<0 
            i(j)=0;
        end
        if c(j)<0
            c(j)=0;
        end
        if a(j)<0
            a(j)=0;
        end
        j=j+1;
    end
  if i(end)<=0 && c(end)<=0 && a(end)<=0
     count=count+1;
     T_allend = [T_allend t(end)];
  else
      T_allend = [T_allend Time];
      time_cut_off = time_cut_off+1;
  end
end
R0
P_ext_ana
P_ext_sim = count/numpath
error = P_ext_ana-P_ext_sim
time_cut_off
end






















