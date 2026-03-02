function pdfpaper2
%--------------------------
%Parameters
N= 1000;
beta= 0.000775; 
alpha= 0.8585; 
gamma= 0.07; %recover fast 2days 
xi= 0.189; %about 80 days
thc= 0.012;  
thi= 0.75; 
eta= 0.1; %10percent of sympomatic
f= 0.2; %error 10percent
sig= 0.5; 
psi=0.01; 
var=0.15; %related with thtaC

 
%Time Step
Time=30000; h=0.01;  n=Time/h;
 
%Initial Condition
x(1,1)=0;
y(1,1)=0;
z(1,1)=0;

%Euler Algorithm
for i=1:n
    Dx(1,i) = (beta*eta*N+xi*gamma+psi+thc)*(((1-alpha)*beta*eta*N*x(1,i)^2+alpha*beta*eta*N*x(1,i)*y(1,i)+xi*gamma+psi*y(1,i)+thc*z(1,i))/(beta*eta*N+xi*gamma+psi+thc)-x(1,i));
    Dy(1,i) = (beta*N+gamma+thi)*(((1-alpha)*beta*N*x(1,i)*y(1,i)+alpha*beta*N*y(1,i)^2+gamma+thi*z(1,i))/(beta*N+gamma+thi)-y(1,i));
    Dz(1,i) = sig*((1-f)+var*f*x(1,i)+(1-var)*f*y(1,i)-z(1,i));
     
    x(1,i+1)= x(1,i)+h*Dx(1,i);
    y(1,i+1)= y(1,i)+h*Dy(1,i);
    z(1,i+1)= z(1,i)+h*Dz(1,i);
   
end
pdfx = (beta*eta*N+xi*gamma+psi+thc)*(((1-alpha)*beta*eta*N*x.^2+alpha*beta*eta*N*x.*y+xi*gamma+psi*y+thc*z)/(beta*eta*N+xi*gamma+psi+thc)-x); 
pdfy = (beta*N+gamma+thi)*(((1-alpha)*beta*N*x.*y+alpha*beta*N*y.^2+gamma+thi*z)/(beta*N+gamma+thi)-y);
pdfz = sig*((1-f)+var*f*x+(1-var)*f*y-z);

format long
P_ext_c = round(max(x),8)
P_ext_i = round(max(y),8)

format


  t=0:h:Time;
  plot(t,pdfx,'LineWidth',1)
% plot(t,1-y,'LineWidth',1)
%  xlabel('{\it t} (day)')
%  ylabel('Probability that the epidemic lasts at least {\it t} days')
 


 

