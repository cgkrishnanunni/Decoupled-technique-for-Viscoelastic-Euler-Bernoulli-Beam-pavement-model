%VEHICLE PAVEMENT UNCOUPLED SYSTEM USING AN INTERATIVE PROCEDURE FOR CONVERGENCE%
%Programmed by Krishnanunni C G at IIT MADRAS                          %
%Programming dates: March 2018                                       %
%Problem: Refer Sec 4.3 Dynamics of Vehicle road coupled system 
% ----------------------------------------------------------------- 



%VEHICLE PROPERTIES
% -----------------------------------------------------------------

Mb=21260; % Body Mass in Kg
Ms=190;   % Axle Mass
Kb=2.06*10^6 ; %Body Stiffness
Cb= 3*10^4; % Body Damping
V1=15;       % Vehicle speed



YY=6.988e9;% Young's modulus
pp=2373; % density
LL=160; %length
L0=5;
k1=97.552;%dimensionless
k3=2.497e6; %dimensionless value
mu=39.263;
V=0.00874 ;%dimensionles
GG=3.2e9;
alfa1=0.1829;
Area=0.3;


beta=6.25e5;


Gp=0.0318;
kf=1.626e3;
cf=2.618e4;
a=0.01;



timepoints=18000; %No of time points

deltat=0.0042918;



deltat1=(deltat*LL)/(sqrt(YY/pp));



r=zeros(timepoints+1,1);

%for i=1:(timepoints+1)      %roughness generation
 %  x1(i)=V1*deltat1*(i-1);
 %  r(i)=a*sin((2*pi*x1(i))/L0);
%end

for i=1:(timepoints+1)
    
r(i)=ACLASS(i);
end

%VEHICLE EQUATION SOLVER

wb=zeros(timepoints+1,1);

wbnew=zeros(timepoints+1,1);


while ((((max(wbnew)-max(wb))/max(wbnew)))>0.01 || max(wbnew)==0)
    
        
    wb=wbnew;
    max(wb)

 delta = 0.5; 
	alfa = 0.25;
	a2 = 1.0/alfa/deltat1;
	a1 = delta*a2;
	a0 = a2/deltat1;
	a4 = delta/alfa - 1;
	a3 = 0.5/alfa - 1;
	a5 = deltat1/2.0*(delta/alfa - 2);
	a6 = deltat1*(1 - delta);
	a7 = delta * deltat1;
    
    Uvehicle=zeros(1,1);
    Uvehiclet=zeros(1,1);
    Uvehiclett=zeros(1,1);
    newUvehicle=zeros(1,1);
    newUvehiclet=zeros(1,1);
    newUvehiclett=zeros(1,1);

    Qvehicle=Kb+a0*Mb+a1*Cb;
    Qvehicle=inv(Qvehicle);
     timevehicle=0;
    

Us=zeros(timepoints+1,1);
Ust=zeros(timepoints+1,1);
Ustt=zeros(timepoints+1,1);
Fz=zeros(timepoints+1,1);
Us=wb+r;

for i=1:timepoints+1
    if i+1<=timepoints+1
    Ust(i)=((wb(i+1)-wb(i))/deltat1)+((r(i+1)-r(i))/deltat1);
    end
    if i-1~=0 & i+1<=timepoints+1
    Ustt(i)=((wb(i+1)-2*wb(i)+wb(i-1))/deltat1^2)+((r(i+1)-2*r(i)+r(i-1))/deltat1^2);
    end
    Fvehicle=[Cb*Ust(i)+Kb*Us(i)];
    
    Fvehiclenew=zeros(1,1);
        Fvehiclenew = Fvehicle + Mb*(a0*Uvehicle + a2*Uvehiclet + a3*Uvehiclett)+ Cb*(a1*Uvehicle + a4*Uvehiclet + a5*Uvehiclett);
    
    newUvehicle = Qvehicle * Fvehiclenew;
		
		newUvehiclett = a0*(newUvehicle - Uvehicle) - a2*Uvehiclet - a3*Uvehiclett;
		newUvehiclet  = Uvehiclet + a6*Uvehiclett + a7*newUvehiclett;
		
		Uvehicle   = newUvehicle;
		Uvehiclet  = newUvehiclet;
		Uvehiclett = newUvehiclett;
        
        Fz(i)=((Mb+Ms)*9.81-(Mb*Uvehiclett+Ms*Ustt(i)))/(Area*YY);
         
         
         t1(i)=timevehicle;
         r1(i)=Uvehicle(1);
         r2(i)=Uvehiclett(1);
         r3(i)=Uvehiclet(1);
         
         timevehicle=timevehicle+ deltat1;
         
end

% GALERKIN TRUNCATION START

N=100;% GALERKIN TRUNCATION POINTS
n=13;

q=zeros(timepoints,N);
L=zeros(timepoints,N);
X=zeros(n);
l=zeros(n,1);
I=zeros(n,1);
x=zeros(n,1);

for i=1:n
    x(i,1)=0.5*(1-cos(((i-1)*pi)/(n-1)));
end

for i=1:n
    l(i,1)=1/i;
end
for i=1:n
    for j=1:n
        
        X(i,j)=(x(j,1))^(i-1);
    end
end
       I=linsolve(X,l);
       t=0;
      for j=1:N
         q(1,j)=0;
          q(2,j)=q(1,j);
      end
      
       for j=1:N
         L(1,j)=0;
          L(2,j)=L(1,j);
       end
      
       for i=1:timepoints
         
           for o=1:N
               P=0;
           for j=1:n
               R=0;
               for p=1:N
           R=R+q(i,p)*sin(p*pi*x(j));
               end
               R=R^3;
           P=P+(I(j)*R*sin(o*pi*x(j)));
           end
           
           
          if (i-1)~=0
              ll=(mu/deltat)+(1/(deltat)^2);
          q(i+1,o)=(1/ll)*(((mu*q(i,o))/deltat)-((q(i-1,o)-2*q(i,o))/(deltat)^2)-(k1+Gp*(o*pi)^2+alfa1*(o*pi)^2)*q(i,o)+(2*Fz(i)*sin(o*pi*V*t))- P*2*k3+ L(i,o)*alfa1*(o*pi));
          
          end
          
           if (i-1)~=0
              lll=(cf/deltat)+(1/(deltat)^2);
          L(i+1,o)=(1/lll)*(((cf*L(i,o))/deltat)-((L(i-1,o)-2*L(i,o))/(deltat)^2)-(kf+(o*pi)^2+beta)*L(i,o)+beta*(o*pi)*q(i,o));
          
          end
           end
       
           t=0+i*deltat;  
       end 
       
       time11=0;
       for i=1:timepoints+1
    wbnew(i)=0;
    for j=1:N
        
      wbnew(i)=wbnew(i)+q(i,j)*sin(j*pi*V*time11);
      
   
    end
   wbnew(i)=wbnew(i)*LL;
    time11=time11+deltat;
   
    
end


end % while loop

w=0;time=0;
 
 k=zeros(1000,1);
      g=zeros(1000,1);
       e=zeros(timepoints,1);
       r123=zeros(timepoints,1);
       time=0;e=0;w=0;

for i=1:timepoints
    w=0;
    for j=1:N
        
      w=w+q(i,j)*sin(j*pi*0.5);
       time=((i*deltat)*LL)/(sqrt(YY/pp));
        
   end
    
    e(i)=w*LL;
    r123(i)=time;
    
end
       
        xx=zeros(1000,1);
       for i=1:1000
          xx(i)=0.5*(1-cos(((i-1)*pi)/(999)));
       end
       
       for i=1:1000
         w=0;
      for j=1:N
           
       w=w+q(13333,j)*sin(j*pi*xx(i));
       end
      
       k(i,1)=w*LL;
       g(i,1)=xx(i)*LL;
       
      end

     

