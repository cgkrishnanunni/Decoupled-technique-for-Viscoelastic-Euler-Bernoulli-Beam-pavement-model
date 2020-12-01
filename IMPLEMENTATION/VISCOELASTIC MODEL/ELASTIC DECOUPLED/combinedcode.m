clear all

%VEHICLE PROPERTIES
% -----------------------------------------------------------------

Mb=21260; % Body Mass in Kg
Ms=190;   % Axle Mass
Kb=2.06*10^6 ; %Body Stiffness
Cb= 3*10^4; % Body Damping
V=16.67;       % Vehicle speed

h1 = 0.2;
h2 = 0.2;



B=1;
L=50;        

MU = 0.3*(10^(6));  

E1 = 19942*(10^(6));
E2= 0.4*(10^(9));

%n=50-----------------------------------------------------------------
kj=1;


 deltat1=0.0001;
timepoints=round((L/V)/deltat1);

       
        
p1=zeros(1,51);
KK = zeros(51,1);
 w1=zeros (51,1);
 q1= zeros(timepoints,51);
 q2= zeros(timepoints,51);

 

 L0=10;
 a=0.0002;
        % Speed
        
        x1=zeros(timepoints,1);
          r=zeros(timepoints,1);
          e=zeros(timepoints,1);
        for i=1:(timepoints)      %roughness generation
  x1(i)=V*deltat1*(i-1);
   r(i,1)=a*sin((2*pi*x1(i))/L0);
end
        
        
        wb=zeros(timepoints,1);

wbnew=zeros(timepoints,1);


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

for i=1:timepoints
    if i+1<=timepoints
    Ust(i)=((wb(i+1)-wb(i))/deltat1)+((r(i+1)-r(i))/deltat1);
    end
    if i-1~=0 & i+1<=timepoints
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
        
        Fz(i)=((Mb+Ms)*9.81-(Mb*Uvehiclett+Ms*Ustt(i)));
         
         
         t1(i)=timevehicle;
         r1(i)=Uvehicle(1);
         r2(i)=Uvehiclett(1);
         r3(i)=Uvehiclet(1);
         
         timevehicle=timevehicle+ deltat1;
         
end

%Galerkin truncation
        
 for k = 1:1:50
    q1(1,k)=0;
    q2(1,k)=0;
end
Cof =0;
 for i=1:15
    xp(i,1)=0.5*(1-cos(((i-1)*pi)/(15-1)));
end

for i=1:15
    l(i,1)=1/i;
end

for i=1:15
    for j=1:15
        
        X(i,j)=(xp(j,1))^(i-1);
    end
end

I=linsolve(X,l);


for m = 1:1:timepoints
    t(m,1)= m*0.0001;
end

 h= 0.0001;       
for i = 1:1:timepoints
    
    for k = 1:1:50  
    ko = h*q2(i,k);
    lo = h*gE(t(i,1),q1(i,k),q2(i,k),k,Cof,Fz(i),V,L,B,h1,E1,E2,MU);
     
    k1 = h*(q2(i,k)+(lo/2));
    l1 = h*gE(t(i,1)+(h/2),q1(i,k)+(ko/2),q2(i,k)+(lo/2),k,Cof,Fz(i),V,L,B,h1,E1,E2,MU);

    k2 = h*(q2(i,k)+(l1/2));
    l2 = h*gE(t(i,1)+(h/2),q1(i,k)+(k1/2),q2(i,k)+(l1/2),k,Cof,Fz(i),V,L,B,h1,E1,E2,MU);
    
    k3 = h*(q2(i,k)+l2);
    l3 = h*gE(t(i,1)+h,q1(i,k)+k2,q2(i,k)+l2,k,Cof,Fz(i),V,L,B,h1,E1,E2,MU);
    
    q1(i+1,k)  = q1(i,k)+((ko+(2*(k1+k2))+k3)/6);
    q2(i+1,k)  = q2(i,k)+((lo+(2*(l1+l2))+l3)/6);
    

    KK(k+1,1)= q1(i,k);
    p1(1,k)= q1(i,k);
    
    
    end
     
    for j =1:1:15
      
        for k11 =1:1:50
            
            FJ(k11,j) = KK(k11,1)*sin(k11*pi*xp(j));
            
        end
        
    end
      
    Fj = transpose(sum(FJ));
    
    for j=1:1:15
        
       COF(j,1) = (Fj(j,1)^(3))*I(j,1)*sin(j*pi*xp(j));
    end
    
    
    
    Cof = sum(COF);
     
        
end

%for h1=1:1:50
   % for x = 0.1:0.1:50;
      %  e=round(10*x);
       
   %w1(h1,e)=q1(15000,h1)*sin(h1*pi*(x/L));
   
    %end
%end
%S1=sum(w1);
%W=transpose(S1);


       


%rat=1;
 %for x = 0.1:0.1:50;
 %value(rat)=x;
%rat=rat+1;
% end
 

time11=0;
       for i=1:timepoints
    wbnew(i)=0;
    for j=1:50
        
      wbnew(i)=wbnew(i)+q1(i,j)*sin((j*pi*V*time11)/L);
      
   
    end
   wbnew(i)=wbnew(i);
    time11=time11+deltat1;

       end
end

for i=1:timepoints
    w=0;
    for j=1:50
        
      w=w+q1(i,j)*sin(j*pi*0.5);
      
        
   end
    
   e(i)=w;   
end


%f1 = figure;
%figure(f1);
%x = 0.1:0.1:50;
   % e=round(10*x);
   % x1=x;
    % plot(x1,W(e,1),'-o');
    %legend('Elastic');
%title('the vertical deflection of the pavement');
%xlabel('X(m)');
%ylabel('Un(X,L/2v)(m)');