clear all

h1 = 0.2;
h2 = 0.2;
H=h1+h2;
B=1;
L=50;                     




E1 = 19.942*(10^(9));
E2 = 0.4*(10^(9));
h0 = ((E1*(h1*h1))+(E2*((2*h1)+h2)*h2))/((2*E1*h1)+(2*E2*h2)); 
D1 = ((h1*((h1*h1*E1)-(3*h0*h1*E1)+(3*h0*h0*E1)))/(3))+((E2*(((h0-h1)^3)-((h0-h1-h2)^3)))/(3)); 
V1 = ((h0^(3))-((h0-h1)^(3)))/(3);

EV1  = 23172*(10^(6));                     %N/m^2
EV2  = 10730*(10^(6));                     %N/m^2
n1= 4313*(10^(6));                  %Ns/m^2
n2  = 2475*(10^(6));


%Viscoelastic parameters

p11 = (n1/EV1)+((n1+n2)/EV2);
p22 = ((n1*n2)/(EV1*EV2));
q11 = n1;
q22 = (n1*n2/EV2);

al = (p11+(((p11*p11)-(4*p22))^(1/2)))/(2*p22);
be = (p11-(((p11*p11)-(4*p22))^(1/2)))/(2*p22);

A = ((al*q22)-q11)/(((p11*p11)-(4*p22))^(1/2));
Be = (q11-(be*q22))/(((p11*p11)-(4*p22))^(1/2));


%n=50-----------------------------------------------------------------
p1=zeros(1,51);
KK = zeros(51,1);
 w1=zeros (51,1);
 q1= zeros(700000,51);
 q2= zeros(700000,51);
 v= 16.67; 
        % Speed
 for k = 1:1:50
    q1(1,k)=0;
    q2(1,k)=0;
       Tu1(1,k)=0;
    Tu2(1,k)=0;
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

for m = 1:1:70000
    t(m,1)= m*0.0001;
end

 h= 0.0001;       
for i = 1:1:15000
    
    for k = 1:1:50  
        
         m11 = -(2*al*Tu1(i,k))-(A*al*q1(i,k));
    m21 = -(2*be*Tu2(i,k))-(Be*be*q1(i,k));
    m12 = -(2*al*(Tu1(i,k)+(h*m11/2)))-(A*al*q1(i,k));
    m22 = -(2*be*(Tu2(i,k)+(h*m21/2)))-(Be*be*q1(i,k));
    m13 = -(2*al*(Tu1(i,k)+(h*m12/2)))-(A*al*q1(i,k));
    m23 = -(2*be*(Tu2(i,k)+(h*m22/2)))-(Be*be*q1(i,k));
    m14 = -(2*al*(Tu1(i,k)+(h*m13)))-(A*al*q1(i,k));
    m24 = -(2*be*(Tu2(i,k)+(h*m23)))-(Be*be*q1(i,k));
   
     Tu1(i+1,k) = Tu1(i,k)+(h*((m11+(2*(m12+m13))+m14)/6));
     Tu2(i+1,k) = Tu2(i,k)+(h*((m21+(2*(m22+m23))+m24)/6));
     
     
    ko = h*q2(i,k);
    lo = h*gE(t(i,1),q1(i,k),q2(i,k),k,Cof,Tu1(i,k),Tu2(i,k));
     
    k1 = h*(q2(i,k)+(lo/2));
    l1 = h*gE(t(i,1)+(h/2),q1(i,k)+(ko/2),q2(i,k)+(lo/2),k,Cof,Tu1(i,k),Tu2(i,k));

    k2 = h*(q2(i,k)+(l1/2));
    l2 = h*gE(t(i,1)+(h/2),q1(i,k)+(k1/2),q2(i,k)+(l1/2),k,Cof,Tu1(i,k),Tu2(i,k));
    
    k3 = h*(q2(i,k)+l2);
    l3 = h*gE(t(i,1)+h,q1(i,k)+k2,q2(i,k)+l2,k,Cof,Tu1(i,k),Tu2(i,k));
    
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

for h1=1:1:50
    for x = 0.1:0.1:50;
        e=round(10*x);
   w1(h1,e)=q1(i,h1)*sin(h1*pi*(x/L));
    end
end
S1=sum(w1);
W=transpose(S1);

f1 = figure;
figure(f1);
x = 0.1:0.1:50;
    e=round(10*x);
    x1=x;
     plot(x1,W(e,1),'-o');
    legend('Elastic');
title('the vertical deflection of the pavement');
xlabel('X(m)');
ylabel('Un(X,L/2v)(m)');