
function gE=gE(t,q1,q2,k,Cof,Fz,V,L,B,h1,E1,E2,MU)
%----BEAM DETAILS---------------------------------------------------------
P= 2373;

h2 = 0.2;
H=h1+h2;



Ar = (H*B);                     % Area




h0 = ((E1*(h1*h1))+(E2*((2*h1)+h2)*h2))/((2*E1*h1)+(2*E2*h2)); 
D1 = ((h1*((h1*h1*E1)-(3*h0*h1*E1)+(3*h0*h0*E1)))/(3))+((E2*(((h0-h1)^3)-((h0-h1-h2)^3)))/(3)); 



%----FOUDATION DETAILS----------------------------------------------------
K1= 100*(10^(6));                   % Linear stiffness
K3= 100*(10^(6));                   % Nonlinear stiffness
              % Viscous damping
GP = 6.66*10^7;            % Shear deformaton coefficient
 


     x=V*t+0.5*0.75*t^2;                        % Speed

b=((pi*k));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
gE=(1/P*Ar)*(-MU*q2-(K1+((GP*b*b)/L^2)+(D1*B*b*b*b*b)/L^4)*q1-(2*K3*Cof)+(2/L)*sin((b*x)/L)*Fz);
end

