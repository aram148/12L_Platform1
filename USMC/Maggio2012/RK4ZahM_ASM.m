%%   Date: 01/03/2021
% 
%    Author: Anand Rampadarath
%    Version 1
%    This model represents the Uterine smooth muscle cell as 
%    developed by Maggio et al 2012. It is a variant of the 
%    Hai-Murphy model including spatian dependence of the 
%    cross-bridges, similar to Mijailovich et al. 
%    We us the DM approximation to reduce the computational complexity.
% 

function [X0,T,R,rLum,force,stiffness,phosphorylation,m10,m20,nmp,nm,P_0,V,Raw,radius,Pab,vv,tau]=RK4ZahM_ASM(trange,N,M,xrange,lambda,rho,s,f,P01,gamma,Pawb)

%N number of time points
%trange specifies the interval of time steps
%T time grid
%X0 space grid
%gamma=25 rho=1 s=airway generation lambda=force activation

T=linspace(trange(1),trange(2),N);
X0=linspace(xrange(1),xrange(2),M);
dt=(trange(2)-trange(1))/(N); % dt must be <=0.1

%% Initial conditions f
   M10=0.309120293063990;
   M11=0.118325676501984;
   M12=0.066855519431672;
   M20=0.134379171029593;
   M21=-0.045132660622671;
   M22=0.101673063419951;
   C=0.6;

R=[M10;M11;M12;M20;M21;M22;C];

%% This run of fsolve finds the equilibria as Initial condition

% g=@(r) mis2(r,lambda,rho,N1(s),N2(s),P1(s),P2(s),rmax_sq(s),Ri_sq(s),P01);
% options = optimset('fsolve');
% [R,~,exitflag]=fsolve(g,R,options);
% R=[0.620322306169257; 0.310009711094628; 0.206570517906060; 0.174751423933002; 0.086923224037493; 0.057784269634018; 0.5];

stiffness = zeros(1,N);
force = zeros(1,N);
phosphorylation = zeros(1,N);
m10 = zeros(1,N);
m20 = zeros(1,N);
nmp = zeros(1,N);
nm = zeros(1,N);
nMpave = zeros(1,N);
nAmpave = zeros(1,N);
%[k1,k2] = rates(N,trange);
for i=1:N %loop over number of timesteps to update then run euler
    
    
    t=trange(1)+(i)*dt;
 
    [F1,~,~,~,~]=DM_funcs(t,R,lambda,f,rho,gamma);
    [F2,~,~,~,~]=DM_funcs(t+(dt/2),R+((dt/2)*F1),lambda,f,rho,N1(s),N2(s),Ri_sq(s),rmax_sq(s),rmax(s),P1(s),P2(s),gamma,Pawb,P01);%,t1,t2);%,k1,k2(i));
    [F3,~,~,~,~]=DM_funcs(t+(dt/2),R+((dt/2)*F2),lambda,f,rho,N1(s),N2(s),Ri_sq(s),rmax_sq(s),rmax(s),P1(s),P2(s),gamma,Pawb,P01);%,t1,t2);%,k1,k2(i));
    [F4,~,~,~,~,P_0(i),V(i),Raw(i),radius(i),Pab(i),vv(i),tau(i)]=DM_funcs(t+(dt),R+((dt)*F3),lambda,f,rho,N1(s),N2(s),Ri_sq(s),rmax_sq(s),rmax(s),P1(s),P2(s),gamma,Pawb,P01);%,t1,t2);%,k1,k2(i));
  
    Rnew = R + (dt/6)*(F1+(2*F2)+(2*F3)+F4);
    R=Rnew;
    
%% A rebuild of the distribution approximation for Namp(M10:M12) and Nam(M20:M22)

%     n1=(R(1)/((sqrt(2*pi))*q1)).*exp(-((X0-p1).^(2))/(2*(q1.^(2))));
%     n2=(R(4)/((sqrt(2*pi))*q2)).*exp(-((X0-p2).^(2))/(2*(q2.^(2))));
%%  ASM stiffness and radii of airway lumen  
    stiffness(i)=R(1)+R(4); 
    rLum(i)=R(8) ;
%% Crossbridge populations
    nmp(i)=1-R(7)-R(1);
    nm(i)=1-R(4)-R(1)-nmp(i);
    m10(i)=R(1);
%     m11(i)=R(2);
%     m21(i)=R(5);
    m20(i)=R(4);
    nMpave(i)=nmp(i);
    nAmpave(i)=R(1);
%%ASM phosphorylation and Force
    phosphorylation(i)=(nMpave(i)+nAmpave(i));
%% Emperical L-T relationship taken from Donovan 2016
    force_l(i)=(sin((pi*R(8))/(2*rmax(s)))).^3.*(R(8)<=2*rmax(s))+0.*(R(8)>2*rmax(s));
%     force_l(i)=1;
    force_a(i)=lambda*(R(2)+R(5));
    force(i)=force_l(i)*force_a(i);
%% plotting cross-bridge distributions
%                    if (mod(i,100)==0);% && i==8000 )
% %                         figure,
%                         plot(X0,n1+n2)
%                         xlim([-5 5])
% %                         title('Namp+Nam distribution plot for DM')
% %                         ylabel('n_{AMp}+n_{AM} state populations')
%                         drawnow
% %                       break
%                         hold on
% %                           pause
%     %     
%     %                      figure(2)
%     %                      plot(X0,n2)
%     %                      xlim([-5 5])
%     %                      title('N_{AM} distribution plot for DM')
%     %                      drawnow
%     %                      hold on
%     %                      pause
%     %
%     %                     figure(3)
%     %                     plot(X0,n1)
%     %                     title('Namp distribution for DM')
%     %                     drawnow
%     %                     hold on
%     
%     %                     pause
%                   end
  
    
end