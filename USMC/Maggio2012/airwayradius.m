function drdt=airwayradius(t,r,P0,kappa)

%% Parameter values used for both ASM-aw and Anafi-Wilson coupled system
Ri_sq=0.0174;
rmax_sq=0.1980;
N1=1;
N2=8;
P1=0.4021;
P2=-33.3394;
rmax=0.4450;
rho=1;
E=254.9;%cmH2O/ml %Elastance of acinus
L = 2.6852;%mm %airway length
% L=2*pi*rmax;
aa = 10.19716213; % Conversion factor for kg*(mm)^(-1)*(s)^(-2) to cmH20
moo = 1.9008e-8; % dynamic viscosity kg(mm)^(-1)s^(-1)

f=0.25;
%  %% Anafi-Wilson type pressure difference
%     P01=P011;
%    
%     Pawbb = Pawb*abs(sin(2*pi*f));
% 
%     Raw = 8*aa*L*moo/(pi*r.^4);
% 
%     Paw = P01 + Pawbb;%pressure in airway
% 
%     Pab = Pawb*E/(E.^2 + (2*pi*f*Raw).^2)^0.5;
% 
%     alph = atan(2*pi*f*Raw/E);
% 
%     Pabb = Pab*abs(sin(2*pi*f*t-alph));%.*(t<a)+((10*Pab))*abs(sin(2*pi*f*t-alph)).*(t>=a).* (t<=b)...
% 
%     Pa = P01 + Pabb; %pressure in acinus
%     vv = 0.2 +0.04*Pa; %Acinus volume
%     xx = 1-(r(8)/vv.^(1/3));
%     P_lumen = (Paw + Pa)/2;
%     %P0=P_lumen + Pa;
%     tau = Pa + Pa*(1.4*xx+2.1*xx.^2);%parenchymal tethering force

 %Transmural pressure as a function of force.
% Ptm=Pmin-((force*r_ref)./r(8))+tau;
    Ptm=((P0)-(kappa./r));

if Ptm<=0
    % a=a0(s).*(1-(Ptm./P1)).^(-N1(s));
    rad=sqrt(Ri_sq*(1-(Ptm/P1)).^-N1);
    
elseif Ptm>=0
    % a=(1-((1-a0(s)).*(1-(Ptm./P2)).^(-N2(s))));
    rad=sqrt(rmax_sq-(rmax_sq-Ri_sq)*(1-(Ptm/P2)).^-N2);
    
end

drdt=rho*(rad-(r));  % this is to be used for the orginal DM coupled system.


