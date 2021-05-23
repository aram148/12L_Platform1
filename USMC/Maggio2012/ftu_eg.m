function [r, PA, PAW]=ftu_eg()

%Parameters
%     k2 = 0.5;
%     k7 = 0.01;
%     k3 = 0.4;
%     k4 = 0.1;
%     k5 = 0.5;
%     Ri_sq = 0.0174;
%     rmax_sq = 0.1980;
%     N1 = 1;
%     N2 = 8;
%     P1 = 0.295768328;
%     P2 = -24.52310023;
%     rho = 1;
%     rmax = 0.4450;
%     P0 = 10;
%     f = 0.25;
%     Pmin = 10;
%     PL = 2;
%     kappa = 1;
%     aa = 10.19716213;
%     moo = 1.9008e-8;
%     E = 186;
%     L = 2.6;
%     Do = 0.0000156;
%     Dc = 0.0000316;
%     fom = 0.21;
%     fcm = 0;
%     VD = 0.151;
%     VT = 0.41;
%     L = 2.6;
%     KT = 10000;
%     Kc = 3600000;
%     sigma = 0.0000014;
%     Vc = 0.071;
%     Th = 2e-6;
%     l2 = 164000;
%     r2 = 0.12;
%     sigmac = 0.000033;
%     z = 0.00000044219;
%     delta = power(10.0000, 1.90000);
%     hin =  1.00000.*power(1.00000,  - 9.00000);
%  
% %Initial conditions
% 
%     M = 1.0;
%     AM = 0.0;
%     Mp = 0.0;
%     AMp = 0.0;
%     r = 0.4450;
%     PA = 10;
%     VA = 1;
%     fo = 0.1368;
%     fc = 0.05263;
%     po = 25;
%     pc = 28;
%     
%     R=[M,AM,Mp,AMp,r,PA,VA,fo,fc,po,pc];
% 
% %Algebraic Terms
%     k1=0.35.*(t<5)+0.06.*(t>=5);
%     k6=k1;
%     df_satdp = ( ( L.*power(1.00000+ KT.*sigma.*po, 4.00000)+power(1.00000+ Kc.*sigma.*po, 4.00000)).*( 3.00000.*L.*power(KT, 2.00000).*power(sigma, 2.00000).*po.*1.00000.*power(1.00000+ KT.*sigma.*po, 2.00000)+ L.*KT.*sigma.*1.00000.*power(1.00000+ KT.*sigma.*po, 3.00000)+ 3.00000.*power(Kc, 2.00000).*power(sigma, 2.00000).*po.*1.00000.*power(1.00000+ Kc.*sigma.*po, 2.00000)+ Kc.*sigma.*1.00000.*power(1.00000+ Kc.*sigma.*po, 3.00000)) -  ( L.*KT.*sigma.*po.*power(1.00000+ KT.*sigma.*po, 3.00000)+ Kc.*sigma.*po.*power(1.00000+ Kc.*sigma.*po, 3.00000)).*( 4.00000.*L.*KT.*sigma.*1.00000.*power(1.00000+ KT.*sigma.*po, 3.00000)+ 4.00000.*Kc.*sigma.*1.00000.*power(1.00000+ Kc.*sigma.*po, 3.00000)))./power( L.*power(1.00000+ KT.*sigma.*po, 4.00000)+power(1.00000+ Kc.*sigma.*po, 4.00000), 2.00000);
%     v = 0.200000+ 0.00400000.*PA;
%     x = 1.00000 - r./power(v, 0.330000);
%     Ptau = PA+ PA.*( 1.40000.*x+ 2.10000.*power(x, 2.00000));
%     Raw = ( 8.00000.*aa.*L.*moo)./(  pi.*power(r, 4.00000));
%     Pab = ( P0.*E)./power((power(E, 2.00000)+power( 2.00000.* pi.*f.*Raw, 2.00000)), 1.0 ./ 2);
%     stress = AMp+AM;
%     Plumen = (Paw+PA)./2.00000;
%     Ptm = (Plumen - ( kappa.*stress)./r)+Ptau;
%     pac =  fc.*(PA - Paw);
%     q = (Pab - PA)./Raw;
%     pao =  fo.*(PA - Paw);
%     QA = q+ 1.00000.*Dc.*(pc - pac)+ 1.00000.*Do.*(po - pao);
%     phos = AMp+Mp;
%     fL = power( sin( ((  pi.*r)./2.00000).*rmax), 3.00000);
%     alpha = atan(( 2.00000.* pi.*f.*VOI.*Raw)./E);
%     P0 = Plumen - PA;
%     
%     if Ptm<=0
%     % a=a0(s).*(1-(Ptm./P1)).^(-N1(s));
%     rad=sqrt(Ri_sq*(1-(Ptm/P1)).^-N1);
%     
%     elseif Ptm>=0
%     % a=(1-((1-a0(s)).*(1-(Ptm./P2)).^(-N2(s))));
%     rad=sqrt(rmax_sq-(rmax_sq-Ri_sq)*(1-(Ptm/P2)).^-N2);
%     
%     end
%    radius=rad;
%    
%     if VT>=VD
%         foi=(fo.*VD+ fom.*(VT - VD))./VT;
%     elseif VT<VD
%         foi=fo;
%     end
% 
%     if VT>=VD
%         fci=(fc.*VD+ fcm.*(VT - VD))./VT;
%     elseif VT<VD
%         fci=fc;
%     end
% Paw =  P0.* sin( 2.00000.* pi.*f.*t);
% dxdt = ((( 0.00132000.*r.*Pab.*E.*QA)./( PA.*1.00000)+dPtau)./( 0.00400000.*PA+0.200000).^1.33000 - rdot./( 0.00400000.*PA+0.200000).^0.330000);
% dPtau =( PA.*( 4.20000.*x.*dxdt+ 1.40000.*dxdt)+( Pab.*E.*QA)./( PA.*1.00000)+ dPtau.*( 1.40000.*x+ 2.10000.*(x).^2);
 