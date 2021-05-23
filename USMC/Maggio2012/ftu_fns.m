function [odes,Ptm,Plumen,Paw] = ftu_fns(t,R)

M=R(1); AM=R(2); Mp=R(3); AMp=R(4);
r=R(5); PA=R(6); VA=R(7); fo=R(8);
fc=R(9);po=R(10); pc=R(11); z=R(12);
dPtau=R(13);

    k2 = 0.5; %s^-1
    k7 = 0.01; %s^-1
    k3 = 0.4; %s^-1
    k4 = 0.1; %s^-1
    k5 = 0.5; %s^-1
%     Ri_sq = 0.0174; %mm^2
    Ri_sq = 0.000174; %cm^2
%     rmax_sq = 0.1980; %mm^2 %
    rmax_sq = 0.00198; %cm^2;
    N1 = 1; %dimensionless
    N2 = 8; %dimensionless
%     P1 = 0.295768328; %mmHg
    P1 = 0.40209999959928; %cmH2O
%     P2 = -24.52310023; %mmHg
    P2 = -33.339399993687; %cmH2O
    rho = 1; %dimensionless
%     rmax = 0.4450; %mm
    rmax = 0.04450; %cm^2
    P0 = 10; %cmH20
    f = 0.25; %s^-1
%     Pmin = 10;
%     PL = 2;
    kappa = 1; %dimensionless
    aa = 10.19716213; % Conversion factor for kg*(mm)^(-1)*(s)^(-2) to cmH20
    moo = 1.9008e-8;  % dynamic viscosity kg(mm)^(-1)s^(-1)
%     E = 186; %mmHg/ml
    E = 252.8688649029042; %cmH2O/ml
%     L = 2.6;
%     Do = 0.0000156; %mols^(-1)mmHg^(-1)
    Do = 1.1475e-05; %mols^(-1)cmH2O(-1)
%     Dc = 0.0000316; %mols^(-1)mmHg^(-1)
    Dc = 2.3244e-05; %mols^(-1)cmH2O(-1)
    fom = 0.21; %dimensionless
    fcm = 0; %dimensionless
%     VD = 0.151; %litres
    VD = 151; %ml
%     VT = 0.41; %litres
    VT = 410; %ml
%     L = 2.6; %mm
    L = 0.26; %cm
    KT = 10000; %dimensionless
    Kc = 3600000; %dimensionless
%     sigma = 0.0000014; %mol L^-1 mmHg
    sigma = 1.9033e-06; %mol L^-1 cmH2O
%     Vc = 0.071; %litre
    Vc = 71; %ml
    Th = 2e-6; %mol l^-1
    l2 = 164000; %dimensionless
    r2 = 0.12; %s^-1
%     sigmac = 0.000033; %mol l^-1 mmHg
    sigmac = 4.4864e-05; %mol L^-1 cmH2O
%     z = 0.00000044219; %dimensionless
    delta = 10^1.9; %dimensionless
    hin =  1.00000.*1^-9; %dimensionless
%     A=5e5;
%     
%     mu = abs(A/(3*r.^4))/2;
%     
    Paw =  P0.* sin( 2.00000.* pi.*f.*t);
    k1=0.35.*(t<5)+0.06.*(t>=5);
    k6=k1;
    df_satdp = ( ( L.*power(1.00000+ KT.*sigma.*po, 4.00000)+power(1.00000+ Kc.*sigma.*po, 4.00000)).*( 3.00000.*L.*power(KT, 2.00000).*power(sigma, 2.00000).*po.*1.00000.*power(1.00000+ KT.*sigma.*po, 2.00000)+ L.*KT.*sigma.*1.00000.*power(1.00000+ KT.*sigma.*po, 3.00000)+ 3.00000.*power(Kc, 2.00000).*power(sigma, 2.00000).*po.*1.00000.*power(1.00000+ Kc.*sigma.*po, 2.00000)+ Kc.*sigma.*1.00000.*power(1.00000+ Kc.*sigma.*po, 3.00000)) -  ( L.*KT.*sigma.*po.*power(1.00000+ KT.*sigma.*po, 3.00000)+ Kc.*sigma.*po.*power(1.00000+ Kc.*sigma.*po, 3.00000)).*( 4.00000.*L.*KT.*sigma.*1.00000.*power(1.00000+ KT.*sigma.*po, 3.00000)+ 4.00000.*Kc.*sigma.*1.00000.*power(1.00000+ Kc.*sigma.*po, 3.00000)))./power( L.*power(1.00000+ KT.*sigma.*po, 4.00000)+power(1.00000+ Kc.*sigma.*po, 4.00000), 2.00000);
%     v = 0.200000+ 0.00400000.*PA;
%     x = 1.00000 - r./power(v, 0.330000);
%     Ptau = 2*mu*(PA+ PA.*( 1.40000.*x+ 2.10000.*power(x, 2.00000)));
%     Ptau = 2*mu.*((rmax - r)/rmax + 1.5*((rmax - r)/rmax).^2);

    
    Raw = ( 8.00000.*aa.*L.*moo)./(  pi.*power(r, 4.00000));
    Ptau = P0 -(0.5.*Raw.*f.*VT).*sin(f.*t)-E.*(2.5-0.5.*VT.*cos(f.*t));
    Pab = ( P0.*E)./power((power(E, 2.00000)+power( 2.00000.* pi.*f.*Raw, 2.00000)), 1.0 ./ 2);
    stress = AMp+AM;
    fL = power( sin( ((  pi.*r)./2.00000).*rmax), 3.00000);
    Plumen = (Paw+PA)./2.00000;
    Ptm = (Plumen - ( fL*kappa.*stress)./r)+Ptau;
    pac =  fc.*(PA - Paw);
    q = (Pab - PA)./Raw;
    pao =  fo.*(PA - Paw);
    QA = q+ 1.00000.*Dc.*(pc - pac)+ 1.00000.*Do.*(po - pao);
%     phos = AMp+Mp;
 
%     alpha = atan(( 2.00000.* pi.*f.*t.*Raw)./E);
%     P0 = Plumen - PA;
    
    if Ptm<=0
    % a=a0(s).*(1-(Ptm./P1)).^(-N1(s));
    rad=sqrt(Ri_sq*(1-(Ptm/P1)).^-N1);
    
    elseif Ptm>=0
    % a=(1-((1-a0(s)).*(1-(Ptm./P2)).^(-N2(s))));
    rad=sqrt(rmax_sq-(rmax_sq-Ri_sq)*(1-(Ptm/P2)).^-N2);
    
    end
%    radius=rad;
   
    if VT>=VD
        foi=(fo.*VD+ fom.*(VT - VD))./VT;
    elseif VT<VD
        foi=fo;
    end

    if VT>=VD
        fci=(fc.*VD+ fcm.*(VT - VD))./VT;
    elseif VT<VD
        fci=fc;
    end
% dPtau = PA.*( 4.20000.*x.*dxdt+ 1.40000.*dxdt)+( Pab.*E.*QA)./( PA.*1.00000)+ dPtau.*( 1.40000.*x+ 2.10000.*(x).^2);
%        =  2*mu*(-(3*(rmax-r)*rho*(rad-r))/rmax.^2-(rho*(rad-r)/rmax))
% dPtau = -0.5*Raw*f^2*VT*cos(f*t)-E*(2.5-0.5*VT*sin(f*t));

% dxdt = ((( 0.00132000.*r.*Pab.*E.*QA)./( PA.*1.00000)+dPtau)./( 0.00400000.*PA+0.200000).^1.33000 - rdot./( 0.00400000.*PA+0.200000).^0.330000);

% dPtau = 2*mu(-(3*(rmax-r)*rho*(rad-r))/rmax.^2-(rho*(rad-r)/rmax));

odes =[ - k1.*M+ k2.*Mp+ k7.*AM; k5.*AMp -  (k6+k7).*AM; ( k4.*AMp+ k1.*M) -  (k2+k3).*Mp; ( k3.*Mp+ k6.*AM) -  (k5+k4).*AMp;
    rho.*(rad - r); ( Pab.*E.*QA)./PA+ dPtau; ((Pab - Ptau) -  VA.*E)./Raw; (1.00000./VA).*(( 1.00000.*Do.*(po - pao)+ (foi - fo).*q) -  fo.*( 1.00000.*Dc.*(pc - pac)+ 1.00000.*Do.*(po - pao)));
    (1.00000./VA).*(( 1.00000.*Dc.*(pc - pac)+ (fci - fc).*q) -  fc.*( 1.00000.*Do.*(po - pao)+ 1.00000.*Dc.*(pc - pac)));
    (Do./( sigma.*Vc)).*power(1.00000+ (( 4.00000.*Th)./sigma).*df_satdp,  - 1.00000).*( fo.*(PA - Paw) - po);
    ( (Dc./( sigmac.*Vc)).*(pac - pc)+ (( 1.00000.*delta.*l2)./sigmac).*hin.*z) -  delta.*r2.*pc;( delta.*r2.*sigmac.*pc)./1.00000 -  delta.*l2.*hin.*z;
   -0.5.*Raw.*f^2.*VT.*cos(f.*t)-E.*(2.5-0.5.*VT.*sin(f.*t))];



