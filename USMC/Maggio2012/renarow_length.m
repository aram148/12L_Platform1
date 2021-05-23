function [F_DI2,F_min2,F11,F22,Bron_dil2,T,Force]=renarow_length(M,xrange,f,lambda,rho,s,P01,gamma,Pmin,t1,t2)
% clear all

trange2=t1+60;
%b=a+2;
a1=60;
b1=62;
N=trange2*100;
r0=0.005;
%close all
 format short
% a=40;
% b=42;
% a1=20;
% b1=22;
% [~,T,~,Force,~,~,~,~,~,~,~,~,~]=RK4HHM1_ASM_rad([0 trange2],N,M,xrange,f,lambda,rho,s,P01,gamma,Pmin,r0,a,b);
% [~,T,~,~,Force,~,~,~,~,~,~,~,~]=RK4HHM1_ASM([0 trange2],N,M,xrange,f,lambda,rho,s,P01,gamma,Pmin,a,b);
[~,T,~,~,Force,~,~,~,~,~,~,~,~,Raw]=RK4ZahM_ASM([0 trange2],N,M,xrange,lambda,rho,s,f,P01,gamma,Pmin,t1,t2);
% [~,T,~,Force,~,~,~,~,~,~,~,~]=RK4ZahM_ASM_rad([0 trange2],N,M,xrange,lambda,rho,s,f,P01,gamma,Pmin,r0,a,b);

[Maxima,MaxIdx] = findpeaks(Force);
% [Max_rad,Max_rad_Idx]=findpeaks(rLum);
% [Max_P(i,j),Max_P_Idx(i,j)]=findpeaks(P0);
DataInv = 1.01*max(Force) - Force;
[Minima_F,MinIdx] = findpeaks(DataInv);
Minima_F = Force(MinIdx);
%% FOR 2 DI AFTER EQUILIBRIA
M=round(length(Maxima)/3);
[val1,ind1]=max(Maxima(1:M));
[val2,ind2]=max(Maxima(M+1:end));
ind2=ind2+M;
[val3,ind3]=min(Minima_F(round(M/2):M));

[F_min2,ind4]=min(Minima_F(M+1:end));
ind4=ind4+M;
%% FOR 2 DI, ONE AT START
% [val1,ind1]=max(Maxima);
Force1_bDI=a1/0.01;
Force1_aDI=b1/0.01;
Force_bDI=(t1/0.01);
Force_aDI=t2/0.01;
% time_axis=1:length(Force(Force_bDI:end));
% % plot(60.*time_axis./time_axis(end),Force(Force_bDI:end));
% plot((T(Force_bDI:end)),(Force(Force_bDI:end)));

% 
%  figure(1)
% plot((T(Force_bDI:end)-T(Force_bDI)),(Force(Force_bDI:end)));

% figure,
% plot((T(Force1_bDI:Force_bDI)-T(Force1_bDI)),(Force(Force1_bDI:Force_bDI))./Force(Force1_bDI));

%% Average force during DI

F_DI1=trapz(T(Force1_bDI:Force1_aDI),Force(Force1_bDI:Force1_aDI));
F_DI2=trapz(T(Force_bDI:Force_aDI),Force(Force_bDI:Force_aDI));

%% F_{min} post DI

% F_min1=Force(Force1_aDI);
% F_min2=Force(Force_aDI);

% DataInv_r = 1.01*max(rLum) - rLum;
% [Minima_r,MinIdx_r] = findpeaks(DataInv_r);
% Minima_r = rLum(MinIdx_r);
% [val1,ind1]=max(Max_rad);

% DataInv_P = 1.01*max(P0) - P0;
% [Minima_P,MinIdx_P] = findpeaks(DataInv_P);
% Minima_P = P0(MinIdx_P);
% [val2,ind2]=max(Max_P);


% Force(MinIdx(ind))
% T(MinIdx(ind))
F1=trapz(T(MinIdx(ind1-2):MinIdx(ind1-1)),Force(MinIdx(ind1-2):MinIdx(ind1-1)))/(T(MinIdx(ind1-1))-T(MinIdx(ind1-2)));%Force before 1st DI
F2=trapz(T(MinIdx(ind1):MinIdx(ind1+1)),Force(MinIdx(ind1):MinIdx(ind1+1)))/(T(MinIdx(ind1+1))-T(MinIdx(ind1)));%Force after 1st DI
% F3=trapz(T(MinIdx(ind1-1):MinIdx(ind1)),Force(MinIdx(ind1-1):MinIdx(ind1)))/T(MinIdx(ind1)-MinIdx(ind1-1));%Force during DI
% F4=trapz(T(MinIdx(end-1):MinIdx(end)),Force(MinIdx(end-1):MinIdx(end)))/(T((MinIdx(end)-MinIdx(end-1))));%Max Force before 2nd DI
 Bron_dil1=abs(F1-F2);
%  Del_F1=F4-F2;
% % %     end
% % % end
% % 
% force_rec= Maxima(ind1+1:end)./Bron_dil1;
% % norm_force_rec=force_rec./force_rec(1);
% norm_force_rec=(force_rec-min(force_rec))./((Maxima(ind1-1)./Bron_dil1)-min(force_rec));
% time_axis_step=1:length(force_rec);
% time_axis=time_axis_step.*(60./length(force_rec));
% % plot(time_axis,Maxima(ind2+1:end),'o-')
%  plot(time_axis,(norm_force_rec).*100,'o-');
 %% 

% figure,
% plot(time_axis,(norm_force_rec).*100,'o-');
% title('post 1st DI')
%% for second DI
F11=trapz(T(MinIdx(ind2-2):MinIdx(ind2-1)),Force(MinIdx(ind2-2):MinIdx(ind2-1)))/(T(MinIdx(ind2-1))-T(MinIdx(ind2-2)));%Force pre DI
F22=trapz(T(MinIdx(ind2):MinIdx(ind2+1)),Force(MinIdx(ind2):MinIdx(ind2+1)))/(T(MinIdx(ind2+1))-T(MinIdx(ind2)));%Force post DI
% F33=trapz(T(MinIdx(ind2-1):MinIdx(ind2)),Force(MinIdx(ind2-1):MinIdx(ind2)))/T(MinIdx(ind2)-MinIdx(ind2-1));%Force during DI
% F44=trapz(T(MinIdx(end-1):MinIdx(end)),Force(MinIdx(end-1):MinIdx(end)))/(T((MinIdx(end)-MinIdx(end-1))));%Max Force recovered after 2mins 
 Bron_dil2=abs(F11-F22);
%  Del_F2=F44-F22;
% 

%% Maximum Force post 2nd DI

[Fmax,indx]=findpeaks((Force(Force_aDI:end)));
% max_force=findpeaks((100.*Force(Force_bDI:end))./Force(Force_bDI));
force_rec1=Fmax./Bron_dil2;
norm_force_rec1=(force_rec1-min(force_rec1))./((Maxima(ind2-1)./Bron_dil2)-min(force_rec1));
%  figure(2)
plot(T(indx),norm_force_rec1.*100,'-o');

% force_rec1=Maxima(ind2-1:end)./Bron_dil2;
% % norm_force_rec1=force_rec1./force_rec1(end);
% norm_force_rec1=(force_rec1-min(force_rec1))./((Maxima(ind2-1)./Bron_dil2)-min(force_rec1));
% time_axis_step1=1:length(force_rec1);
% time_axis1=time_axis_step1.*(60./length(force_rec1));
% % % figure,
% % % plot(time_axis1,(norm_force_rec1).*100,'o-');
% plot(time_axis1,Maxima(ind2-1:end),'o-')
% 
% 

% title('post 2nd DI')
% plot(120.*T(MinIdx(ind2-1):MinIdx(end))./(T(MinIdx(end))),Force(MinIdx(ind2-1):MinIdx(end)));
% plot(120.*T(MinIdx(ind2-1):MinIdx(end))./(T(MinIdx(end))),Force(MinIdx(ind2-1):MinIdx(end)));
% figure,
% plot(Max_P_Idx(ind2+1:end),Max_P(ind2+1:end))
% figure,
% plot(Max_rad(ind1+1:end),'-o')
% figure,
% plot(Maxima(ind+2:end),'-*')
% figure,
% plot(Force((MinIdx(ind+1):MinIdx(end))))
% figure,
% plot(T,rLum./rmax)
% ylabel('normalised renarrowing of airway radius')
% xlabel('Time(s)')
% figure,
% plot(T,Force)
% ylabel('Induced radius  for DI')
% xlabel('Time(s)')