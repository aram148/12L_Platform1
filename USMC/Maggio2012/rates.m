function [k1,k2,t] = rates(N,trange)

T=linspace(trange(1),trange(2),N);
%define constants
k1a = 0.5962;
k1b = 1.35;
tau = 156.9;
k2b = 242.14;
kon1 = 0.00125;
kon2 = 0.8988;
koff1 = 0.4629;
koff2 = 20.035;
P0 = 0;

c = 0.1.*(T<1/0.01)+0.65.*(T>=1/0.01 & T<=10/0.01)+1.3.*(T>10/0.01 & T<=18/0.01)+ 0.65.*(T>18/0.01);
a = 0.*(T<22/0.01)+6.*(T>=22/0.01 & T<30/0.01)+0.*(T>=30);

kon = kon1 + (c.^2/(kon2^2+c.^2));
koff = koff1 + (a*koff2)/(1+a);


% F = (1/tau)*kon*(1-P)-koff*P;

[t,P] = ode15s(@(t,P) (1/tau)*kon*(1-P)-koff*P, T, P0);
k2 = k2b*P.^2;
k1 = k1a*c.^4/(k1b^4+c.^4);
k1=k1';
k2=k2';