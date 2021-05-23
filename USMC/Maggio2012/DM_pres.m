function P0=DM_pres(P01,f,t,t1,t2,Pmin,a,b)

P0=P01*sin(2*pi*f*t).*(t<a)+((7*P01))*sin(2*pi*f*t).*(t>=a).* (t<=b)+P01*sin(2*pi*f*t).*(t>b).*(t<t1)+((7*P01))*sin(2*pi*f*t).*(t>=t1).*(t<=t2)+(P01)*sin(2*pi*f*t).*(t>t2)+Pmin;

% P0=((7*P01))*sin(2*pi*f*t).*(t>=a).* (t<=b)+P01*sin(2*pi*f*t).*(t>b).*(t<t1)+((7*P01))*sin(2*pi*f*t).*(t>=t1).*(t<=t2)+(P01)*sin(2*pi*f*t).*(t>t2)+Pmin;

% P0=((7*P01))*sin(2*pi*f*t).*(t>=a).* (t<=b)+0.*(t>b).*(t<t1)+((7*P01))*sin(2*pi*f*t).*(t>=t1).*(t<=t2)+0.*(t>t2)+Pmin;

% P0=0.*(t<a)+((7*P01))*sin(2*pi*f*t).*(t>=a).* (t<=b)+0.*(t>b).*(t<t1)+((7*P01))*sin(2*pi*f*t).*(t>=t1).*(t<=t2)+0.*(t>t2)+Pmin;