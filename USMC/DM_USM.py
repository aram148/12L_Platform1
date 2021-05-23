import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
# import opencor as oc
# simulation = oc.open_simulation(r"C:\Users\aram148\Desktop\Platform1\EP3\USMC\ca2.cellml")

def DM_funcs_USM(t,R):
    
    
    # data = simulation.data()
    # data.set_starting_point(t)
    # data.set_ending_point(trange2)
    # data.set_point_interval(1)
    # simulation.run()
    # results=simulation.results()
    # P = results.states()['calcium/P'].values() 
    f = 0.25
    # if t<=32:
    #     C_cai = 0*np.sin(2*np.pi*f*t)
    # elif t>32 and t< 34:
    #     C_cai = 10*np.sin(2*np.pi*f*t)
    # else: 
    #     C_cai= 0*np.sin(2*np.pi*f*t)*(t>=34)
    C_cai = 0*np.sin(2*np.pi*f*t)*(t<=32)+3*np.sin(2*np.pi*f*t)*(t>=32 and t<=34)+0*np.sin(2*np.pi*f*t)*(t>34 and t<=60)+\
    3*np.sin(2*np.pi*f*t)*(t>=60 and t<=62)+0*np.sin(2*np.pi*f*t)*(t>62)  
    n = 8.7613
    Ca_mlck = 256.98
    k1 = (C_cai**n)/((Ca_mlck**n)+(C_cai**n))

    k2 = 1.2387
    g1 = 0.0756
    gp1 = 0.0709
    fp1 = 0.2838
    
    r = R
    gam = 100
    
    # Mean for cdf
    p1 = r[1]/r[0]
    p2 = r[4]/r[3]
    
    # Std dev for cdf
    q1 = np.sqrt((r[2]/r[0])-(r[1]/r[0])**2)
    q2 = np.sqrt((r[5]/r[3])-(r[4]/r[3])**2)
    
    # r, phi and I for M1_lambda
    r0 = -p1/q1
    r1 = (1-p1)/q1
    
    phi0 = 0.5*(1+math.erf((r0-p1)/(q1*np.sqrt(2))))
    phi1 = 0.5*(1+math.erf((r1-p1)/(q1*np.sqrt(2))))
    
    I0 = -(np.exp(-((-p1/q1)**2)/2))/(np.sqrt(2*np.pi))
    I1 = -(np.exp(-((((1-p1)/q1))**2)/2))/(np.sqrt(2*np.pi));
    
    # r, phi and I for M2_lambda
    r20 = -p2/q2
    r21 = (1-p2)/q2
    
    phi20 = 0.5*(1+math.erf((r20-p2)/(q2*np.sqrt(2))))
    phi21 = 0.5*(1+math.erf((r21-p2)/(q2*np.sqrt(2))))
    
    I20 = -(np.exp(-((-p2/q2)**2)/2))/(np.sqrt(2*np.pi))
    I21 = -(np.exp(-((((1-p2)/q2))**2)/2))/(np.sqrt(2*np.pi))
    
    # Functions for the RHS of M1 PDE
    J0 = phi0
    J10 = ((p1*phi0) + (q1*I0))
    J11 = (p1*phi1) + (q1*I1)
    J20=((p1**(2))*phi0)+((2*p1*q1)*I0)+((q1**(2))*(phi0+(r0*I0)))
    J21=((p1**(2))*phi1)+((2*p1*q1)*I1)+((q1**(2))*(phi1+(r1*I1)))
    J30=(p1**(3)*phi0)+((3*p1**(2)*q1)*I0)+((3*p1*q1**(2))*((phi0)+(r0*I0)))+((q1**(3)*(2+r0**(2))*I0))
    J31=(p1**(3)*phi1)+((3*p1**(2)*q1)*I1)+((3*p1*q1**(2))*(phi1+(r1*I1)))+(q1**(3)*(2+(r1**(2)))*I1)
    
    #Functions defined for the RHS of the second PDE M2_lambda
    K0 = phi20
    K01 = phi21
    K10 = (p2*phi20) + (q2*I20)
    K11 = (p2*phi21) +( q2*I21)
    K20=((p2**(2))*phi20) + ((2*p2*q2)*I20) + ((q2**(2))*(phi20 + (r20*I20)))
    K21=((p2**(2))*phi21) + ((2*p2*q2)*I21) + ((q2**(2))*(phi21 + (r21*I21)))
    K30=(p2**(3)*phi20) + ((3*p2**(2)*q2)*I20) + ((3*p2*q2**(2))*((phi20) + (r20*I20))) + ((q2**(3)*(2+r20**(2))*I20))
    K31=(p2**(3)*phi21) + ((3*p2**(2)*q2)*I21) + ((3*p2*q2**(2))*(phi21 + (r21*I21))) + (q2**(3)*(2+(r21**(2)))*I21)
    
    # Components for the matrix F that will represent each moment,
    # M1_lambda and M2_lambda
  
    A0 = ((fp1*(1-r[6]))/1)-(fp1*(J11-J10)*r[0])
    A1 = ((fp1*(1-r[6]))/2)-(fp1*(J21-J20)*r[0])
    A2 = ((fp1*(1-r[6]))/3)-(fp1*(J31-J30)*r[0])
    
    B0=(3*(fp1+gp1)*J0)+gp1*(J11-J10)+(4*gp1)*(p1-J11)
    B1=(3*(fp1+gp1)*J10)+gp1*(J21-J20)+(4*gp1)*((p1**(2)+q1**(2))-J21)
    B2=(3*(fp1+gp1)*J20)+gp1*(J31-J30)+(4*gp1)*((p1**(3)+3*p1*q1**(2))-J31)
    
    C0=k1
    C1=k1*p2
    C2=k1*(p2**(2)+q2**(2))
    
    D0=k2
    D1=k2*p1
    D2=k2*(p1**(2)+q1**(2))
    
    E0=(20*g1*K0)+g1*(K11-K10)+(g1)*(1-K01)
    E1=(20*g1*K10)+g1*(K21-K20)+(g1)*((p2)-K11)
    E2=(20*g1*K20)+g1*(K31-K30)+(g1)*(p2**(2)+q2**(2)-K21)
    
    V = gam*(A1-E0-B0)/(1+gam*(D0+C0))
    
    F=np.array([A0-B0*r[0]+C0*r[3]-k2*r[0], A1-B1*r[0]+C1*r[3]-k2*r[1]-V*r[0], A2-B2*r[0]+C2*r[3]-k2*r[2]-2*V*r[1], D0*r[0]-E0*r[3]-k1*r[3], D1*r[0]-E1*r[3]-k1*r[4]-V*r[3], D2*r[0]-E2*r[3]-k1*r[5]-2*V*r[4], -k1*r[6]+(1-r[6])*k2])
    return F, C_cai
    
def USM_DM(trange,N,M,xrange):
    T = np.linspace(trange[0], trange[1], N)
    X0 = np.linspace(xrange[0], xrange[1], M)
    dt=(trange[1]-trange[0])/(N)
    M10 = 0.309120293063990
    M11 = 0.118325676501984
    M12 = 0.066855519431672
    M20 = 0.134379171029593
    M21 = -0.045132660622671
    M22 = 0.101673063419951
    C = 0.6
 #   R = np.zeros(7)
 # Calculating the Ca2+ for the same trange
    # simulation = oc.open_simulation(r"C:\Users\aram148\Desktop\Platform1\EP3\USMC\ca2.cellml")
    # data = simulation.data()
    # data.set_starting_point(trange[0])
    # data.set_ending_point(trange[1])
    # data.set_point_interval(dt)
    # simulation.run()
    # results=simulation.results()
    # C_cai = results.states()['calcium/P'].values() 
    C_cai = np.zeros([len(T),1])
    R = np.array([M10, M11, M12, M20, M21, M22, C])
    Rnew = np.zeros([len(T),len(R)])
    Rnew[0,:] = R
    # i = 0

    
    for i in range(1,len(T)):
        t = T[i-1]
        F1, C_cai[i] =DM_funcs_USM(t,Rnew[i-1,:])
        F2, C_cai[i]=DM_funcs_USM(t+(dt/2),Rnew[i-1,:]+(dt/2)*F1)
        F3, C_cai[i]=DM_funcs_USM(t+(dt/2),Rnew[i-1,:]+(dt/2)*F2)
        F4, C_cai[i]=DM_funcs_USM(t+(dt),Rnew[i-1,:]+(dt)*F3)
        Rnew[i,:] = Rnew[i-1,:]+ (dt/6)*(F1 + (2*F2) + (2*F3) + F4)
        # Rnew[i,:] = Rnew[i-1,:] + dt*F1
    return Rnew, C_cai
        
        
        
trange = [0, 100]
xrange = [-3, 3]
N = 5000
M = 1000
kappa = 5

T = np.linspace(trange[0], trange[1], N)
Rnew, C_cai = USM_DM(trange, N, M, xrange)
force = (Rnew[:,1]+Rnew[:,4])
stiffness = Rnew[:,0]+Rnew[:,3]
nmp = 1-Rnew[:,6]-Rnew[:,0]
phosphorylation = nmp + Rnew[:,0]
plt.plot(T,force)
# plt.plot(T,C_cai)
plt.show()




    
    
        
    


