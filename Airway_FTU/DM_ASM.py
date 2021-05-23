# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 20:02:48 2021

@author: aram148
"""
import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
import pandas as pd

A = pd.ExcelFile(r"airways1.xlsx")
sheet1 = A.parse(0)
Ri_sq = sheet1['Ri ^2 (mm^2)']
rmax_sq = sheet1['Rimax^2 (mm^2)']
N1 = sheet1['n 1']
N2 = sheet1['n 2']
P1 = sheet1['P1 (cm H2O)']
P2 = sheet1['P2 (cm H2O)']
rmax = sheet1['r imax (mm)']
Rm = sheet1['R m (mm)']
Rw = sheet1['R w (mm)']
Ri = sheet1['R i (mm)']
rsmax = sheet1['r smax(mm)']

def DM_funcs(t,R,lam,s):
    
    c = 0.1*(t<10)+(0.055*(t-10)+0.1)*(t>=10 and t<20)+ 0.65*(t>=20 and t<100)+0.1*(t>=100)
    a = 0*(t<=50)+ 6*(t>=50 and t<=75)+ 0*(t>75)
    k_on = 0.000301
    k_ia = 0.5962
    k_ib = 1.35 #uM
    k_2b = 76.23
    k_off1 = 0.4629
    k_off2 = 20.035
    k1 = k_ia*(c**4/(k_ib**4 + c**4))
    k_off = k_off1+(k_off2*a/(a+1))
    k2 = k_2b*(k_on/(k_on+k_off))
    # k1=0.35*(t<5)+0.06*(t>=5)
    # k1 = 0.06
    k2=0.1
    k7=0.005
    g1=2.*k7
    g2=20.*g1
    g3=3.*g1
    fp1=0.88
    gp1=0.22
    gp2=4.*(fp1+gp1)
    gp3=3.*gp1
    rho = 1
    f = 0.25
    P_t = 5*np.sin(2*np.pi*f*t)
    Pmin = 10+ P_t
    gamma = 25
    nu = 1.5
    moo = 15.09
    r = R
    
    #ASM force
    force_a=lam*(r[1]+r[4])
    
    
    Ptm = Pmin - (force_a/r[7])
  
   
    
    if Ptm<=0:
        rad = np.sqrt(Ri_sq[s-1]*(1-(Ptm/P1[s-1]))**-N1[s-1])
    elif Ptm>=0:
        rad = np.sqrt(rmax_sq[s-1]-(rmax_sq[s-1]-Ri_sq[s-1])*(1-(Ptm/P2[s-1]))**-N2[s-1])
    else:
        print('There is an error, probably with your erf fn')
    
    # Thickness of SMC layer and airway wall
    e_w = (Rw[s-1]-Ri[s-1])/Ri[s-1]
    e_m = (Rm[s-1]-Rw[s-1])/Ri[s-1]
    # Radii of the doifferent layers
    r_w = Ri[s-1]*np.sqrt((1+e_w)**2 + (rad/Ri[s-1])**2 - 1)
    r_m = Ri[s-1]*np.sqrt((e_m+e_w+1)**2 + (rad/Ri[s-1])**2 - 1)
    r_s = 0.5*(r_w+r_m)
    
    # Expt data fit for AM force at latch state
    if r_s<=2*rsmax[s-1]:
        f_l = np.sin(np.pi*r_s/(2*rsmax[s-1]))**3
    else:
        f_l = 0
    
    #total active force
    fa = f_l+force_a
    # Radial Stress 
    D_Rm = (Rm[s-1]-r_m)/Rm[s-1]
    Sigma_rrm = 2*moo*(D_Rm + nu*D_Rm**2) + P_t
    Sigma_rrw = Sigma_rrm-fa*(Rm[s-1]-Rw[s-1])/(0.5*(r_m+r_w))
    
    # Hoop stress
    Sigma_ttm = Sigma_rrm + fa*(Rm[s-1]-Rw[s-1])/(r_m+r_w)
    Sigma_ttw = Sigma_rrw + fa*(Rm[s-1]-Rw[s-1])/(r_m+r_w)
        

    # rad = (np.sqrt(Ri_sq[s-1]*(1-(Ptm/P1[s-1]))**-N1[s-1]))*(Ptm<=0) + (np.sqrt(rmax_sq[s-1]-(rmax_sq[s-1]-Ri_sq[s-1])*(1-(Ptm/P2[s-1]))**-N2[s-1]))*(Ptm>=0)

    
    drdt = rho*(rad-r[7])
    
    # Velocity of contraction of SMC
    v_smc = (1/(r_w-r_m))*drdt*rad*np.log(r_w/r_m)
    V = -((gamma)/(2*np.pi*rmax[s-1]))*v_smc
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
    J20=((p1**(2))*phi0)+((2*p1*q1)*I0)+((q1**(2))*(2*phi0+(r0*I0)))
    J21=((p1**(2))*phi1)+((2*p1*q1)*I1)+((q1**(2))*(2*phi1+(r1*I1)))
    J30=(p1**(3)*phi0)+((3*p1**(2)*q1)*I0)+((3*p1*q1**(2))*((2*phi0)+(r0*I0)))+((q1**(3)*(3+r0**(2))*I0))
    J31=(p1**(3)*phi1)+((3*p1**(2)*q1)*I1)+((3*p1*q1**(2))*(2*phi1+(r1*I1)))+(q1**(3)*(3+(r1**(2)))*I1)
    
    #Functions defined for the RHS of the second PDE M2_lambda
    K0 = phi20
    K01 = phi21
    K10 = (p2*phi20) + (q2*I20)
    K11 = (p2*phi21) +( q2*I21)
    K20=((p2**(2))*phi20) + ((2*p2*q2)*I20) + ((q2**(2))*(2*phi20 + (r20*I20)))
    K21=((p2**(2))*phi21) + ((2*p2*q2)*I21) + ((q2**(2))*(2*phi21 + (r21*I21)))
    K30=(p2**(3)*phi20) + ((3*p2**(2)*q2)*I20) + ((3*p2*q2**(2))*((2*phi20) + (r20*I20))) + ((q2**(3)*(3+r20**(2))*I20))
    K31=(p2**(3)*phi21) + ((3*p2**(2)*q2)*I21) + ((3*p2*q2**(2))*(2*phi21 + (r21*I21))) + (q2**(3)*(3+(r21**(2)))*I21)
    
    # Components for the matrix F that will represent each moment,
    # M1_lambda and M2_lambda
  
    A0 = ((fp1*(1-r[6]))/1)-(fp1*(J11-J10)*r[0])
    A1 = ((fp1*(1-r[6]))/2)-(fp1*(J21-J20)*r[0])
    A2 = ((fp1*(1-r[6]))/3)-(fp1*(J31-J30)*r[0])
    
    B0=(gp2*J0)+gp1*(J11-J10)+(gp1+gp3)*(p1-J11);
    B1=(gp2*J10)+gp1*(J21-J20)+(gp1+gp3)*((p1**(2)+q1**(2))-J21);
    B2=(gp2*J20)+gp1*(J31-J30)+(gp1+gp3)*((p1**(3)+3*p1*q1**(2))-J31);
    
    C0=k1
    C1=k1*p2
    C2=k1*(p2**(2)+q2**(2))
    
    D0=k2
    D1=k2*p1
    D2=k2*(p1**(2)+q1**(2))
    
    E0=(g2*K0)+g1*(K11-K10)+(g1+g3)*(p2-K11);
    E1=(g2*K10)+g1*(K21-K20)+(g1+g3)*((p2**(2)+q2**(2))-K21);
    E2=(g2*K20)+g1*(K31-K30)+(g1+g3)*(p2**(3)+(3*p2*q2**(2))-K31);
    
    # Time derivative of force
    dfdt = -V*lam*(D0+C0)+A1+B0+E0
    
    
    F=np.array([A0-B0*r[0]+C0*r[3]-k2*r[0], A1-B1*r[0]+C1*r[3]-k2*r[1]-V*r[0], A2-B2*r[0]+C2*r[3]-k2*r[2]-2*V*r[1], 
                D0*r[0]-E0*r[3]-k1*r[3], D1*r[0]-E1*r[3]-k1*r[4]-V*r[3], D2*r[0]-E2*r[3]-k1*r[5]-2*V*r[4], -k1*r[6]+(1-r[6])*k2, rho*(rad-r[7])])
    return F,Ptm,r_w, r_m, Sigma_rrm, Sigma_rrw, Sigma_ttm, Sigma_ttw,dfdt

def USM_DM(trange,N,M,xrange,lam,s):
    T = np.linspace(trange[0], trange[1], N)
    X0 = np.linspace(xrange[0], xrange[1], M)
    dt=(trange[1]-trange[0])/(N)
    M10=0.005
    M11=0.001
    M12=0.01
    M20=0.005
    M21=0.001
    M22=0.01
    C=1
    r_Lum=0.4
 #   R = np.zeros(7)
 # Calculating the Ca2+ for the same trange
    
    R = np.array([M10, M11, M12, M20, M21, M22, C, r_Lum])
    Rnew = np.zeros([len(T),len(R)])
    Rnew[0,:] = R
    Ptm = np.zeros([len(T),1])
    r_w = np.zeros([len(T),1])
    r_m = np.zeros([len(T),1])
    Sigma_rrm = np.zeros([len(T),1])
    Sigma_rrw = np.zeros([len(T),1])
    Sigma_ttm = np.zeros([len(T),1])
    Sigma_ttw = np.zeros([len(T),1])
    dfdt = np.zeros([len(T),1])
    a = np.zeros([len(T),1])

    
    for i in range(1,len(T)):
        t = T[i-1]
        F1,Ptm[i],r_w[i], r_m[i], Sigma_rrm[i], Sigma_rrw[i], Sigma_ttm[i], Sigma_ttw[i], dfdt[i] = DM_funcs(t,Rnew[i-1,:],lam,s)
        F2,Ptm[i],r_w[i], r_m[i], Sigma_rrm[i], Sigma_rrw[i], Sigma_ttm[i], Sigma_ttw[i], dfdt[i] = DM_funcs(t+(dt/2),Rnew[i-1,:]+(dt/2)*F1,lam,s)
        F3,Ptm[i],r_w[i], r_m[i], Sigma_rrm[i], Sigma_rrw[i], Sigma_ttm[i], Sigma_ttw[i], dfdt[i] = DM_funcs(t+(dt/2),Rnew[i-1,:]+(dt/2)*F2,lam,s)
        F4,Ptm[i],r_w[i], r_m[i], Sigma_rrm[i], Sigma_rrw[i], Sigma_ttm[i], Sigma_ttw[i], dfdt[i] = DM_funcs(t+(dt),Rnew[i-1,:]+(dt)*F3,lam,s)
        Rnew[i,:] = Rnew[i-1,:]+ (dt/6)*(F1 + (2*F2) + (2*F3) + F4)
        # Rnew[i,:] = Rnew[i-1,:] + dt*F1
    return Rnew,Ptm,r_w, r_m, Sigma_rrm, Sigma_rrw, Sigma_ttm, Sigma_ttw,dfdt


        
        
        
trange = [0, 120]
xrange = [-30, 30]
N = 5000
M = 1000
# lam = [0,2,4,8,10]
s = 7
lam = 8

T = np.linspace(trange[0], trange[1], N)
# for i in range(len(lam)):
#     Rnew,Ptm,r_w, r_m, Sigma_rrm, Sigma_rrw, Sigma_ttm, Sigma_ttw, dfdt= USM_DM(trange, N, M, xrange,lam[i],s)
#     radius = Rnew[:,7]
#     plt.plot(T,radius)
    
# plt.show()
Rnew,Ptm,r_w, r_m, Sigma_rrm, Sigma_rrw, Sigma_ttm, Sigma_ttw, dfdt = USM_DM(trange, N, M, xrange,lam,s)
# force = (Rnew[:,1]+Rnew[:,4])
# stiffness = Rnew[:,0]+Rnew[:,3]
radius = Rnew[:,7]
# nmp = 1-Rnew[:,6]-Rnew[:,0]
# phosphorylation = nmp + Rnew[:,0]
# # plt.plot(T,a)
plt.plot(T,radius,T,r_w,T,r_m)
plt.show()