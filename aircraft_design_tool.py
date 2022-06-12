#!/usr/bin/env python3

#Wrapper for running like, everything
import sys
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
import tkinter as tk
import numpy as np
#from read_aircraft_txt import *
global V_LD_max
global g
global R
global M_entry
global T_entry
global S_entry
global pi
global cd0
global rho_sl
global viscosity
#Default Constants used a LOT
g=9.81;         #m/s^2
R=287;          #J/(Kg K)
pi=np.pi        #It's PI my GUY
Cd0=0.00984
rho_sl=1.225    #kg/m^3
viscosity=0.0000148     #m^2/s

def close_figures():
    plt.close('all')
    #plt.pause(.01)
    
def surface_area():
    global AR
    global S
    global mean_cord
    wingspan=float(wingspan_entry.get())
    root_cord=float(root_cord_entry.get())
    taper_ratio=float(taper_ratio_entry.get())
    sweep=float(sweep_entry.get())
    mean_cord=wingspan/6*((1+2*taper_ratio)/(1+taper_ratio))
    tip_cord=taper_ratio*root_cord
    S=wingspan*((tip_cord+root_cord)/2)
    AR=wingspan**2/S

    ans_S=tk.DoubleVar()
    ans_S.set(S)
    ans_AR=tk.DoubleVar()
    ans_AR.set(AR)
    tk.Label(text='Wing Area (m^2)').grid(row=2,column=5)
    tk.Label(textvariable=ans_S).grid(row=2,column=6)
    tk.Label(text='Aspect Ratio').grid(row=3,column=5)
    tk.Label(textvariable=ans_AR).grid(row=3,column=6)
    return S, AR, mean_cord

def altitude_characteristics():
    global rho_h
    global T_h
    h=float(h_entry.get())
    a_1=-6.5e-3;    #k/m
    T_sl=288.16;    #k
    rho_sl=1.225;   #g/m^3
    g=9.81;         #m/s^2
    R=287;          #J/(Kg K)
    if h<11000:
        T_h=T_sl+a_1*h;   #K
        rho_h=rho_sl*(T_h/T_sl)**(-1-g/(a_1*R))
    else:
        T_h=-56.5                        #K
        p=22.65*np.exp(1.73-0.000157*h)    #Kpa
        rho_h=p/(.2869*(T_h+273))           #Kg/m^3 MAYBE T_h+273.1
    print('The temperature is '+str( T_h)+'K')
    print('The density is '+str(rho_h)+'Kg/m^3')
    return rho_h, T_h
def altitude_characteristics_loop(h):
    a_1=-6.5e-3;    #k/m
    T_sl=288.16;    #k
    rho_sl=1.225;   #g/m^3
    g=9.81;         #m/s^2
    R=287;          #J/(Kg K)
    if h<11000:
        T_h_loop=T_sl+a_1*h;   #K
        rho_h_loop=rho_sl*(T_h_loop/T_sl)**(-1-g/(a_1*R))
    else:
        T_h=_loop-56.5                        #K
        p=22.65*np.exp(1.73-0.000157*h)    #Kpa
        rho_h_loop=p/(.2869*(T_h_loop+273))           #Kg/m^3 MAYBE T_h+273.1
    return rho_h_loop, T_h_loop
    
def lift_drag_calculator():
    global L
    global W
    global LD_max
    global CL_max
    global V_LD_max
    global TR_LD
    global v
    global CL
    global CL_needed
    global CD
    global LD
    global TR
    global k
    h=float(h_entry.get())
    M=float(M_entry.get())
    vc=float(vc_entry.get())
    v_min=float(v_min_entry.get())
    v_max=float(v_max_entry.get())
    CL=[]
    CD=[]
    LD=[]
    TR=[]
    k=(4/3)/(pi*6.16*.8)
    Cd0=0.00984
    #First, calculate 3D Lift and Drag
    L=M*g            #N
    W=L                 #N
    #CL=L/(.5*rho_h*vc**2*S)
    #CD=Cd0+k*CL**2
    #LD=CL/CD;
    v=np.linspace(v_min,v_max,num=100)
    #print('Coeff. 3D Lift is '+str(CL))
    #print('Coeff. 3D Drag is '+str(CD))
    #print('Lift/Drag Ratio is '+str(LD))
    #Lift to Drag Plots and Optimization
    ii=0
    for ii, jj in enumerate(v, start=0):
        CL_calc=L/(0.5*rho_h*v[ii]**2*S)
        CL=np.append(CL,CL_calc)
        CD=np.append(CD,(Cd0+k*CL[ii]**2))
        LD=np.append(LD,CL[ii]/CD[ii])
        TR=np.append(TR,(.5*rho_h*v[ii]**2*S*Cd0+(k*S/(.5*rho_h*v[ii]**2))*(W/S)**2)/1000)
        ii+=1
    LD_max=1/np.sqrt(4*Cd0*k)
    V_LD_max=np.sqrt((2/rho_h)*(W/S)*np.sqrt(k/Cd0)); #m/s
    #CL=np.delete(CL,0)
    #CD=np.delete(CD,0)
    #LD=np.delete(LD,0)
    TR_LD=.5*rho_h*V_LD_max**2*S*Cd0+(k*S)/(.5*rho_h*V_LD_max**2)*(W/S)**2
    q_infinity=0.5*rho_h*(V_LD_max**2)
    CL_needed=W/(q_infinity*S)
    Cl_needed=CL_needed*(1+(2/AR))
    reynolds_num_at_max_LD=(V_LD_max*mean_cord*rho_h)/viscosity
    print('L/D max is at '+str(LD_max))
    print('Airspeed for L/D max is '+str(V_LD_max))
    print('Thrust required at L/D max is '+str(TR_LD))
    plt.ion()
    plt.show()
    
    plt.figure(num=1)
    plt.title("Lift to Drag vs. Velocity")
    plt.xlabel("Velocity (m/s)")
    plt.ylabel("L/D Ratio")
    plt.plot(v,LD)
    
    plt.figure(num=2)
    plt.title("Thrust Needed vs. Velocity")
    plt.xlabel("Velocity (m/s)")
    plt.ylabel("Thrust (kN)")
    plt.plot(v,TR)
    
    plt.show()
    #Max LD Ratio
    ans_LD_max=tk.DoubleVar()
    ans_LD_max.set(LD_max)
    tk.Label(text='Max LD Ratio').grid(row=4,column=5)
    tk.Label(textvariable=ans_LD_max).grid(row=4,column=6)
    #V at LD Max
    ans_V_LD_max=tk.DoubleVar()
    ans_V_LD_max.set(V_LD_max)
    tk.Label(text='LD Max Speed (m/s)').grid(row=5,column=5)
    tk.Label(textvariable=ans_V_LD_max).grid(row=5,column=6)
    #CL
    ans_3d_cl_needed=tk.DoubleVar()
    ans_3d_cl_needed.set(CL_needed)
    tk.Label(text='CL at Max LD').grid(row=6,column=5)
    tk.Label(textvariable=ans_3d_cl_needed).grid(row=6,column=6)
    #Cl
    ans_2d_cl_needed=tk.DoubleVar()
    ans_2d_cl_needed.set(Cl_needed)
    tk.Label(text='Cl at Max LD').grid(row=7,column=5)
    tk.Label(textvariable=ans_2d_cl_needed).grid(row=7,column=6)
    print('Needed Cl is ' + str(Cl_needed))
    #Needed Thrust
    ans_TR_LD=tk.DoubleVar()
    ans_TR_LD.set(TR_LD)
    tk.Label(text='Thrust at Max LD (N)').grid(row=8,column=5)
    tk.Label(textvariable=ans_TR_LD).grid(row=8,column=6)
    #Reynolds
    ans_reynolds_num_at_max_LD=tk.DoubleVar()
    ans_reynolds_num_at_max_LD.set(reynolds_num_at_max_LD)
    tk.Label(text='Reynolds Number at Max LD').grid(row=9,column=5)
    tk.Label(textvariable=ans_reynolds_num_at_max_LD).grid(row=9,column=6)
    return L,W,LD_max,V_LD_max,TR_LD,v,CL,CD,LD,TR,k,Cd0

def climb_calc():
    global THsl
    T=float(T_entry.get())
    THsl=T*9.81;  #N Thrust at Sea Level
    M=float(M_entry.get())
    TW=T/W  #Thrust to Weight
    WS=W/S  #Weight to Area (Wing Loading)
    RC_max_rate=()  #m/s
    RC_v_max=[]     #m/s
    gamma=[]        #rad
    RC_at_v=[]      #m/s
    Z=1+np.sqrt(1+3/(LD_max**2*TW**2))
    RC_max_rate=np.sqrt((WS*Z)/(3*rho_h*Cd0))*TW**(3/2)*(1-Z/6-3/(2*TW**2*LD_max**2*2)) #m/s RC Max Rate
    RC_v_max=np.sqrt((TW*WS*Z)/(3*rho_h*Cd0))       #m/s Max V
    gamma=np.arcsin(RC_max_rate/RC_v_max)   #rad
    for ii,jj in enumerate(v,start=0):
        Q=.5*rho_h*v[ii]**2
        RC_at_v=np.append(RC_at_v,v[ii]*(TW-Q*Cd0/WS-WS*(k/Q))); #RC at velocity
    print('Max climb rate '+str(RC_max_rate)+' m/s')
    print('Airspeed at max climb rate '+str(RC_v_max)+' m/s')
    plt.figure(num=3)
    plt.title("Max Climb Rate vs. Velocity")
    plt.xlabel("Velocity (m/s)")
    plt.ylabel("Max Climb Rate (m/s)")
    plt.plot(v,RC_at_v)
    plt.show()
    return RC_max_rate,RC_v_max,gamma,RC_at_v

def climb_calc_loop(v,W,S,rho_h,T_h,LD_max,TW,k,Cd0):
    WS=W/S
    RC_max_rate=()  #m/s
    RC_v_max=[]     #m/s
    gamma=[]        #rad
    RC_at_v=[]      #m/s
    Z=1+np.sqrt(1+3/(LD_max**2*TW**2))
    RC_max_rate=np.sqrt((WS*Z)/(3*rho_h*Cd0))*TW**(3/2)*(1-Z/6-3/(2*TW**2*LD_max**2*2)) #m/s RC Max Rate
    RC_v_max=np.sqrt((TW*WS*Z)/(3*rho_h*Cd0))       #m/s Max V
    gamma=np.arcsin(RC_max_rate/RC_v_max)   #rad
    for ii,jj in enumerate(v,start=0):
        Q=.5*rho_h*v[ii]**2
        RC_at_v=np.append(RC_at_v,v[ii]*(TW-Q*Cd0/WS-WS*(k/Q))); #RC at velocity
    return RC_max_rate,RC_v_max,gamma,RC_at_v

def ceiling_calc_loop(v,W,S,rho_h,T_h,LD_max,TW,k,Cd0,h_floor,h_max,THsl):
        RC_max_rate_app=[]
        RC_v_max_app=[]
        gamma_app=[]
        rho_sl=1.225    #kg/m^3
        h_ceiling_range=np.linspace(h_floor,h_max,num=1000)
        for ii, jj in enumerate(h_ceiling_range,start=0):
                h=h_ceiling_range[ii]
                rho_h,T_h=altitude_characteristics(h)
                TW=THsl*(rho_h/rho_sl)/W
                RC_max_rate,RC_v_max,gamma,RC_at_v=h=climb_calc(v,W,S,rho_h,T_h,LD_max,TW,k,Cd0)
                RC_max_rate_app=np.append(RC_max_rate_app,RC_max_rate)
                RC_v_max_app=np.append(RC_v_max_app,RC_v_max)
                gamma_app=np.append(gamma_app,gamma)
        a=np.argmax(RC_max_rate_app<.508)
        max_service_ceiling=h[a]
        b=np.argmax(RC_max_rate_app<.01)
        abs_service_ceiling=h[b]
        return max_service_ceiling, abs_service_ceiling, RC_max_rate_app, RC_v_max_app, h_ceiling_range 


def ceiling_calc():
    h_max=float(h_max_entry.get())
    h_floor=float(h_floor_entry.get())
    RC_max_rate_app=[]
    RC_v_max_app=[]
    gamma_app=[]
    rho_sl=1.225    #kg/m^3
    W=L 
    h_ceiling_range=np.linspace(h_floor,h_max,num=1000)
    for ii, jj in enumerate(h_ceiling_range,start=0):
        h=h_ceiling_range[ii]
        rho_h_loop,T_h_loop=altitude_characteristics_loop(h)
        TW=THsl*(rho_h_loop/rho_sl)/W
        RC_max_rate,RC_v_max,gamma,RC_at_v=h=climb_calc_loop(v,W,S,rho_h,T_h,LD_max,TW,k,Cd0)
        RC_max_rate_app=np.append(RC_max_rate_app,RC_max_rate)
        RC_v_max_app=np.append(RC_v_max_app,RC_v_max)
        gamma_app=np.append(gamma_app,gamma)
    plt.figure(num=4)
    plt.title('RCmax Versus Altitude')
    plt.xlabel('Altitude (m)')
    plt.ylabel('RC Max (m/s)')
    plt.plot(h_ceiling_range,RC_max_rate_app)
    plt.show()
    a=np.argmax(RC_max_rate_app<.508)
    max_service_ceiling=h[a]
    b=np.argmax(RC_max_rate_app<.01)
    abs_service_ceiling=h[b]
    print('Max Service Ceiling is '+str(max_service_ceiling)+' m')
    print('Abs Service Ceiling is '+str(abs_service_ceiling)+' m')
    return max_service_ceiling, abs_service_ceiling, RC_max_rate_app, RC_v_max_app, h_ceiling_range 

def turn_radius_calc(W,h,THsl,v,Cd0,S,k,g_factor):
        M=float(M_entry.get())
        rho_sl=1.225    #Kg/m^3
        g=9.81          #m/s^2
        rho_h,T_h=altitude_characteristics_loop(h)
        T=THsl*(rho_h/rho_sl)
        Q=.5*rho_h*v**2         #kPa
        CL_max=(M*g)/(Q*S)      #CL needed. 
        n_max_s=g_factor        #n_max_s
        n_max_t=np.sqrt(np.abs(Q/(k*(W/S))*((T/W)-((Q*Cd0)/(W/S)))))    # REMOVE ABSn_max_t
        wing_load=W/S
        n_max_alpha=(Q*CL_max)/(W/S)                            #n_max_alpha
        R_min=(4*k*(W/S))/(g*rho_h*(T/W)*np.sqrt(1-(4*k*Cd0)/(T/W)**2)) #m
        V_R_min=np.sqrt((4*k*(W/S))/(rho_h*(T/W)))                      #m/s
        return n_max_s,n_max_t,n_max_alpha,R_min,V_R_min
    
def varying_v_turn_radius_calc():
        M=float(M_entry.get())
        T=float(T_entry.get())
        h=float(h_entry.get())
        v_min=float(v_min_entry.get())
        v_max=float(v_max_entry.get())
        phi_factor=float(phi_factor_entry.get())
        phi=np.deg2rad(phi_factor) #Rad
        g_factor=1/np.cos(phi) 
        THsl=T*9.81;  #N Thrust at Sea Level
        V=np.linspace(v_min,v_max,num=100)
        g=9.81          #m/s^2
        R_t=[]
        R_a=[]
        for ii, jj in enumerate(V,start=0):
                v=V[ii] #m/s
                n_max_s,n_max_t,n_max_alpha,R_min,V_R_min=turn_radius_calc(W,h,THsl,v,Cd0,S,k,g_factor)
                R_t=np.append(R_t,v**2/(g*np.sqrt(n_max_t**2-1)))
                R_a=np.append(R_a,v**2/(g*np.sqrt(n_max_alpha**2-1)))
                print(str(n_max_alpha))
        print('I NEED TO FIND OUT WHAT CL MAX IS. MAX FOR AIRFOIL OR MAX FOR VELOCITY?')
        print(str(R_t))
        print(str(R_a))
        print(str(V))
        plt.figure(num=5)
        plt.title('Turn Radius with Varying Velocity at '+str(h)+' (m)')
        plt.xlabel('Velocity (m/s)')
        plt.ylabel('Turn Radius (m)')
        plt.legend(['Thrust Turn Radius','Alpha Turn Radius'])
        plt.plot(V,R_t)
        plt.plot(V,R_a)
        plt.show()
        return R_t,R_a,V
def range_equation_elect_calc():
        M=float(M_entry.get())
        motor_draw=float(motor_draw_entry.get())    #A
        bat_cap=float(Ah_entry.get())               #Ah
        time_of_flight=(bat_cap*3600)/motor_draw*.8
        print('Time of Flight: ' + str(time_of_flight)+ ' Seconds')
        ans_time_of_flight=tk.DoubleVar()
        ans_time_of_flight.set(time_of_flight)
        tk.Label(text='Time of Flight (s)').grid(row=10,column=5)
        tk.Label(textvariable=ans_time_of_flight).grid(row=10,column=6)
        ans_time_of_flight=tk.DoubleVar()   
        ans_time_of_flight.set(time_of_flight)
        tk.Label(text='Time of Flight (s)').grid(row=10,column=5)
        tk.Label(textvariable=ans_time_of_flight).grid(row=10,column=6)
        return time_of_flight

def pitch_stability_calc():
        # Needs CL_aw,CLdet,St,S,CL_at,K_ea,hcm,hac,lt,C_mean,CM_ac_w,CLw_0,i_trim,e0,AOA,h,v,W
        wingspan=float(wingspan_entry.get()) #m
        M=float(M_entry.get())  #Kg
        W=M*g                   #N
        h=float(h_entry.get())  #m
        Cvt=float(vert_tail_coeff_entry.get())   #Vertical Tail Volume Coefficent
        Cht=float(horz_tail_coeff_entry.get())   #Horizontal Tail Volume Coefficent
        v=V_LD_max              #m/s
        C_mean=mean_cord        #m, global
        lth=float(length_to_horizontal_tail_entry.get())  #m horizontal
        ltv=float(length_to_vertical_tail_entry.get())  #m Vertical Tail
        CL_at=3.1615 #1/rad (Constant)
        CL_aw=4.698   #1/rad
        CLw_0=float(coef_lift_naught_entry.get())         #Airfoil characteristic
        CM_ac_w =float(coef_mom_wing_entry.get())    #Airfoil characteristic
        hac =float(hac_entry.get())            #m User defined
        hcm=float(hcm_entry.get())             #m User defined
        i_trim=np.deg2rad(float(itrim_entry.get())) # rad tail trim angle
        e0=0        #Constant
        K_ea=.1     #Constant
        AOA=np.deg2rad(float(AOA_entry.get()))
        CLdet=float(CLdet_entry.get())
        rho_h,T_h=altitude_characteristics_loop(h)
        St=Cht*mean_cord*S/lth          #Horizontal Tail Area (m^2)
        Sht=St                          #Horizontal Tail Area (m^2)
        Svt=Cvt*wingspan*S/ltv          #Vertical Tail Area (m^2)
        CL_a=CL_aw+Sht/S*CL_at*(1-K_ea)  #Lift as vs. AoA
        CL_0=CL_at*(-e0-i_trim)+CLw_0   #Initial Lift
        CM_0=CM_ac_w+CL_0*(hcm-hac)-Cht*CL_at*(i_trim-e0)        #Moment 0
        CM_a=CL_a*(hcm-hac)-Cht*CL_at*(1-K_ea) #Cm f(AoA)
        CL_trim=CL_0+CL_a*(-(CM_0/CM_a))
        CM=CM_0+CM_a*AOA                #Cm Total
        hnp=hac+(CL_at/CL_a)*Cht*(1-K_ea)        #Neutral Point
        static_margin=hnp-hcm           #Static Margin
        CLde=Sht/S*CLdet                 #CLde
        CMde=CLde*(hcm-hac)-CLdet*Cht    #CMde
        # Tail and Elevator Trim
        A=np.array([[CL_a,CLde],[CM_a,CMde]]) #A Matrix
        n=np.array([[CL_needed-CL_0],[-CM_0]])
        T=np.matmul(np.linalg.inv(A),n)
        tail_trim=T[0]
        elevator_trim=T[1]
        tail_trim=np.rad2deg(tail_trim)
        elevator_trim=np.rad2deg(elevator_trim)
        v_trim=np.sqrt(W/(.5*rho_h*S*CL_trim))
        print('The static margin is '+ str(static_margin))
        print('Horizontal Tail area is '+str(Sht)+' m^2')
        print('Vertical Tail Area is '+str(Svt)+' m^2')
        print('Trim Angle of Tail is ' +str(tail_trim) +' Deg')
        print('Trim Angle of Elevator is ' +str(elevator_trim) +' Deg')
        print('Velocity of Trim is ' + str(v_trim) +' m/s')
        ans_Sht=tk.DoubleVar()
        ans_Sht.set(Sht)
        tk.Label(text='Horizontal Tail Area (m^2)').grid(row=11,column=5)
        tk.Label(textvariable=ans_Sht).grid(row=11,column=6)
        ans_Svt=tk.DoubleVar()
        ans_Svt.set(Svt)
        tk.Label(text='Vertical Tail Area (m^2)').grid(row=12,column=5)
        tk.Label(textvariable=ans_Svt).grid(row=12,column=6)
        ans_static_stability=tk.DoubleVar()
        ans_static_stability.set(static_margin)
        tk.Label(text='Static Margin').grid(row=13,column=5)
        tk.Label(textvariable=ans_static_stability).grid(row=13,column=6)


#Defining Window
mwr=tk.Tk()
mwr.title("Aircraft Design Tool")
mwr.grid()
#Defining Input Labels and inputs
row_num=1
row_ans_num=1
row_but_num=2
tk.Label(text='Aircraft Characteristics').grid(row=row_num,column=2)
row_num=row_num+1
#Cruise Altitude (m)
default_h=tk.StringVar()
default_h.set('10')
h_entry=tk.Entry(mwr, textvariable=default_h)
h_entry.grid(row=row_num,column=3)
tk.Label(text='Height (h) in m').grid(row=row_num,column=2)
row_num=row_num+1
#Mass (Kg)
default_M=tk.StringVar()
default_M.set('2')
M_entry=tk.Entry(mwr, textvariable=default_M)
M_entry.grid(row=row_num,column=3)
tk.Label(text='Mass (Kg)').grid(row=row_num,column=2)
row_num=row_num+1
#Wingspan (m)
tk.Label(text='Wingspan (m)').grid(row=row_num,column=2)
default_wingspan=tk.StringVar()
default_wingspan.set('1.5')
wingspan_entry=tk.Entry(mwr,textvariable=default_wingspan)
wingspan_entry.grid(row=row_num,column=3)
row_num=row_num+1
#Root Cord (m)
tk.Label(text='Root Cord (m').grid(row=row_num,column=2)
default_root_cord=tk.StringVar()
default_root_cord.set('.15')
root_cord_entry=tk.Entry(mwr,textvariable=default_root_cord)
root_cord_entry.grid(row=row_num,column=3)
row_num=row_num+1
#Taper Ratio (N/A)
tk.Label(text='Taper Ratio (1 for None)').grid(row=row_num,column=2)
default_taper_ratio=tk.StringVar()
default_taper_ratio.set('1')
taper_ratio_entry=tk.Entry(mwr,textvariable=default_taper_ratio)
taper_ratio_entry.grid(row=row_num,column=3)
row_num=row_num+1
#Sweep (Deg)
tk.Label(text='Sweep (Deg)').grid(row=row_num,column=2)
default_sweep=tk.StringVar()
default_sweep.set('0')
sweep_entry=tk.Entry(mwr,textvariable=default_sweep)
sweep_entry.grid(row=row_num,column=3)
row_num=row_num+1
#Cruise Velocity (m/s)
tk.Label(text='Desired Cruise Velocity (m/s)').grid(row=row_num,column=2)
default_vc=tk.StringVar()
default_vc.set('22')
vc_entry=tk.Entry(mwr, textvariable=default_vc)
vc_entry.grid(row=row_num,column=3)
row_num=row_num+1
#Min Velocity (m/s)
tk.Label(text='Minimum Velocity (m/s)').grid(row=row_num,column=2)
default_v_min=tk.StringVar()
default_v_min.set('5')
v_min_entry=tk.Entry(mwr, textvariable=default_v_min)
v_min_entry.grid(row=row_num,column=3)
row_num=row_num+1
#Max Velocity (m/s)
tk.Label(text='Maximum Velocity (m/s)').grid(row=row_num,column=2)
default_v_max=tk.StringVar()
default_v_max.set('30')
v_max_entry=tk.Entry(mwr, textvariable=default_v_max)
v_max_entry.grid(row=row_num,column=3)
row_num=row_num+1
#Thrust (N)
tk.Label(text='Thrust (N)').grid(row=row_num,column=2)
default_T=tk.StringVar()
default_T.set('1.5')
T_entry=tk.Entry(mwr, textvariable=default_T)
T_entry.grid(row=row_num,column=3)
row_num=row_num+1
#h floor (m)
tk.Label(text='Min Floor (m)').grid(row=row_num,column=2)
default_h_floor=tk.StringVar()
default_h_floor.set('15')
h_floor_entry=tk.Entry(mwr, textvariable=default_h_floor)
h_floor_entry.grid(row=row_num,column=3)
row_num=row_num+1
#h max (m)
tk.Label(text='Max Floor (m)').grid(row=row_num,column=2)
default_h_max=tk.StringVar()
default_h_max.set('1000')
h_max_entry=tk.Entry(mwr, textvariable=default_h_max)
h_max_entry.grid(row=row_num,column=3)
row_num=row_num+1
#h max (m)
tk.Label(text='phi_factor (Deg)').grid(row=row_num,column=2)
default_phi_factor=tk.StringVar()
default_phi_factor.set('45')
phi_factor_entry=tk.Entry(mwr, textvariable=default_phi_factor)
phi_factor_entry.grid(row=row_num,column=3)
row_num=row_num+1
#Battery Amp Hours(Ah)
tk.Label(text='Battery Amp Hours (Ah)').grid(row=row_num,column=2)
default_Ah=tk.StringVar()
default_Ah.set('5')
Ah_entry=tk.Entry(mwr, textvariable=default_Ah)
Ah_entry.grid(row=row_num,column=3)
row_num=row_num+1
#Motor Draw Amperes (A)
tk.Label(text='Motor Draw (A)').grid(row=row_num,column=2)
default_motor_draw=tk.StringVar()
default_motor_draw.set('12')
motor_draw_entry=tk.Entry(mwr, textvariable=default_motor_draw)
motor_draw_entry.grid(row=row_num,column=3)
row_num=row_num+1
#Length to  Horizontal Tail (m)
tk.Label(text='Distance to  Hor. Tail (m)').grid(row=row_num,column=2)
default_length_to_tail=tk.StringVar()
default_length_to_tail.set('.66')
length_to_horizontal_tail_entry=tk.Entry(mwr, textvariable=default_length_to_tail)
length_to_horizontal_tail_entry.grid(row=row_num,column=3)
row_num=row_num+1
#Length to Vertical Tail (m)
tk.Label(text='Distance to  Vert. Tail (m)').grid(row=row_num,column=2)
default_length_to_v_tail=tk.StringVar()
default_length_to_v_tail.set('.66')
length_to_vertical_tail_entry=tk.Entry(mwr, textvariable=default_length_to_v_tail)
length_to_vertical_tail_entry.grid(row=row_num,column=3)
row_num=row_num+1
#Horizontal Tail Volume Coefficent
tk.Label(text='Horz. Tail Volume Coeff.').grid(row=row_num,column=2)
default_horz_tail_coeff=tk.StringVar()
default_horz_tail_coeff.set('.5')
horz_tail_coeff_entry=tk.Entry(mwr, textvariable=default_horz_tail_coeff)
horz_tail_coeff_entry.grid(row=row_num,column=3)
row_num=row_num+1
#Vertical Tail Volume Coefficent
tk.Label(text='Vert. Tail Volume Coeff.').grid(row=row_num,column=2)
default_vert_tail_coeff=tk.StringVar()
default_vert_tail_coeff.set('.04')
vert_tail_coeff_entry=tk.Entry(mwr, textvariable=default_vert_tail_coeff)
vert_tail_coeff_entry.grid(row=row_num,column=3)
row_num=row_num+1
#Coefficent of Lift Slope: Tail
tk.Label(text='Coef. Lift Tail (1/rad)').grid(row=row_num,column=2)
default_coef_lift_slope_tail=tk.StringVar()
default_coef_lift_slope_tail.set('3.615')
coef_lift_slope_tail_entry=tk.Entry(mwr, textvariable=default_coef_lift_slope_tail)
coef_lift_slope_tail_entry.grid(row=row_num,column=3)
row_num=row_num+1
#Coefficent Lift Slope: Wing
tk.Label(text='Coef. Lift Wing (1/rad)').grid(row=row_num,column=2)
default_coef_lift_slope_wing=tk.StringVar()
default_coef_lift_slope_wing.set('4.698')
coef_lift_slope_wing_entry=tk.Entry(mwr, textvariable=default_coef_lift_slope_wing)
coef_lift_slope_wing_entry.grid(row=row_num,column=3)
row_num=row_num+1
#Coefficent Wing Naught
tk.Label(text='Coef. Lift Wing Naught').grid(row=row_num,column=2)
default_coef_lift_naught=tk.StringVar()
default_coef_lift_naught.set('0.7')
coef_lift_naught_entry=tk.Entry(mwr, textvariable=default_coef_lift_naught)
coef_lift_naught_entry.grid(row=row_num,column=3)
row_num=row_num+1
#Coefficent of Moment Wing
tk.Label(text='Coef. Mom. Wing').grid(row=row_num,column=2)
default_coef_mom_wing=tk.StringVar()
default_coef_mom_wing.set('-.1')
coef_mom_wing_entry=tk.Entry(mwr, textvariable=default_coef_mom_wing)
coef_mom_wing_entry.grid(row=row_num,column=3)
row_num=row_num+1
#hac center of effort
tk.Label(text='hac').grid(row=row_num,column=2)
default_hac=tk.StringVar()
default_hac.set('.25')
hac_entry=tk.Entry(mwr, textvariable=default_hac)
hac_entry.grid(row=row_num,column=3)
row_num=row_num+1
#hcm center of mass relative to wing front
tk.Label(text='hcm').grid(row=row_num,column=2)
default_hcm=tk.StringVar()
default_hcm.set('.25')
hcm_entry=tk.Entry(mwr, textvariable=default_hcm)
hcm_entry.grid(row=row_num,column=3)
row_num=row_num+1
#i_trim
tk.Label(text='Trim Angle (Deg)').grid(row=row_num,column=2)
default_itrim=tk.StringVar()
default_itrim.set('-.2')
itrim_entry=tk.Entry(mwr, textvariable=default_itrim)
itrim_entry.grid(row=row_num,column=3)
row_num=row_num+1
#AOA for Trim
tk.Label(text='Angle of Attack (Deg)').grid(row=row_num,column=2)
default_AOA=tk.StringVar()
default_AOA.set('0')
AOA_entry=tk.Entry(mwr, textvariable=default_AOA)
AOA_entry.grid(row=row_num,column=3)
row_num=row_num+1
#Long. Control Derivative
tk.Label(text='CLdet').grid(row=row_num,column=2)
default_CLdet=tk.StringVar()
default_CLdet.set('.5')
CLdet_entry=tk.Entry(mwr, textvariable=default_CLdet)
CLdet_entry.grid(row=row_num,column=3)
row_num=row_num+1

# Run Default Calculations so System doesn't get mad
##default_calc()
#Defining Buttons
#Close Plots
lift_drag_calculator_button = tk.Button(mwr,text='Close Plots',command=close_figures)
lift_drag_calculator_button.grid(row=20,column=6)

#Update Area and Aspect Ratio surface_area
button_altitude_characteristics = tk.Button(mwr,text='Update Wing Area and AR',command=surface_area)
button_altitude_characteristics.grid(row=row_but_num,column=4)
row_but_num=row_but_num+1
#Atmosphere Characteristics
button_altitude_characteristics = tk.Button(mwr,text='Get Altitude Characteristics',command=altitude_characteristics)
button_altitude_characteristics.grid(row=row_but_num,column=4)
row_but_num=row_but_num+1
#List and Drag Calc
lift_drag_calculator_button = tk.Button(mwr,text='Lift Drag Calculator',command=lift_drag_calculator)
lift_drag_calculator_button.grid(row=row_but_num,column=4)
row_but_num=row_but_num+1
#List and Drag Calc
climb_calc_button = tk.Button(mwr,text='Climb Calculator',command=climb_calc)
climb_calc_button.grid(row=row_but_num,column=4)
row_but_num=row_but_num+1
#ceiling_calc
ceiling_calc_button = tk.Button(mwr,text='Ceiling Calculator',command=ceiling_calc)
ceiling_calc_button.grid(row=row_but_num,column=4)
row_but_num=row_but_num+1
#Varying Turn Calc
turn_calc_button = tk.Button(mwr,text='Turn Calculator',command=varying_v_turn_radius_calc)
turn_calc_button.grid(row=row_but_num,column=4)
row_but_num=row_but_num+1
#Range Calc
range_calc_button = tk.Button(mwr,text='Range Calculator',command=range_equation_elect_calc)
range_calc_button.grid(row=row_but_num,column=4)
row_but_num=row_but_num+1
#Stability Calc
static_stability_calc_button = tk.Button(mwr,text='Static Stability Calculator',command=pitch_stability_calc)
static_stability_calc_button.grid(row=row_but_num,column=4)
row_but_num=row_but_num+1
#Control Derivatives
Control_calc_button = tk.Button(mwr,text='Control Derivative Calculator',command=pitch_stability_calc)
Control_stability_calc_button.grid(row=row_but_num,column=4)
row_but_num=row_but_num+1



mwr.mainloop()
