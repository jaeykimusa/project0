import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.patches import Rectangle, Polygon
import matplotlib.patches as mpatches

# Import the simulation class
import sys

from rk4 import RK4

# class OliveOil:
#     def __init__(
#             this,
#     ):
#         this.

class OliveOilSystem:
    def __init__(
                this, 
                functions: callable, 
                d_c: float,
                d_ce: float,
                h_c0: float,
                h_c_max: float,
                d_fe: float,
                h_f0: float,
                h_f_max: float,
                angle_f: float,
                C_d: float,
                g: float,
                rho_oil: float,
                dt: float
            ):
        # problem system variables
        # this.f1 = functions[0]
        # this.f2 = functions[1]
        this.d_c = d_c
        this.d_ce = d_ce
        this.h_c = h_c0
        this.h_c_max = h_c_max
        this.d_f = 0
        this.d_fe = d_fe
        this.h_f = h_f0
        this.h_f_max = h_f_max
        this.angle_f = angle_f
        this.C_d = C_d
        this.g = g
        this.rho_oil = rho_oil
        this.dt = dt

        # fluid system variables
        this.v_c_out = this.C_d*np.sqrt(2*this.g*this.h_c)
        this.A_c = 1/4*np.pi*d_c**2
        this.A_ce = 1/4*np.pi*d_ce**2
        this.Q_c_out = this.A_ce * this.v_c_out

        this.v_f_out = this.C_d*np.sqrt(2*this.g*this.h_f)
        this.A_f = 1/4*np.pi*this.d_f**2
        this.A_fe = 1/4*np.pi*this.d_fe**2
        this.Q_f_in = this.Q_c_out
        this.Q_f_out = this.A_fe * this.v_f_out

        # define RK4 problem
        this.problem = RK4(functions, this.h_c, this.h_f, 0, 3, this.dt)


    def step(this):
        this.problem.solve()
        # this.problem.step()
        # this.problem.printStep()


def main():

    # cylindrical container variables
    d_c = 0.25 # diameter (m)
    d_ce = 0.01 # exit hole diameter (m)
    h_c0 = 0.25 # initial height (m)
    h_c_max = 0.30 # height (m) ***custom variable***

    # cone funnel variables
    d_fe = 0.01 # exit hole diameter (m)
    h_f0 = 0.001 # initial height (m)
    h_f_max = 0.15 # hieght (m)
    d_fe2 = 0.02 # exit hole diameter for scenario 2 (m)
    angle_f = 45 * np.pi/180 # angle (rad)

    # other constants
    C_d = 0.65 # discharge coefficient
    g = 9.81 # gravity constant (g)
    rho_water = 997 # water density (kg/m^3)
    rho_oil = 0.90 * rho_water # olive oil density (kg/m^3)
    dt = 0.01 # time step (sec)

    # system variables
    A_c = 1/4*np.pi*d_c**2
    A_ce = (1/4)*np.pi*(d_ce**2)
    A_fe = (1/4)*np.pi*(d_fe**2)

    def computeFunnelArea(h_f, angle_f):
        return np.pi * (h_f/np.tan(angle_f))**2

    # this is the first equation in the system
    def func1(t, h_c, h_f):
        return -1*((d_ce**2)/(d_c**2))*C_d*np.sqrt(2*g*h_c)
    
    # this is the second equation in the system
    def func2(t, h_c, h_f):
        # A_f = computeFunnelArea(h_f, angle_f)
        # return (1/A_f) * (A_ce*C_d*np.sqrt(2*g*h_c) - A_fe*C_d*np.sqrt(2*g*h_f))
        # return (A_ce*C_d*np.sqrt(2*g*h_c) - A_fe*C_d*np.sqrt(2*g*h_f)) / (np.pi * ((h_f_max-h_f)*np.tan(angle_f)**2))
        # return ((np.tan(angle_f)**2)/(np.pi*(h_f**2)))*(A_ce*C_d*np.sqrt(2*g*h_c) - A_fe*C_d*np.sqrt(2*g*h_f))
        # return (C_d/(4*h_f**2))*((d_ce**2)*np.sqrt(2*g*h_c)-(d_fe**2)*np.sqrt(2*g*h_f))
        # return (1/(4*(h_f**2)))*((0.01**2)*0.65*np.sqrt(2*9.81*h_c)-(0.01**2)*np.sqrt(2*9.81*h_f))

        return (np.sqrt(2*9.81)/(4*(h_f**2)))*(((0.01**2)*0.65*np.sqrt(h_c))-(0.01**2)*np.sqrt(h_f))

    
    rkrk = RK4(
                [func1, func2],
                y1_0=h_c0,
                y2_0=h_f0,
                x_lb=0,
                x_ub=200,
                h=0.001
            )
    rkrk.solve()
    rkrk.plot()
        # problem = RK4(
        #             [func1, func2], 
        #             y1_0=-np.pi, 
        #             y2_0=-np.pi, 
        #             x_lb=0, 
        #             x_ub=np.pi, 
        #             h=np.pi/4
        #         )

    # test = OliveOilSystem(
    #                 [func1, func2], 
    #                 d_c=d_c,
    #                 d_ce=d_ce,
    #                 h_c0=h_c0,
    #                 h_c_max=h_c_max,
    #                 d_fe=d_fe,
    #                 h_f0=h_f0,
    #                 h_f_max=h_f_max,
    #                 angle_f=angle_f,
    #                 C_d=C_d,
    #                 g=g,
    #                 rho_oil=rho_oil,
    #                 dt=dt
    #             )
    # test.step()

if __name__ == "__main__":
    main()


