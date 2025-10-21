import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
import matplotlib.animation as animation
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
    
    def animate_funnel(rkrk, d_c, h_c_max, h_f_max, angle_f):
        # Create figure and axes
        fig, ax = plt.subplots(figsize=(6, 8))
        ax.set_xlim(-0.2, 0.2)
        ax.set_ylim(0, h_c_max + h_f_max + 0.1)
        ax.set_aspect('equal')
        ax.set_title("Olive Oil Draining Simulation", fontsize=12)
        ax.set_xlabel("Width (m)")
        ax.set_ylabel("Height (m)")
        
        # Cylinder geometry
        cyl_left = -d_c / 2
        cyl_right = d_c / 2
        gap = 0.02  
        cyl_bottom = h_f_max + gap
        cyl_top = cyl_bottom + h_c_max

        # Funnel geometry
        funnel_half_width_top = h_f_max / np.tan(angle_f)
        funnel_bottom_y = 0
        funnel_top_y = h_f_max
        funnel_points = [
            [0, funnel_bottom_y],
            [-funnel_half_width_top, funnel_top_y],
            [funnel_half_width_top, funnel_top_y]
        ]

        # Draw outlines
        ax.plot([cyl_left, cyl_left, cyl_right, cyl_right, cyl_left],
                [cyl_bottom, cyl_top, cyl_top, cyl_bottom, cyl_bottom],
                color='black')
        funnel_outline = Polygon(funnel_points, fill=False, color='black')
        ax.add_patch(funnel_outline)

        # Create patches for the oil in the cylinder and funnel
        oil_cylinder = Rectangle((cyl_left, cyl_bottom), d_c, 0, color='goldenrod', alpha=0.7)
        oil_funnel = Polygon([[0, 0]], color='goldenrod', alpha=0.7)
        ax.add_patch(oil_cylinder)
        ax.add_patch(oil_funnel)

        # Stream between cylinder and funnel
        stream, = ax.plot([0, 0], [0, 0], color='goldenrod', lw=2, alpha=0.6)

        # Add live text display (time, heights)
        info_text = ax.text(
            0.05, h_c_max + h_f_max + 0.05,
            "", fontsize=10, color='black', va='bottom'
        )

        def update(frame):
            # Extract current data
            h_c = rkrk.y1_values[frame]
            h_f = rkrk.y2_values[frame]
            t = rkrk.x_values[frame]

            # Update cylinder oil height
            oil_cylinder.set_height(h_c)
            oil_cylinder.set_y(cyl_bottom)

            # Update funnel shape
            top_width = h_f / np.tan(angle_f)
            funnel_liquid = np.array([
                [0, funnel_bottom_y],
                [-top_width, funnel_bottom_y + h_f],
                [top_width, funnel_bottom_y + h_f]
            ])
            oil_funnel.set_xy(funnel_liquid)

            # Update stream (only when draining)
            if h_c > 0.001:
                stream.set_data([0, 0], [funnel_top_y + gap, cyl_bottom])
            else:
                stream.set_data([], [])

            # Update live info text
            info_text.set_text(
                f"Time: {t:6.2f} s\n"
                f"h_c (cylinder): {h_c:6.3f} m\n"
                f"h_f (funnel): {h_f:6.3f} m"
            )

            return oil_cylinder, oil_funnel, stream, info_text

        ani = animation.FuncAnimation(
            fig, update, frames=range(0, len(rkrk.x_values), 1000),
            interval=0.5, blit=False, repeat=False
        )

        plt.tight_layout()
        plt.show()

        return ani



    rkrk.solve()
    rkrk.plot()
    ani = animate_funnel(rkrk, d_c, h_c_max, h_f_max, angle_f)
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