# rk4.py

# Solves second-order ODEs using fourth-order Runge-Kutta method.
# Author: Jaey Kim
# Date: 10/01/2025


import numpy as np
import matplotlib.pyplot as plt  # <-- import matplotlib for plotting

class RK4:

    # ... your existing __init__ method, but add storage lists:
    def __init__(
            this, 
            functions: callable, 
            y1_0: float, 
            y2_0: float,
            x_lb: float, 
            x_ub: float, 
            h: float
        ):
        # existing initialization...
        this.f1 = functions[0]
        this.f2 = functions[1]
        this.y1 = y1_0
        this.y2 = y2_0
        this.x = x_lb
        this.x_lb = x_lb
        this.x_ub = x_ub
        this.h = h
        this.total_steps = int((x_ub - x_lb) / h)
        this.current_step = 0
        
        # Initialize lists to store values for plotting
        this.x_values = [x_lb]
        this.y1_values = [y1_0]
        this.y2_values = [y2_0]

    def step(this):
        # existing step code...
        h = this.h
        x = this.x
        y1 = this.y1
        y2 = this.y2

        f1 = this.f1
        f2 = this.f2

        k11 = f1(x, y1, y2)
        k12 = f2(x, y1, y2)
        k21 = f1(x + h/2, y1 + (1/2)*k11*h, y2 + (1/2)*k12*h)
        k22 = f2(x + h/2, y1 + (1/2)*k11*h, y2 + (1/2)*k12*h)
        k31 = f1(x + h/2, y1 + (1/2)*k21*h, y2 + (1/2)*k22*h)
        k32 = f2(x + h/2, y1 + (1/2)*k21*h, y2 + (1/2)*k22*h)
        k41 = f1(x + h, y1 + k31*h, y2 + k32*h)
        k42 = f2(x + h, y1 + k31*h, y2 + k32*h)

        this.x += h
        this.current_step += 1
        this.y1 += (1/6)*(k11 + 2*k21 + 2*k31 + k41)*h
        this.y2 += (1/6)*(k12 + 2*k22 + 2*k32 + k42)*h

        # Store values after step
        this.x_values.append(this.x)
        this.y1_values.append(this.y1)
        this.y2_values.append(this.y2)


    def getValues(this):
        return this.y1, this.y2

    def printStep(this):
        '''
        Prints the result of a single RK4 step.

        Parameters:
            None
        Returns:
            None
        '''
        print("-"*15, "Step", this.current_step, "-"*15)
        print(f" x_{this.current_step} = {this.x:.9f}")
        print(f" y1(x_{this.current_step}) = {this.y1:.9f}")
        print(f" y2(x_{this.current_step}) = {this.y2:.9f}\n")
    

    def solve(this):
        '''
        Solves to compute and print the problem results.

        Parameters:
            None
        Returns:
            None
        '''
        this.printStep()
        while(this.current_step < this.total_steps):
            this.step()
            this.printStep()
    
    def plot(this):
        """
        Plots x vs y1 and x vs y2 on two separate subplots.
        """
        plt.figure(figsize=(10, 5))

        # Plot x vs y1
        plt.subplot(1, 2, 1)
        plt.plot(this.x_values, this.y1_values, label='y1(x)', color='b')
        plt.xlabel('Time (s)')
        plt.ylabel('Olive oil height in container (m)')
        plt.title('Plot of y1 vs x')
        plt.grid(True)
        plt.legend()

        # Plot x vs y2
        plt.subplot(1, 2, 2)
        plt.plot(this.x_values, this.y2_values, label='y2(x)', color='r')
        plt.xlabel('Time (s)')
        plt.ylabel('Plive oil height in funnel (m)')
        plt.title('Plot of y2 vs x')
        plt.grid(True)
        plt.legend()

        plt.tight_layout()
        # plt.show()
        plt.savefig("funnel0.02.png")

    def getMaxFunnelHeight(this):
        return max(this.y2_values)