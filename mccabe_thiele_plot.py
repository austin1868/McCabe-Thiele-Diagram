import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

def mccabe_thiele_plot():
    # Collect equilibrium data
    num_points = int(input("Enter the number of equilibrium data points: "))
    x_values = []
    y_values = []
    for i in range(num_points):
        x = float(input(f"Enter x value for point {i + 1}: "))
        y = float(input(f"Enter y value for point {i + 1}: "))
        x_values.append(x)
        y_values.append(y)
    
    # Collect values for Rectifying Section
    R = float(input("Enter the reflux ratio (R): "))
    x_D = float(input("Enter the distillate composition (x_D): "))
    
    # Collect values for Stripping Section
    L_prime = float(input("Enter the L' value (ratio of liquid flow rate in stripping section to vapor flow rate from reboiler): "))
    z_F = float(input("Enter the feed composition (z_F): "))
    x_B = float(input("Enter the bottoms composition (x_B): "))
    
    # Calculate operating lines
    x = np.linspace(0, 1, 100)
    y_rectifying = R * x + ((R - 1) * x_D) / (R + 1)
    y_stripping = L_prime * x + (z_F - L_prime * x_B)
    
    # Plotting
    plt.figure(figsize=(10, 8))
    plt.plot(x_values, y_values, label="Equilibrium Curve", color="blue")
    plt.plot(x, y_rectifying, label="Rectifying Operating Line", color="red")
    plt.plot(x, y_stripping, label="Stripping Operating Line", color="green")
    plt.plot([z_F, z_F], [0, z_F], label="Feed Line", color="purple", linestyle="--")
    
    # Labels and title
    plt.xlabel("x (Liquid phase mole fraction)")
    plt.ylabel("y (Vapor phase mole fraction)")
    plt.title("McCabe-Thiele Diagram")
    plt.legend()
    plt.grid(True)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.show()

# For now, we won't execute the function to avoid interactive inputs in this environment.
# But, you can run the function `mccabe_thiele_plot()` in a local environment to interactively provide inputs and see the plot.

mccabe_thiele_plot()