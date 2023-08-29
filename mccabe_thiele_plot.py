import numpy as np
import matplotlib.pyplot as plt

def mccabe_thiele_plot(alpha, zf, q, xd, xb, RR, eta):
    #inputs

    # find xp
    if q == 1:
        xp = zf
    elif q == 0:
        xp = zf / (alpha * zf * (alpha - 1))
    else:
        xp = (((alpha-1)*(zf+q)-alpha)+np.sqrt(((alpha-1)*(zf+q)-alpha)**2+4*zf*(alpha-1)*q))/(2*(alpha-1)*q)
    
    # find yp
    if q == 1:
        yp = alpha*zf/(1+zf*(alpha-1))
    elif q == 0:
        yp = zf
    else:
        yp = alpha*xp/(1+xp*(alpha-1))
    
    # find RRmin
    RRmin = ((xd-yp)/(yp-xp))

    # find Xf
    if q == 1:
        xf = zf
    else:
        xf = (xd/(1+RR)+zf/(q-1))/(q/(q-1)-RR/(1+RR))

    # find Yf
    yf = (xd+xf*RR)/(1+RR)

    # caculate stage matrix 
    # Initialize lists to store stage values
    stages = []
    xi_values = [xd]  # Start with x0 = xD
    yi_values = [xd]  # Start with y0 = xD

    # Iteratively calculate stage values
    for i in range(999):  
        # Calculate xi
        xi = xi_values[-1] - eta * (xi_values[-1] - yi_values[-1] / (alpha - yi_values[-1] * (alpha - 1)))
        # Calculate yi based on the condition
        if xi > xf:
            yi = (xd + xi * RR) / (1 + RR)
        else:
            yi = ((yf - xb) / (xf - xb)) * (xi - xb) + xb
        
        # Append to lists
        xi_values.append(xi)
        yi_values.append(yi)
        stages.append(i+1)
        
        # Break if yi is less than xb
        if yi < xb:
            break

    stages, xi_values, yi_values



    # find N (number of stages)--not filled in for now (need stage matrix)
    # Given the stages list, we'll calculate the number of stages based on the conditions
    def calculate_stages(stages, RR, RRmin, xb):
        # Check if RR is greater than RRmin
        if RR > RRmin:
            max_stage = max(stages)
            # Get the value from the xi_values list corresponding to the max stage
            xi_max_stage_minus_1 = xi_values[max_stage - 1]
            xi_max_stage = xi_values[max_stage]
            # Calculate the number of stages using the provided formula
            num_stages = max_stage - 1 + (xi_max_stage_minus_1 - xb) / (xi_max_stage_minus_1 - xi_max_stage)
            return num_stages
        else:
            return "INF"

    # Calculate the number of stages
    num_stages = calculate_stages(stages, RR, RRmin, xb)
    num_stages

    # find Nf (feed stage)--not filled in for now (need stage matrix)
    def calculate_feed_stage(stages, xi_values, RR, RRmin, xf):
        # Check if RR is greater than RRmin
        if RR > RRmin:
            # Find the index of the value in xi_values closest to but not greater than xf
            index = next((i for i, x in enumerate(xi_values) if x > xf), None)
            
            if index is None:
                return "No feed stage found"
            
            # Get the stage from the stages list corresponding to the index and add 1
            feed_stage = stages[index] - 1
            return feed_stage
        else:
            return "INF"


    # Calculate the feed stage
    feed_stage = calculate_feed_stage(stages, xi_values, RR, RRmin, xf)
    feed_stage

    # Define the range of mole fractions
    x = np.linspace(0, 1, 1000)

    # Generate the equilibrium curve
    x = np.linspace(0, 1, 1000)
    y = alpha * x / (1 + (alpha - 1) * x)
    
    # Rectifying Operating Line--correct
    S_R = (yf - xd) / (xf - xd)
    I_R = xd - S_R * xd
    y_rectifying = S_R * x + I_R
    
    # Stripping Operating Line--correct
    S_R = (yf - xb) / (xf - xb)
    I_R = xb - S_R * xb
    y_stripping = S_R * x + I_R
    
    # Feed Line--correct
    y_feed = (q*x-zf)/(q-1)

    # 45 Line
    y_45 = x

    outputs = {
        "xi_values": xi_values,
        "yi_values": yi_values,
        "stages": stages,
        "xp": xp,
        "yp": yp,
        "RRmin": RRmin,
        "xf": xf,
        "yf": yf,
        "N": num_stages,
        "Nf": feed_stage
    }

    # Filter x and y values for the Rectifying Operating Line based on the updated conditions
    x_rect_filtered = [x_val for x_val, y_val_rect, y_val_strip, y_45_val in zip(x, y_rectifying, y_stripping, y_45) if y_val_rect > y_45_val and y_val_rect < y_val_strip]
    y_rect_filtered = [y_val_rect for y_val_rect, y_val_strip, y_45_val in zip(y_rectifying, y_stripping, y_45) if y_val_rect > y_45_val and y_val_rect < y_val_strip]

    # Filter x and y values for the Stripping Operating Line based on the updated conditions
    x_strip_filtered = [x_val for x_val, y_val_rect, y_val_strip, y_45_val in zip(x, y_rectifying, y_stripping, y_45) if y_val_strip > y_45_val and y_val_strip < y_val_rect]
    y_strip_filtered = [y_val_strip for y_val_rect, y_val_strip, y_45_val in zip(y_rectifying, y_stripping, y_45) if y_val_strip > y_45_val and y_val_strip < y_val_rect]

    # Filter x and y values for the Feed Line to ensure it is plotted only below the equilibrium curve and above 45
    x_feed_filtered = [x_val for x_val, y_val_feed, y_val_eq, y_45_val in zip(x, y_feed, y, y_45) if y_val_feed > y_45_val and y_val_feed < y_val_eq]
    y_feed_filtered = [y_val_feed for y_val_feed, y_val_eq, y_45_val in zip(y_feed, y, y_45) if y_val_feed > y_45_val and y_val_feed < y_val_eq]

    # Plotting
    plt.figure(figsize=(10, 8))

    # Plot Equilibrium Curve
    plt.plot(x, y, label="Equilibrium Curve", color="blue")

    # Plot the filtered Rectifying Operating Line segment
    plt.plot(x_rect_filtered, y_rect_filtered, label="Rectifying Operating Line", color="red")

    # Plot the filtered Stripping Operating Line segment
    plt.plot(x_strip_filtered, y_strip_filtered, label="Stripping Operating Line", color="green")

    # Plot Feed Line
    plt.plot(x_feed_filtered, y_feed_filtered, label="Feed Line", color="purple", linestyle="--")

    # Plot 45-degree line
    plt.plot(x, y_45, label="45-degree line", color="black")

    # Plot the stages and the lines connecting them to the equilibrium curve
    for i in range(1, len(stages)+1):
        # Plot the stage point
        plt.scatter(xi_values[i], yi_values[i], color='black')
        # Vertical line up to equilibrium curve
        y_intersect = alpha * xi_values[i] / (1 + (alpha - 1) * xi_values[i])
        plt.plot([xi_values[i], xi_values[i]], [yi_values[i], y_intersect], color='grey')
        # Horizontal line to the next stage
        plt.plot([xi_values[i], xi_values[i-1]], [y_intersect, y_intersect], color='grey')
    
    # Plot the last stage point
    plt.scatter(xi_values[0], yi_values[0], color='black')

    # Labels, title, and other aesthetics
    plt.xlabel("x (Liquid phase mole fraction)")
    plt.ylabel("y (Vapor phase mole fraction)")
    plt.title("McCabe-Thiele Diagram")
    plt.legend()
    plt.grid(True)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.show()


    
    # Return the outputs dictionary
    return outputs


import tkinter as tk
from tkinter import ttk

def calculate():
    # Retrieve values from input fields
    alpha = float(alpha_entry.get())
    zf = float(zf_entry.get())
    q = float(q_entry.get())
    xd = float(xd_entry.get())
    xb = float(xb_entry.get())
    RR = float(RR_entry.get())
    eta = float(eta_entry.get())

    # Call the mccabe_thiele_plot function
    outputs = mccabe_thiele_plot(alpha, zf, q, xd, xb, RR, eta)

    # Clear the table (Treeview)
    for row in table.get_children():
        table.delete(row)

    # Populate the Treeview with data
    for stage, xi, yi in zip(outputs["stages"], outputs["xi_values"], outputs["yi_values"]):
        table.insert("", "end", values=(stage, xi, yi))

    # Display the other outputs in the text box
    output_text.delete(1.0, tk.END)  # Clear previous outputs
    for label, value in outputs.items():
        if label not in ["xi_values", "yi_values", "stages"]:
            output_text.insert(tk.END, f"{label}: {value}\n")

# Set up the main window
root = tk.Tk()
root.title("McCabe-Thiele Plot")

# Input labels and entry boxes
labels = ["alpha", "zf", "q", "xd", "xb", "RR", "eta"]
entries = [alpha_entry, zf_entry, q_entry, xd_entry, xb_entry, RR_entry, eta_entry] = [ttk.Entry(root) for _ in labels]

for i, label in enumerate(labels):
    ttk.Label(root, text=label).grid(row=i, column=0, padx=10, pady=5)
    entries[i].grid(row=i, column=1, padx=10, pady=5)

# Calculate button
ttk.Button(root, text="Calculate", command=calculate).grid(row=len(labels), column=0, columnspan=2, pady=10)

# Output text box
output_text = tk.Text(root, wrap=tk.WORD, height=10, width=50)
output_text.grid(row=len(labels)+1, column=0, columnspan=2, pady=10)

# Treeview for table display
table = ttk.Treeview(root, columns=("Stages", "xi_values", "yi_values"), show="headings")
table.heading("Stages", text="Stages")
table.heading("xi_values", text="xi_values")
table.heading("yi_values", text="yi_values")
table.grid(row=len(labels)+2, column=0, columnspan=2, pady=10)

root.mainloop()
