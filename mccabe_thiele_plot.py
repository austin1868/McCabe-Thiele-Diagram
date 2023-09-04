import numpy as np
import matplotlib.pyplot as plt
import customtkinter as ctk
from customtkinter import cttk


def mccabe_thiele_plot(alpha, zf, q, xd, xb, RR, eta):
    # Find xp
    if q == 1:
        xp = zf
    elif q == 0:
        xp = zf / (alpha * zf * (alpha - 1))
    else:
        xp = (((alpha-1)*(zf+q)-alpha)+np.sqrt(((alpha-1)*(zf+q)-alpha)**2+4*zf*(alpha-1)*q))/(2*(alpha-1)*q)
    
    # Find yp
    if q == 1:
        yp = alpha*zf/(1+zf*(alpha-1))
    elif q == 0:
        yp = zf
    else:
        yp = alpha*xp/(1+xp*(alpha-1))

    # Calculate RRmin
    RRmin = (xd-yp) / (yp-xp)

    # Find Xf
    if q == 1:
        xf = zf
    else:
        xf = (xd/(1+RR)+zf/(q-1)) / (q/(q-1)-RR/(1+RR))

    # Find Yf
    yf = (xd+xf*RR) / (1+RR)

    # Calculate stage matrix
    stages = []
    xi_values = [xd]  # Start with x0 = xD
    yi_values = [xd]  # Start with y0 = xD
    for i in range(999):
        xi = xi_values[-1] - eta * (xi_values[-1] - yi_values[-1] / (alpha - yi_values[-1] * (alpha - 1)))
        if xi > xf:
            yi = (xd + xi * RR) / (1 + RR)
        else:
            yi = ((yf - xb) / (xf - xb)) * (xi - xb) + xb

        xi_values.append(xi)
        yi_values.append(yi)
        stages.append(i+1)

        if yi < xb:
            break

    # Function to calculate number of stages
    def calculate_stages():
        if RR <= RRmin:
            return "INF"
        max_stage = max(stages)
        xi_max_stage_minus_1 = xi_values[max_stage - 1]
        xi_max_stage = xi_values[max_stage]
        return max_stage - 1 + (xi_max_stage_minus_1 - xb) / (xi_max_stage_minus_1 - xi_max_stage)

    num_stages = calculate_stages()

    # Function to calculate feed stage
    def calculate_feed_stage():
        if RR <= RRmin:
            return "INF"
        index = next((i for i, x in enumerate(xi_values) if x > xf), None)
        if index is None:
            return "No feed stage found"
        return stages[index] - 1

    feed_stage = calculate_feed_stage()

    # Define mole fractions range
    x = np.linspace(0, 1, 1000)

    # Equilibrium curve
    y = alpha * x / (1 + (alpha - 1) * x)
    
    # Rectifying Operating Line
    S_R_rect = (yf - xd) / (xf - xd)
    I_R_rect = xd - S_R_rect * xd
    y_rectifying = S_R_rect * x + I_R_rect

    # Stripping Operating Line
    S_R_strip = (yf - xb) / (xf - xb)
    I_R_strip = xb - S_R_strip * xb
    y_stripping = S_R_strip * x + I_R_strip

    # Feed Line
    y_feed = (q * x - zf) / (q - 1)

    # 45 Line
    y_45 = x

    # Filter lines based on conditions
    filter_func = lambda val, *conditions: all(condition(val) for condition in conditions)
    x_rect_filtered = [x_val for x_val in x if filter_func(x_val, 
                                                           lambda val: S_R_rect * val + I_R_rect > val,
                                                           lambda val: S_R_rect * val + I_R_rect < S_R_strip * val + I_R_strip)]
    y_rect_filtered = [S_R_rect * val + I_R_rect for val in x_rect_filtered]

    x_strip_filtered = [x_val for x_val in x if filter_func(x_val, 
                                                            lambda val: S_R_strip * val + I_R_strip > val,
                                                            lambda val: S_R_strip * val + I_R_strip < S_R_rect * val + I_R_rect)]
    y_strip_filtered = [S_R_strip * val + I_R_strip for val in x_strip_filtered]

    x_feed_filtered = [x_val for x_val in x if filter_func(x_val, 
                                                           lambda val: (q * val - zf) / (q - 1) > val,
                                                           lambda val: (q * val - zf) / (q - 1) < alpha * val / (1 + (alpha - 1) * val))]
    y_feed_filtered = [(q * val - zf) / (q - 1) for val in x_feed_filtered]

    # Plotting
    plt.figure(figsize=(10, 8))
    plt.plot(x, y, label="Equilibrium Curve", color="blue")
    plt.plot(x_rect_filtered, y_rect_filtered, label="Rectifying Operating Line", color="red")
    plt.plot(x_strip_filtered, y_strip_filtered, label="Stripping Operating Line", color="green")
    plt.plot(x_feed_filtered, y_feed_filtered, label="Feed Line", color="purple", linestyle="--")
    plt.plot(x, y_45, label="45-degree line", color="black")

    # Plot stages and connecting lines
    for i in range(1, len(stages) + 1):
        plt.scatter(xi_values[i], yi_values[i], color='black')
        y_intersect = alpha * xi_values[i] / (1 + (alpha - 1) * xi_values[i])
        plt.plot([xi_values[i], xi_values[i]], [yi_values[i], y_intersect], color='grey')
        plt.plot([xi_values[i], xi_values[i-1]], [y_intersect, y_intersect], color='grey')

    plt.scatter(xi_values[0], yi_values[0], color='black')
    plt.xlabel("x (Liquid phase mole fraction)")
    plt.ylabel("y (Vapor phase mole fraction)")
    plt.title("McCabe-Thiele Diagram")
    plt.legend()
    plt.grid(True)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.show()

    return {
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

def calculate():
    output_text.delete(1.0, ctk.END)
    print(f"Alpha entry value before error check: '{alpha_entry.get()}'")
    try:
        alpha_val = alpha_entry.get()
        if not alpha_val:
            output_text.insert(ctk.END, "Alpha input is empty.\n")
            return
        alpha = float(alpha_val)
    except ValueError:
        output_text.insert(ctk.END, f"Problem with alpha input: {alpha_val}\n")
        return

    try:
        zf = float(zf_entry.get())
    except ValueError:
        output_text.insert(ctk.END, f"Problem with zf input: {zf_entry.get()}\n")
        return

    try:
        q = float(q_entry.get())
    except ValueError:
        output_text.insert(ctk.END, f"Problem with q input: {q_entry.get()}\n")
        return

    try:
        xd = float(xd_entry.get())
    except ValueError:
        output_text.insert(ctk.END, f"Problem with xd input: {xd_entry.get()}\n")
        return

    try:
        xb = float(xb_entry.get())
    except ValueError:
        output_text.insert(ctk.END, f"Problem with xb input: {xb_entry.get()}\n")
        return

    try:
        RR = float(RR_entry.get())
    except ValueError:
        output_text.insert(ctk.END, f"Problem with RR input: {RR_entry.get()}\n")
        return

    try:
        eta = float(eta_entry.get())
    except ValueError:
        output_text.insert(ctk.END, f"Problem with eta input: {eta_entry.get()}\n")
        return
    
    try:
        outputs = mccabe_thiele_plot(alpha, zf, q, xd, xb, RR, eta)
    except Exception as e:
        output_text.insert(ctk.END, f"Error in mccabe_thiele_plot: {str(e)}\n")
        return
    
    # Clear the table
    for row in table.get_children():
        table.delete(row)

    # Populate table with data
    for stage, xi, yi in zip(outputs["stages"], outputs["xi_values"], outputs["yi_values"]):
        table.insert("", "end", values=(stage, xi, yi))

    # Display other outputs
    output_text.delete(1.0, ctk.END)
    for label, value in outputs.items():
        if label not in ["xi_values", "yi_values", "stages"]:
            output_text.insert(ctk.END, f"{label}: {value}\n")

def on_slider_change(value, entry):
    """Update the entry when the slider changes."""
    global updating_entry
    print(f"Slider changed to: {value}")  # Debugging line
    updating_entry = True
    entry.delete(0, ctk.END)
    entry.insert(0, str(round(float(value), 2)))
    root.update_idletasks()
    updating_entry = False

def on_entry_change(entry_var, slider):
    """Update the slider when the entry changes."""
    print(f"Entry changed to: {entry_var.get()}")  # Debugging line
    if not updating_entry:
        try:
            slider.set(float(entry_var.get()))
        except ValueError:
            pass

# GUI Setup
root = ctk.Tk()
root.title("McCabe-Thiele Plot")

# Dark Mode Colors
COLORS = {
    "bg": '#2E2E2E',
    "entry_bg": '#333333',
    "entry_fg": '#FFFFFF',
    "text": '#FFFFFF',
    "button": '#555555',
    "button_text": '#FFFFFF',
    "input_text": '#2E2E2E'
}

# Configure window and styles
root.configure(bg=COLORS["bg"])
style = cttk.Style()
style.configure("TButton", background=COLORS["button"], foreground=COLORS["input_text"])
style.configure("TLabel", background=COLORS["bg"], foreground=COLORS["text"])
style.configure("TEntry", foreground=COLORS["input_text"], fieldbackground=COLORS["entry_bg"])
style.configure("Treeview", background=COLORS["entry_bg"], foreground=COLORS["entry_fg"], fieldbackground=COLORS["entry_bg"])
style.map('Treeview', background=[('selected', '#444444')], foreground=[('selected', 'white')])

# Input labels and entry boxes
label_font = ("Arial", 12, "bold")
labels = ["alpha", "zf", "q", "xd", "xb", "RR", "eta"]
entries = [alpha_entry, zf_entry, q_entry, xd_entry, xb_entry, RR_entry, eta_entry] = [cttk.Entry(root) for _ in labels]


# Input labels, entry boxes, and sliders
labels = ["alpha", "zf", "q", "xd", "xb", "RR", "eta"]
ranges = [(1, 10), (0, 1), (0, 1), (0, 1), (0, 1), (0, 10), (0, 1)]  # Example ranges for each slider; adjust as needed

entries = []
sliders = []

for i, (label, (min_val, max_val)) in enumerate(zip(labels, ranges)):
    cttk.Label(root, text=label, font=label_font).grid(row=i, column=0, padx=10, pady=5)
    
    entry_var = ctk.StringVar()
    entry = cttk.Entry(root, textvariable=entry_var)
    entry.grid(row=i, column=1, padx=10, pady=5)
    entries.append(entry)
    
    slider = cttk.Scale(root, from_=min_val, to=max_val, orient=ctk.HORIZONTAL, command=lambda value, entry=entry: on_slider_change(value, entry))
    slider.grid(row=i, column=2, padx=10, pady=5)
    sliders.append(slider)

    # Link slider and entry
    entry_var.trace_add("write", lambda *args, slider=slider, entry_var=entry_var: on_entry_change(entry_var, slider))

# Set default values for sliders and entries (this ensures they aren't empty on start)
defaults = [5, 0.5, 0.5, 0.5, 0.5, 5, 0.5]
for entry, slider, default in zip(entries, sliders, defaults):
    entry.insert(0, str(default))
    slider.set(default)

cttk.Button(root, text="Calculate", command=calculate).grid(row=len(labels), column=0, columnspan=3, pady=10)
output_text = ctk.Text(root, wrap=ctk.WORD, height=10, width=50, bg=COLORS["entry_bg"], fg=COLORS["entry_fg"])
output_text.grid(row=len(labels)+1, column=0, columnspan=3, pady=10)
table = cttk.Treeview(root, columns=("Stages", "xi_values", "yi_values"), show="headings")
table.heading("Stages", text="Stages")
table.heading("xi_values", text="xi_values")
table.heading("yi_values", text="yi_values")
table.grid(row=len(labels)+2, column=0, columnspan=3, pady=10)

root.mainloop()
