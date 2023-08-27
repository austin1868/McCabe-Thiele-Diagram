import numpy as np
import matplotlib.pyplot as plt
alpha = 4
zf = 0.7
q = 0.4
xd = 0.95
xb = 0.1
RR = 1.3
eta = 1

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

    print(xi_values)
    print(yi_values)
    print(stages)

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

    print("N:", num_stages)

    # find Nf (feed stage)--not filled in for now (need stage matrix)
    def calculate_feed_stage(stages, xi_values, RR, RRmin, xf):
        # Check if RR is greater than RRmin
        if RR > RRmin:
            # Find the index of the value in xi_values closest to but not greater than xf
            index = next(i for i, x in enumerate(xi_values) if x > xf)
            # Get the stage from the stages list corresponding to the index and add 1
            feed_stage = stages[index] - 1
            return feed_stage
        else:
            return "INF"

    # Calculate the feed stage
    feed_stage = calculate_feed_stage(stages, xi_values, RR, RRmin, xf)
    feed_stage

    print("Nf:", feed_stage)

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

 
    # Plotting
    plt.figure(figsize=(10, 8))
    plt.plot(x, y, label="Equilibrium Curve", color="blue")
    plt.plot(x, y_rectifying, label="Rectifying Operating Line", color="red")
    plt.plot(x, y_stripping, label="Stripping Operating Line", color="green")
    plt.plot(x, y_feed, label="Feed Line", color="purple", linestyle="--")
    plt.plot(x, y_45, label="", color="black")
    
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
    plt.scatter(xi_values[-1], yi_values[-1], color='black')

    # Labels, title, and other aesthetics
    plt.xlabel("x (Liquid phase mole fraction)")
    plt.ylabel("y (Vapor phase mole fraction)")
    plt.title("McCabe-Thiele Diagram")
    plt.legend()
    plt.grid(True)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.show()

# To use the function, provide the required parameters and call the function:

mccabe_thiele_plot(alpha, zf, q, xd, xb, RR, eta)
