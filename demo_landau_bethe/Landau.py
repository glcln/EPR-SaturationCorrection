import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.integrate import simps

def landau_distribution(x, mu, c):
    """
    Calculate the Landau distribution function.
    
    Parameters:
    x (array-like): The input array of x-values where the distribution is evaluated.
    mu (float): The location parameter (mean).
    c (float): The scale parameter (sigma).
    
    Returns:
    array-like: The Landau distribution evaluated at the input x-values.
    """

    Landau_end = 100
    nPoint = Landau_end * 100
    result = 0

    for i in range(1,nPoint):
        result += np.exp(-i/100)*np.cos(i/100 * (x-mu)/c + 2*i/(100*np.pi) * np.log(i/(100*c)))  # OK!!

    return 1/(np.pi*c)*result

# Parameters
end_plot = 15
mu = 0  # Location parameter (mean)
c = 1   # Scale parameter (sigma)
size = 1000  # Number of random samples to generate

# Generate random samples from the Landau distribution
x = np.linspace(-4, end_plot, size)
y = landau_distribution(x, mu, c)
y2 = landau_distribution(x, mu, c+0.5)
y3 = landau_distribution(x, mu, c+1)

# Calculate the integral and normalize the PDF
integral = simps(y, x) 
y_normalized = y / integral 
integral2 = simps(y2, x) 
y2_normalized = y2 / integral2 
integral3 = simps(y3, x)  
y3_normalized = y3 / integral3 

# Integral after MPV
x_int = np.linspace(x[np.argmax(y)], end_plot, size)
x_int2 = np.linspace(x[np.argmax(y2)], end_plot, size)
x_int3 = np.linspace(x[np.argmax(y3)], end_plot, size)

value1 = simps(y, x_int)
value2 = simps(y2, x_int2)
value3 = simps(y3, x_int3)

print("Integral after max c=1: ", value1)
print("Integral after max c=1.5: ", value2)
print("Integral after max c=2: ", value3)

# Activate LaTeX text rendering
plt.rc('text', usetex=True)

# Plot the Landau distribution
fig, ax = plt.subplots()
plt.plot(x, y_normalized, color='blue', label=r'$\mu = 0, c = 1, \textrm{MPV} = ' + "{:.1f}".format(x[np.argmax(y)]) + r', \int_{\geq \textrm{MPV}} = $'+" {:.1f}".format(value1)+r' \%')
plt.plot(x, y2_normalized, color='red', label=r'$\mu = 0, c = 1.5, \textrm{MPV} = ' + "{:.1f}".format(x[np.argmax(y2)]) + r', \int_{\geq \textrm{MPV}} = $'+" {:.1f}".format(value2) + r' \%')
plt.plot(x, y3_normalized, color='green', label=r'$\mu = 0, c = 2, \textrm{MPV} = ' + "{:.1f}".format(x[np.argmax(y3)]) + r', \int_{\geq \textrm{MPV}} = $'+" {:.1f}".format(value3) + r' \%')
plt.title(r'$\textrm{Distribution de Landau}$', fontsize=14)
plt.xlabel(r'$x$', fontsize=12)
plt.ylabel(r'$\textrm{Densit\'e de Probabilit\'e}$', fontsize=12)
plt.legend()
ax.xaxis.set_label_coords(0.5,-0.07)
ax.yaxis.set_label_coords(-0.07,0.5)
ax.grid(which='minor', alpha=0.2)
ax.grid(which='major', alpha=0.5)

plt.savefig("landau_distribution.pdf")

