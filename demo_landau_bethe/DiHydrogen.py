import numpy as np
import matplotlib.pyplot as plt
from scipy.special import expi
import matplotlib.ticker as mticker

E = 13.6         # Fondamental state of hydrogen
a0 = 0.529       # Hydrogen radius
gamma = 0.577    # Euler-Mascheroni constant

R = np.linspace(0.01, 5, 1000)
Ei = np.zeros(2)

S = []
Sp = []
eps = []
v = []
J = []
K = []

Ep = []
Em = []
Dp = []
Dm = []

for i in range (len(R)):
    rho = R[i]/a0

    # exponential integral
    x = [-4 * rho, -2 * rho]
    Ei[0] = expi(x[0])
    Ei[1] = expi(x[1])

    S.append(np.exp(-rho) * (1 + rho + 1./3 * rho**2))
    Sp.append(np.exp(rho) * (1 - rho + 1./3 * rho**2))
    eps.append(-E - 2*E/rho*(1 - np.exp(-2*rho)*(1 + rho)))
    v.append(-E - 2*E*np.exp(-rho)*(1 + rho))
    J.append(2*E * (1./rho - np.exp(-2*rho)*(1./rho + 11./8 + 3./4*rho + 1./6*rho**2)))
    K.append(2./5 * E * (np.exp(-2*rho)*(25./8 - 23./4*rho - 3*rho**2 - 1./3*rho**3) + 6./rho*(S[i]**2*(gamma + np.log(rho)) + Sp[i]**2*Ei[0] - 2*S[i]*Sp[i]*Ei[1])) )
    
    Ep.append(2*E/rho + (2*(eps[i] - S[i]*v[i]) + J[i] - K[i])/(1 - S[i]**2) )
    Em.append(2*E/rho + (2*(eps[i] + S[i]*v[i]) + J[i] + K[i])/(1 + S[i]**2) )
    Dp.append(Ep[i]  + 2*E)
    Dm.append(Em[i] + 2*E)


# Activate LaTeX text rendering
plt.rc('text', usetex=True)
plt.rcParams.update({
    'text.usetex': True,
    'font.family': 'serif',
})

Ep = [x/E for x in Ep]
Em = [x/E for x in Em]
Dp = [x/E for x in Dp]
Dm = [x/E for x in Dm]
rho = [x/a0 for x in R]

# Plot
fig, ax = plt.subplots()

ax.plot(rho, Dp, color='blue', label=r'$\Delta E_{+}$')
ax.plot(rho, Dm, color='red', label=r'$\Delta E_{-}$')
plt.title(r'$\textrm{Dihydrogen - Energies curves}$', fontsize=14)
plt.xlabel(r'$\frac{R}{a_0}$',fontsize=15, fontweight='bold')
plt.ylabel(r'$\frac{E}{|E_I|}$',fontsize=15, fontweight='bold')
ax.set_ylim(-1, 4)
ax.grid(which='minor', alpha=0.2)
ax.grid(which='major', alpha=0.5)



plt.legend()

minpos = rho[Em.index(min(Em))]
plt.plot(minpos, min(Dm),marker='o', markerfacecolor='black', markeredgecolor='gray', markeredgewidth=2, markersize=10)
ax.text(minpos, min(Dm) + 0.3, r'$($'+r'${:.2f}$'.format(minpos)+r'$,$'+r'${:.2f}$'.format(min(Dm))+r'$)$', color='black', fontsize=12, horizontalalignment='center')

R_eq = "{:.2f}".format(minpos * a0)
E_eq = "{:.2f}".format(min(Em) * E + 2 * E)
textstr = rf"$R_{{eq}} = {R_eq}$ \AA" "\n" rf"$E_{{eq}} = {E_eq}$ eV"
ax.text(0.75, 0.75, textstr, transform=ax.transAxes, fontsize=12, verticalalignment='center', horizontalalignment='left', bbox=dict(facecolor='white', alpha=0.5))

plt.savefig("DiHydrogen.pdf")

print("Distance d'equilibre: R = ", "{:.2f}".format(minpos*a0), " Angstroms")
print("Energie d'equilibre: E = ", "{:.2f}".format(min(Em)*E+2*E), " eV")