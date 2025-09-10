import seafreeze.seafreeze as sf
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

def fp_finder_func(t, p, m, ice_phase):
    """ Dummy argument to use as the callable for root_scalar """
    tup = np.empty((1,), dtype=object)
    tup[0] = (p, t, m)
    exp = sf.getProp(tup, phase='NaClaq').muw[0]/0.0180153 - sf.getProp(tup, phase=ice_phase).G[0]
    return exp

# Ih freezing points
Np, Nm = 51, 37
p_Ih = np.linspace(0.1, 251.1, Np)
m_Ih = np.linspace(0., 3.6, Nm)
Tf_Ih = np.zeros((Np, Nm))
for i in range(Np):
    for j in range(Nm):
        func = lambda t: fp_finder_func(t, p_Ih[i], m_Ih[j], ice_phase='Ih')
        Tf_Ih[i,j] = root_scalar(func, x0=270).root

# III freezing points
Np, Nm = 11, 37
p_III = np.linspace(200, 360, Np)
m_III = np.linspace(0., 3.6, Nm)
Tf_III = np.zeros((Np, Nm))
for i in range(Np):
    for j in range(Nm):
        func = lambda t: fp_finder_func(t, p_III[i], m_III[j], ice_phase='III')
        Tf_III[i,j] = root_scalar(func, x0=255).root

# III freezing points
Np, Nm = 31, 37
p_V = np.linspace(340, 640, Np)
m_V = np.linspace(0., 3.6, Nm)
Tf_V = np.zeros((Np, Nm))
for i in range(Np):
    for j in range(Nm):
        func = lambda t: fp_finder_func(t, p_V[i], m_V[j], ice_phase='V')
        Tf_V[i,j] = root_scalar(func, x0=255).root

# VI freezing points
Np, Nm = 26, 37
p_VI = np.linspace(600, 1000, Np)
m_VI = np.linspace(0, 3.6, Nm)
Tf_VI = np.zeros((Np, Nm))
for i in range(Np):
    for j in range(Nm):
        func = lambda t: fp_finder_func(t, p_VI[i], m_VI[j], ice_phase='VI')
        Tf_VI[i,j] = root_scalar(func, x0=280).root

# Plot the data to see how best to fit it
fig, axes = plt.subplots(2, 2, figsize=(12, 12), layout='compressed')
Ih_ax = axes[0,0]
III_ax = axes[0,1]
V_ax = axes[1,0]
VI_ax = axes[1,1]
ps = [p_Ih, p_III, p_V, p_VI]
ts = [Tf_Ih, Tf_III, Tf_V, Tf_VI]
ax_list = [Ih_ax, III_ax, V_ax, VI_ax]
phase_list = ['Ih', 'III', 'V', 'VI']
for i in range(len(ps)):
    ax = ax_list[i]
    cm = ax.contourf(ps[i], m_VI, ts[i], levels=np.arange(240, 301, 3), cmap='plasma')
    ax.set_xlabel('Pressure (MPa)')
    ax.set_ylabel('Salt molality (mol/kg)')
    ax.set_title(phase_list[i])
fig.colorbar(cm, 'Freezing point (K)')
plt.savefig('freezing_points.png', format='png')
