import sys
import numpy as np
from math import sqrt

from scipy.optimize import minimize

from matplotlib import pyplot as plt
plt.style.use("./mystyle.mplstyle")
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import cm
from matplotlib import colors
from matplotlib import rc
rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern'],'size':14})

from mpl_toolkits.axes_grid1 import Divider, Size

from potpaths_crit import h1path, h2path, a2path, spath, thetapath

lightgrey, grey, darkgrey, red = '#CCD7DB','#96ABB3','#586D75','#C42021'
vir = cm.get_cmap('viridis')
# print(colors.to_hex(vir(0.5)))
# sys.exit()


def pot(h1, h2, a2, s):
    return (29071. * a2**2 + 0.114897 * a2**4 - 30.7127 * a2 * h1 + 27490.4 * h1**2 + 0.229795 * a2**2 * h1**2 + 0.114897 * h1**4 - 60000. * h1 * h2 + 29071. * h2**2 + 0.229795 * a2**2 * h2**2 - 0.100697 * h1**2 * h2**2 + 0.114897 * h2**4 - 4809.53 * s**2 + 0.375 * a2**2 * s**2 - 0.0575 * a2 * h1 * s**2 + 0.125 * h1**2 * s**2 + 0.375 * h2**2 * s**2 + 0.25 * s**4)


nx,ny = 300, 300
# print(pot_r1_S(*[100.,100.]))
potEWmin = minimize(lambda x: pot(*x),[100., 100., 0., 0.])
potSmin = minimize(lambda x: pot(*x),[0., 0., 0., 100.])
Vmin = potEWmin.fun
vh1c = potEWmin.x[0]
vh2c = potEWmin.x[1]
tbc = 1.
vr1c = vh1c/sqrt(1+tbc**2) + vh2c/sqrt(1+tbc**(-2))
va2c = potEWmin.x[2]
vsc = potSmin.x[3]
#print(vh1c)
#print(vh2c)
#print(va2c)
#print(tbc)
#print(vr1c)
#print(vsc)

r1path = h1path/sqrt(1+tbc**2) + h2path/sqrt(1+tbc**(-2))
r2path = h1path/sqrt(1+tbc**(-2)) - h2path/sqrt(1+tbc**2)
#print(h1path[0])
#print(h2path[0])
#print(r1path[0])
#print(r2path[0])


barrier = (pot(h1path,h2path,a2path,spath)-Vmin)/(1e7)

def pot_r1_S(h1, h2, s, a2c=va2c):
    return pot(h1,h2,a2c,s)
    
# def pot_r1_S(r1, S):
#     return pot(r1,0.5489951542860402,S,3.141620544211027)

r1range = np.linspace(-np.round(0.1*vr1c,10),np.round(1.2*vr1c,10),num=nx)
srange = np.linspace(-np.round(0.1*vsc,10),np.round(1.2*vsc,10),num=ny)
r1mesh, smesh = np.meshgrid(r1range,srange)

Zpot = np.nan_to_num( np.log10(pot_r1_S(r1mesh/(1./sqrt(1+tbc**2) + 1./tbc/sqrt(1+tbc**(-2))), r1mesh/(tbc/sqrt(1+tbc**2) + 1./sqrt(1+tbc**(-2))), smesh) - Vmin) , nan=0.)

zmin, zmax = np.amin(Zpot), np.amax(Zpot)
#print(zmin)
#print(zmax)

# plot & subplot dimensions
hh = [Size.Fixed(1.), Size.Fixed(4.), Size.Fixed(1)]
vv = [Size.Fixed(0.7), Size.Fixed(1.), Size.Fixed(0.1), Size.Fixed(1.), Size.Fixed(0.1), Size.Fixed(1.),Size.Fixed(0.1), Size.Fixed(4.), Size.Fixed(0.4)]

# figure
wid = sum([i.get_size(0)[1] for i in hh])
hei = sum([i.get_size(0)[1] for i in vv])
fig = plt.figure(figsize=[wid,hei])

divider = Divider(fig, (0.0, 0.0, 1., 1.), hh, vv, aspect=False)

# subplots
axes = []

# axes
axes.append(fig.add_axes( (0.0, 0.0, 1., 1.), axes_locator=divider.new_locator(nx=1, ny=7)))
axes.append(fig.add_axes( (0.0, 0.0, 1., 1.), axes_locator=divider.new_locator(nx=1, ny=5)))
axes.append(fig.add_axes( (0.0, 0.0, 1., 1.), axes_locator=divider.new_locator(nx=1, ny=3)))
axes.append(fig.add_axes( (0.0, 0.0, 1., 1.), axes_locator=divider.new_locator(nx=1, ny=1)))

# plot
axes[0].contour(r1mesh, smesh, Zpot, levels=10, cmap='viridis', vmin=zmin, vmax=zmax)
axes[0].plot(r1path, spath,'--k')
axes[0].scatter(r1path[0], spath[0], marker='*', s=50, c=colors.to_hex(vir(0.5)),zorder=5)
axes[0].scatter(r1path[-1], spath[-1], marker='*', s=50, c=colors.to_hex(vir(0.5)),zorder=5)

axes[1].plot(r1path, r2path,'--k')
axes[1].scatter(r1path[0], r2path[0], marker='*', s=50, c=colors.to_hex(vir(0.5)),zorder=5)
axes[1].scatter(r1path[-1], r2path[-1], marker='*', s=50, c=colors.to_hex(vir(0.5)),zorder=5)

axes[2].plot(r1path, barrier,'--k')
axes[2].scatter(r1path[0], barrier[0], marker='*', s=50, c=colors.to_hex(vir(0.5)),zorder=5)
axes[2].scatter(r1path[-1], barrier[-1], marker='*', s=50, c=colors.to_hex(vir(0.5)),zorder=5)

axes[3].plot(r1path, thetapath*100,'--k')
axes[3].scatter(r1path[0], thetapath[0]*100., marker='*', s=50, c=colors.to_hex(vir(0.5)),zorder=5)
axes[3].scatter(r1path[-1], thetapath[-1]*100., marker='*', s=50, c=colors.to_hex(vir(0.5)),zorder=5)


xmin, xmax = -np.round(0.1*vr1c,10),np.round(1.2*vr1c,10)
ymin, ymax = -np.round(0.1*vsc,10),np.round(1.2*vsc,10)

# settings
for i in range(4):

    axes[i].xaxis.set_major_locator(MultipleLocator(40.))
    axes[i].xaxis.set_minor_locator(MultipleLocator(4.))
    axes[i].set_xlim(xmin, xmax)
    axes[i].tick_params(axis='both',which='both',direction='in')
    
    axes[i].plot([xmin, xmax],[0.,0.], ls='dotted', lw=1., c=grey)
    # axes[i].plot([xmin, xmax],[vsc,vsc], ls='dotted', lw=0.5, c=grey)
    axes[i].plot([0.,0.],[ymin, ymax], ls='dotted', lw=1., c=grey)

    if i!=0: axes[i].plot([vr1c, vr1c],[ymin, ymax], ls='dotted', lw=1., c=grey)

    if i!=3: axes[i].xaxis.set_ticklabels([])

#axes[0].text(4.,2.,
#r'$t_\beta=1$'
##'\n'
##r'$m_{a_1}=90\,\mathrm{GeV}$'
##'\n'
##r'$v_S=120\,\mathrm{GeV}$'
#'\n'
#r'$\lambda_{\beta}=0.5$'
#'\n'
#r'$\lambda_S=1$'
#'\n'
#r'$\lambda_{S_3}^I=0.115$'
#'\n'
#r'$m_s=130\,\mathrm{GeV}$'
##'\n'
##r'$\sin\theta_a=-0.047$'
#, fontsize=12)

axes[0].yaxis.set_major_locator(MultipleLocator(20.))
axes[0].yaxis.set_minor_locator(MultipleLocator(2.))
axes[0].set_ylim(bottom=ymin,top=ymax)
axes[0].set_title(r'$\log\!\left((V-V_{\mathrm{min}})\,\mathrm{[GeV]^4}\right);\,T_C=80.06\,\mathrm{GeV}$', pad=10., fontsize=16)
#axes[0].set_title(r'$\log\!\left((V-V_{\mathrm{min}})\,\,\mathrm{[GeV]^4}\right);\,'+r'a_2={:.2f}$ GeV'.format(va2c), pad=10., fontsize=16)
axes[0].set_ylabel(r'$s$ [GeV]', labelpad=5., fontsize=16)
axes[0].text(0., vsc*1.04, r'\textbf{CPV}', fontsize=16, fontweight='bold', horizontalalignment='center')
axes[0].text(vr1c*1.04, 0., r'\textbf{EW}', fontsize=16, fontweight='bold', horizontalalignment='left',verticalalignment='center')

axes[1].yaxis.set_major_locator(MultipleLocator(1.))
axes[1].yaxis.set_minor_locator(MultipleLocator(0.2))
axes[1].set_ylim(bottom=-0.4,top=3.4)
axes[1].set_ylabel(r'$\rho_2$ [GeV]', labelpad=18., fontsize=16)

axes[2].yaxis.set_major_locator(MultipleLocator(1.))
axes[2].yaxis.set_minor_locator(MultipleLocator(0.2))
axes[2].set_ylim(bottom=-0.4,top=1.9)
axes[2].set_ylabel(r'$(V-V_{\mathrm{min}})$''\n'r'$\times 10^{-7}$[GeV]$^4$', labelpad=12, fontsize=12)
#axes[2].set_ylabel(r'$(V-V_{\mathrm{min}})$''\n'r'$\times 10^{-7}$[GeV]$^4$', labelpad=15, fontsize=12, rotation=0., verticalalignment='center')
#axes[2].annotate("",
#            xy=(44., 0.), xycoords='data',
#            xytext=(44., 3.3), textcoords='data',
#            arrowprops=dict(arrowstyle="<->",
#                            color=darkgrey, lw=1,
#                            connectionstyle="arc3"),
#            )
#axes[2].text(48., 1.65, r"$\Delta V$", verticalalignment='center', c=darkgrey)


axes[3].yaxis.set_major_locator(MultipleLocator(1))
axes[3].yaxis.set_minor_locator(MultipleLocator(0.2))
axes[3].set_ylim(bottom=-2.6,top=0.4)
axes[3].set_ylabel(r'$\theta\times 100$', labelpad=6., fontsize=16)
axes[3].set_xlabel(r'$\rho_1$ [GeV]', labelpad=5., fontsize=16)
axes[3].annotate("",
            xy=(0., 0.), xycoords='data',
            xytext=(0., thetapath[-1]*100+0.05), textcoords='data',
            arrowprops=dict(arrowstyle="<->",
                            color=darkgrey, lw=1,
                            connectionstyle="arc3"),
            )
axes[3].text(4., thetapath[-1]*50, r"$\Delta\theta$", verticalalignment='center', c=darkgrey)

# plt.show()
plt.savefig('potential_plot_crit.pdf')
