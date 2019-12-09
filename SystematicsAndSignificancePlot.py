import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

angular_size = 1.0
energy_threshold = 200
bkg_rate = 5000.*pow(energy_threshold/200.,-2.0) # per hour per 2 degree FoV (gamma-ray cut at MSCW=1 and MSCL=1)
crab_rate = 1e-7*3.75*pow(energy_threshold/1000.,-2.2);

def staterr2hour(x):
        x = x*0.01
        return 1./(pow(x,2)*bkg_rate*angular_size*angular_size/(2.*2.))
def hour2staterr(x):
        return 1./pow(x*bkg_rate*angular_size*angular_size/(2.*2.),0.5)
def Sensitivity(p_stat,p_syst):
        return pow(p_stat*5.,2)/2. + pow( pow(pow(p_stat*5.,2)/2.,2) + pow(p_stat*5.,2) + pow(p_syst*5.,2) ,0.5)

exposure_hours = [0.5, 2.0, 10.0, 100.0]
legends = ['0.5 hrs', '2 hrs', '10 hrs', '100 hrs']
stat_err = []
for i in range(0,len(exposure_hours)):
    stat_err += [hour2staterr(exposure_hours[i])*100.]

delta = 0.01
x = np.arange(0, stat_err[0], delta) # relative stat err
y = np.arange(0, 2.*stat_err[0], delta) # relative syst err
X, Y = np.meshgrid(x, y)
Z = Sensitivity(X*0.01,Y*0.01)

fig, ax = plt.subplots()
CS1 = ax.contour(Y, X, Z)
plt.plot(Y, 0.5*Y, color='r')
for i in range(0,len(exposure_hours)):
    plt.scatter(2.*stat_err[i], stat_err[i], color='red')
    plt.text(2.*stat_err[i]+0, stat_err[i]+0, legends[i], fontsize=9)
ax.clabel(CS1, inline=1, fontsize=10)
ax.set_title('Lowest signal-to-background ratio detected at $5\sigma$')
ax.axis('on')
ax.set_ylabel("Statistical error ($\Delta_{stat}/N$) %")
ax.set_xlabel('Systematic error ($\Delta_{syst}/N$) %')
#ax2 = ax.twinx()
#new_tick_locations = np.array([1,2,5])
#ax2.set_ylim(ax.get_ylim())
#ax2.set_yticks(new_tick_locations)
#ax2.set_yticklabels(staterr2hour(new_tick_locations))
#ax2.set_ylabel('exposure hours')
plt.savefig("output_plots/SysytematicsAndSignificance_%sGeV.png"%(energy_threshold))

delta = 0.001
exposure = 10.
x = np.arange(0.1, 1.0, delta) # angular size
N_CR = exposure*bkg_rate*x*x/(2.*2.)
X = 1./pow(N_CR,0.5)
Y0 = 0
Z0 = Sensitivity(X,Y0*0.01)
Y1 = 4
Z1 = Sensitivity(X,Y1*0.01)

fig, ax = plt.subplots()
plt.plot(x, Z0, color='b', label='$\Delta_{syst}/N$=0%')
plt.plot(x, Z1, color='r', label='$\Delta_{syst}/N$=4%')
ax.set_title('For an exposure of %s hours at %s GeV'%(exposure,energy_threshold))
ax.axis('on')
ax.set_xlabel('Angular size [degree]')
ax.set_ylabel('Lowest signal-to-background ratio detected at $5\sigma$')
ax.legend(loc='best')
plt.savefig("output_plots/AngularSizeAndSignificance_%sGeV.png"%(energy_threshold))

delta = 1.0
exposure = 10.
angular_size = 1.0
x = np.arange(200, 4000, delta) # energy threshold
bkg_rate = 5000.*pow(x/200.,-2.0) # per hour per 2 degree FoV (gamma-ray cut at MSCW=1 and MSCL=1)
N_CR = exposure*bkg_rate*angular_size*angular_size/(2.*2.)
X = 1./pow(N_CR,0.5)
Y0 = 0
Z0 = Sensitivity(X,Y0*0.01)
Y1 = 4
Z1 = Sensitivity(X,Y1*0.01)

fig, ax = plt.subplots()
plt.plot(x, Z0, color='b', label='$\Delta_{syst}/N$=0%')
plt.plot(x, Z1, color='r', label='$\Delta_{syst}/N$=4%')
ax.set_title('For an exposure of %s hours with $%s^{\circ}$ FoV'%(exposure,angular_size))
ax.axis('on')
ax.set_xlabel('Energy [GeV]')
ax.set_ylabel('Lowest signal-to-background ratio detected at $5\sigma$')
ax.set_xscale('log')
ax.legend(loc='best')
plt.savefig("output_plots/EnergyThresholdAndSignificance.png")


syst_a = 0.0
syst_b = 0.02
crab_percent = 0.2
eff_area = 100.*100.
angular_size_0 = 0.1
angular_size_1 = 1.0
delta = 0.1
x = np.arange(0.1, 100, delta) # exposure hours
bkg_rate = 5000.*pow(200./200.,-2.0) # per hour per 2 degree FoV (gamma-ray cut at MSCW=1 and MSCL=1)
N_SR = x*3600.*eff_area*crab_rate*crab_percent
N_SR_stat_err = pow(N_SR,0.5)
N_CR_0 = x*bkg_rate*angular_size_0*angular_size_0/(2.*2.)
N_CR_0_stat_err = pow(N_CR_0,0.5)
N_CR_0_syst_err_a = N_CR_0*syst_a
N_CR_0_syst_err_b = N_CR_0*syst_b
FoM_0_a = N_SR/pow(N_SR_stat_err*N_SR_stat_err+N_CR_0_stat_err*N_CR_0_stat_err+N_CR_0_syst_err_a*N_CR_0_syst_err_a,0.5)
FoM_0_b = N_SR/pow(N_SR_stat_err*N_SR_stat_err+N_CR_0_stat_err*N_CR_0_stat_err+N_CR_0_syst_err_b*N_CR_0_syst_err_b,0.5)
N_CR_1 = x*bkg_rate*angular_size_1*angular_size_1/(2.*2.)
N_CR_1_stat_err = pow(N_CR_1,0.5)
N_CR_1_syst_err_a = N_CR_1*syst_a
N_CR_1_syst_err_b = N_CR_1*syst_b
FoM_1_a = N_SR/pow(N_SR_stat_err*N_SR_stat_err+N_CR_1_stat_err*N_CR_1_stat_err+N_CR_1_syst_err_a*N_CR_1_syst_err_a,0.5)
FoM_1_b = N_SR/pow(N_SR_stat_err*N_SR_stat_err+N_CR_1_stat_err*N_CR_1_stat_err+N_CR_1_syst_err_b*N_CR_1_syst_err_b,0.5)
sigma_5 = 5.*x/x

fig, ax = plt.subplots()
plt.plot(x, FoM_0_a, color='b', label='point source $\\theta=%s^{\circ}$, syst = %0.0f%%'%(angular_size_0,syst_a*100.))
plt.plot(x, FoM_1_a, color='r', label='extended source $\\theta=%s^{\circ}$, syst = %0.0f%%'%(angular_size_1,syst_a*100.))
plt.plot(x, FoM_0_b, linestyle='-.', color='b', label='point source $\\theta=%s^{\circ}$, syst = %0.0f%%'%(angular_size_0,syst_b*100.))
plt.plot(x, FoM_1_b, linestyle='-.', color='r', label='extended source $\\theta=%s^{\circ}$, syst = %0.0f%%'%(angular_size_1,syst_b*100.))
plt.plot(x, sigma_5, color='grey')
ax.set_title('For a source of %s Crab unit'%(crab_percent))
ax.axis('on')
ax.set_xlabel('exposure hours')
ax.set_ylabel('significance')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(loc='best')
plt.savefig("output_plots/SignificanceVsExposure.png")
