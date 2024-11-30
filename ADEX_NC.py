import matplotlib.pyplot as plt 
import numpy as np  
from brian2 import *
import pylab as p
import matplotlib.ticker as ticker


plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = "24"


fig=plt.figure(num=None, figsize=(16, 8), dpi=80, facecolor='w', edgecolor='k')

ax3 = plt.subplot2grid((2,2), (0,0), rowspan=2) #, colspan=3)
ax1 = plt.subplot2grid((2,2), (0,1))#, colspan=2)
ax2 = plt.subplot2grid((2,2), (1, 1))
ax1.locator_params(nbins=4)
ax2.locator_params(nbins=4)
ax3.locator_params(nbins=4)


filename='ADEX_NC'



Par={'gl': 10e-09, 'El': -0.065, 'Vt': -0.055, 'Tsyn': 0.005, 'b': 20.e-12, 'Ee': 0.0, 'Cm': 200e-12, 'Dt': 0.005, 'a': 2.0e-9, 'Is':12.0e-11 , 'refractory': 0.005, 'Ei': -0.08, 'Vreset': -0.052, 'tau_w': 0.5}

with open(filename+'.txt', 'w') as f:
    print(Par, file=f)


def bin_array(array, BIN, time_array):
    N0 = int(BIN/(time_array[1]-time_array[0]))
    N1 = int((time_array[-1]-time_array[0])/BIN)
    return array[:N0*N1].reshape((N1,N0)).mean(axis=1)


start_scope()
DT=0.1
defaultclock.dt = DT*ms

N2 = 1#8000
#tau = 10*ms
#v0_max = 3.
TotTime=2000
duration = TotTime*ms



eqs='''
dv/dt = (-GsynE*(v-Ee)-GsynI*(v-Ei)-gl*(v-El)+ gl*Dt*exp((v-Vt)/Dt)-w + Is)/Cm : volt (unless refractory)
dw/dt = (a*(v-El)-w)/tau_w:ampere
dGsynI/dt = -GsynI/Tsyn : siemens
dGsynE/dt = -GsynE/Tsyn : siemens
Is:ampere
Cm:farad
gl:siemens
El:volt
tau_w:second
Dt:volt
Vt:volt
Ee:volt
Ei:volt
a:siemens
Tsyn:second
'''#% neuron_params



######################################################################################
######################################################################################

Vinit1=-60.e-3
Winint1=0.e-10 #Par['ga']/(1.0+exp((Par['Vc']-Vinit1)/Par['Vd']))
initC=[Vinit1,Winint1] #[-42.25*mV,3.5*nS ]#[-48.68*mV,6.5550*nS ]

# Population 2 - RS

b2 = Par['b']*amp
Vre= Par['Vreset']*volt
G2 = NeuronGroup(N2, eqs, threshold='v > -40.0*mV', reset='v = Vre; w += b2', refractory='5*ms',  method='heun')
G2.v = initC[0]*volt
G2.w = initC[1]*amp
G2.GsynI=0.0*nS
G2.GsynE=0.0*nS
G2.Cm     = Par['Cm']*farad
G2.gl     = Par['gl']*siemens
G2.El     = Par['El']*volt
G2.Vt     = Par['Vt']*volt
G2.Dt     = Par['Dt']*volt
G2.tau_w  = Par['tau_w']*second
G2.Is     = Par['Is']*amp  #2.50*nA #[0.0 for i in range(N2)]*nA
G2.Ee     = Par['Ee']*volt
G2.Ei     = Par['Ei']*volt
G2.Tsyn   = Par['Tsyn']*second
G2.a      = Par['a']*siemens

# external drive--------------------------------------------------------------------------

#P_ed_exc=PoissonGroup(8000, 4*Hz)
#P_ed_inh=PoissonGroup(2000, 4*Hz)

# connections-----------------------------------------------------------------------------

#Qi=5.*nS
#Qe=1.5*nS
#prbC=0.0
#S_edG2_ex = Synapses(P_ed_exc, G2, on_pre='GsynE_post+=Qe')
#S_edG2_ex.connect(p=prbC)

#S_edG2_in = Synapses(P_ed_inh, G2, on_pre='GsynI_post+=Qi')
#S_edG2_in.connect(p=prbC)

#S_edG2 = Synapses(P_ed, G2, on_pre='GsynE_post+=Qe')
#S_edG2.connect(p=1)


Vtt=1

M1G2 = SpikeMonitor(G2)
M2G2 = StateMonitor(G2, 'v', record=range(Vtt))
M3G2 = StateMonitor(G2, 'w', record=range(Vtt))
M4G2 = StateMonitor(G2, 'GsynE', record=range(Vtt))
M5G2 = StateMonitor(G2, 'GsynI', record=range(Vtt))
FRG2 = PopulationRateMonitor(G2)


print('--##Start simulation##--')
run(duration)
print('--##End simulation##--')

Cm     = Par['Cm']*farad
gl     = Par['gl']*siemens
El     = Par['El']*volt
Vt     = Par['Vt']*volt
Dt     = Par['Dt']*volt
tau_w  = Par['tau_w']*second
Is     = Par['Is']*amp  #2.50*nA #[0.0 for i in range(N2)]*nA
Ee     = Par['Ee']*volt
Ei     = Par['Ei']*volt
Tsyn   = Par['Tsyn']*second
a      = Par['a']*siemens
LwV=[]
Lww=[]

Vlim=[-85e-3,-30e-3] #
Glim=[-.1e-2, 0.25e-2] #[min(wsim)-2.e-3, max(wsim)+2.e-3] # [-0.5e-2, 5.5e-2] #[-2e-2, 2e-2] #[min(wsim)-2.e-3, max(wsim)+2.e-3]#

V=[(i+Vlim[0]*1000)*mV for i in range(int(1000*(Vlim[1]-Vlim[0]))+1)]# range(12000)]

for v in V:
    wV = (-gl*(v-El)+ gl*Dt*exp((v-Vt)/Dt) + Is) #-GsynE*(v-Ee)-GsynI*(v-Ei)
    ww = a*(v-El)
    LwV.append(wV)
    Lww.append(ww)

Vsim=array(M2G2.v[0])
wsim=array(M3G2.w[0]*1.0e7 )


LwV2=[i*1.0e7 for i in LwV]
Lww2=[i*1.0e7 for i in Lww]
ax3.plot(V/volt, LwV2,'b', linewidth=3)
ax3.plot(V/volt, Lww2,'k', linewidth=3)
ax3.set_ylim(Glim[0],Glim[1])
ax3.set_xlim(Vlim[0], Vlim[1])



RasG2 = np.array([M1G2.t/ms, M1G2.i])
#print(M1G2.t/ms)
LVs=[]
Lws=[]
lv=[]
lw=[]
TS=list(M1G2.t/ms)
VmNew=[]
for i in range(len(M1G2.t/ms)+1):
    if i==0:
        ai = 0
    else: 
        ai = TS[i-1]

    if i>=len(M1G2.t/ms):
        bi = TotTime
    else:
        bi = TS[i]        
  
    j=int(ai/DT)
    k=int(bi/DT)
    Vm=Vsim[j:k].copy()
    Vmem=Vm[-1]
    Vm[-1]=0.0*mV    
    VmNew.append(Vm)
    Vs=Vsim[(j+2):k].copy()
   # Vs[-1]=-40.0*mV
    ws=wsim[(j+2):k].copy()
    LVs.append(Vs)
    Lws.append(ws)

Vmf=[item for sublist in VmNew for item in sublist]
Vmf[-1]=Vmem

#plt.figure()

ax1.plot(Vmf, '0.3', c='r', linewidth=3)
ax2.plot(wsim, '0.3', c='r', linewidth=3)

#LVs[0]=list(LVs[0]).insert(0, initC[0])
#Lws[0]=list(Lws[0]).insert(0, initC[1])
#print(LVs[0])
LVs[0]=[initC[0]]+list(LVs[0])
Lws[0]=[initC[1]*1.0e7]+list(Lws[0])
#print(LVs, Lws)
for i in range(len(LVs)):
    ax3.plot(LVs[i],Lws[i], '0.3',c='r', linewidth=3)

    if i!=0:
        ax3.scatter(LVs[i][0],Lws[i][0], marker='x', c='r', s=100)
        ax3.plot([LVs[i][0],LVs[i-1][-1]],[Lws[i][0],Lws[i-1][-1]], c='#e38717', ls='--')
ax3.scatter(initC[0],initC[1]*1.0e7 , marker='o', c='r', s=200, zorder=10)



######################################################################################
######################################################################################
#Background: #16344d
#Arrows: #2b6899
#V-nullcline: #ffffff
#ga-nullcline: #95362c
#trajectory-continous: #e3b917
#trajectory-reset: #e38717
#reset-mark: #e38717


#ax1.set_ylim(Vlim[0], 10.0e-3)

ax1.set_ylabel('$V (mV)$')
ax2.set_ylabel('$w (pA)$')
ax2.set_xlabel('$time (ms)$')
ax3.set_xlabel('$V (mV)$')
ax3.set_ylabel('$w (pA)$')
#ax1.set_facecolor('#16344d')
#ax2.set_facecolor('#16344d')
#ax3.set_facecolor('#16344d')
ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*DT))
ax1.xaxis.set_major_formatter(ticks_x)
ax2.xaxis.set_major_formatter(ticks_x)

scaleV=1.e3
ticks_V = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*scaleV))
ax1.yaxis.set_major_formatter(ticks_V)
ax3.xaxis.set_major_formatter(ticks_V)

scalew=1.e2
ticks_w = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*scalew))
ax2.yaxis.set_major_formatter(ticks_w)
ax3.yaxis.set_major_formatter(ticks_w)

#ax1.set_yticks()

#-------------------------------------------------------
#-------------------------------------------------------

def dX_dt(X, t=0):
    return array([ (-gl*(X[0]-El)+ gl*Dt*exp((X[0]-Vt)/Dt)-X[1] + Is)/Cm,
                  ((a*(X[0]-El))-X[1])/tau_w ])
#-------------------------------------------------------

Cm     = Par['Cm']
gl     = Par['gl']
El     = Par['El']
Vt     = Par['Vt']
Dt     = Par['Dt']
tau_w  = Par['tau_w']
Is     = Par['Is']
Ee     = Par['Ee']
Ei     = Par['Ei']
Tsyn   = Par['Tsyn']
a      = Par['a']

# define a grid and compute direction at each point

nb_points   = 20

x = linspace(Vlim[0], Vlim[1], nb_points)
y = linspace(Glim[0]*1.0e-7,Glim[1]*1.0e-7, nb_points)

X1 , Y1  = meshgrid(x, y)                       # create a grid
DX1, DY1 = dX_dt([X1, Y1])                      # compute growth rate on the gridt
Y1 = Y1*1.0e7
DY1 = DY1*1.0e7
M = (hypot(DX1, DY1))                           # Norm of the growth rate 
M[ M == 0] = 1.                                 # Avoid zero division errors 
DX1 /= M                                        # Normalize each arrows
DY1 /= M

#print(M2)
#-------------------------------------------------------
# Drow direction fields, using matplotlib 's quiver function
# I choose to plot normalized arrows and to use colors to give information on
# the growth speed
#p.title('Trajectories and direction fields')

Q = ax3.quiver(X1, Y1, DX1, DY1, pivot='mid',color='#16344d', alpha=0.5)#, units="xy", scale=5e3) #, linewidths=M2, cmap=p.cm.jet)


#plt.savefig(filename+'.svg', format='svg', dpi=1000)
#plt.savefig(filename+'.eps', format='eps', dpi=1000)

plt.show()
