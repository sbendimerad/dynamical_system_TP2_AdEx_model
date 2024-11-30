from brian2 import *
import random

#########################
#Simulations conditions #-------------------------------------------------------------------------
#########################
start_scope()

DT=0.1 # time step
defaultclock.dt = DT*ms
Nb1 = 0 # number of inhibitory neudel
Nb2 = 2 # number of excitatory neudel 

TotTime=4000 #Simulation duration (ms)
duration = TotTime*ms


###############
# Equations ###----------------------------------------------------------------------------------
###############
eqs='''
#Neudel:
dv/dt = (-GsynE*(v-Ee)-GsynI*(v-Ei)-gl*(v-El)+ gl*Dt*exp((v-Vt)/Dt)-w + Is+Ig)/Cm : volt (unless refractory)
dw/dt = (a*(v-El)-w)/tau_w:ampere
#model of synapses:
dGsynI/dt = -GsynI/Tsyn : siemens
dGsynE/dt = -GsynE/Tsyn : siemens
#parameters:
Is:ampere
Ig:ampere
Cm:farad
gl:siemens
El:volt
a:siemens
tau_w:second
Dt:volt
Vt:volt
Ee:volt
Ei:volt
Tsyn:second
'''
##################
# Populations ####-----------------------------------------------------------------------------
##################
if Nb1!=0:
    # Inhibitory

    G_inh = NeuronGroup(Nb1, eqs, threshold='v > -40.*mV', reset='v = -52*mV; w += b1', refractory='5*ms', method='heun')
    #init variables:
    for k in range(len(G_inh.v)):
        G_inh.v[k] = random.uniform(-55, -65)*mV
    
    G_inh.w = 0.0*pA
    G_inh.GsynI=0.0*nS
    G_inh.GsynE=0.0*nS
    #parameter values:
    b1 = 0.0*pA
    G_inh.Cm = 200.*pF
    G_inh.gl = 10.*nS
    G_inh.El = -65.*mV
    G_inh.Vt = -55.*mV
    G_inh.Dt = 5.0*mV
    G_inh.tau_w = 500.0*ms
    G_inh.a = 0.0*nS
    G_inh.Is = .06*nA 

    G_inh.Ee=0.*mV
    G_inh.Ei=-80.*mV
    G_inh.Tsyn=5.*ms
    # Recording tools

    M2G1 = StateMonitor(G_inh, 'v', record=range(Nb1))
    M3G1 = StateMonitor(G_inh, 'w', record=range(Nb1))
    M4G1 = StateMonitor(G_inh, 'GsynE', record=range(Nb1))
    M5G1 = StateMonitor(G_inh, 'GsynI', record=range(Nb1))
    M1G1 = SpikeMonitor(G_inh)
    FRG1 = PopulationRateMonitor(G_inh)

if Nb2!=0:
    # Excitatory

    G_exc = NeuronGroup(Nb2, eqs, threshold='v > -40.0*mV', reset='v = -52*mV; w += b2', refractory='5*ms',  method='heun')

    #init variables:
    for k in range(len(G_exc.v)):
        G_exc.v[k] = random.uniform(-55, -65)*mV

    G_exc.w = 0.0*pA
    G_exc.GsynI=0.0*nS
    G_exc.GsynE=0.0*nS

    #parameter values
    b2 = 10.*pA
    G_exc.Cm = 200.*pF
    G_exc.gl = 10.*nS
    G_exc.El = -65.*mV
    G_exc.Vt = -55.*mV
    G_exc.Dt = 5.*mV
    G_exc.tau_w = 500.*ms
    G_exc.a = 2.*nS
    G_exc.Is = .120*nA 


    G_exc.Ee=0.*mV
    G_exc.Ei=-80.*mV
    G_exc.Tsyn=5.*ms

    # Recording tools
    M2G2 = StateMonitor(G_exc, 'v', record=range(Nb2))
    M3G2 = StateMonitor(G_exc, 'w', record=range(Nb2))
    M4G2 = StateMonitor(G_exc, 'GsynE', record=range(Nb2))
    M5G2 = StateMonitor(G_exc, 'GsynI', record=range(Nb2))
    M1G2 = SpikeMonitor(G_exc)
    FRG2 = PopulationRateMonitor(G_exc)



##################
# Network ########-----------------------------------------------------------------------------
##################

# quantal increment in synaptic conductances:
Qi=1.5*nS
Qe=.5*nS

#create interaction through model synapses

if Nb1!=0:
    S_11 = Synapses(G_inh, G_inh, on_pre='GsynI_post+=Qi')
    S_11.connect(p=1)


if Nb2!=0:
    S_22 = Synapses(G_exc, G_exc, on_pre='GsynE_post+=Qe')
    S_22.connect(p=1)


if Nb1!=0 and Nb2!=0:
    S_21 = Synapses(G_exc, G_inh, on_pre='GsynE_post+=Qe')
    S_21.connect(p=1)
    S_12 = Synapses(G_inh, G_exc, on_pre='GsynI_post+=Qi') 
    S_12.connect(p=1)

##################
# Run simulation # -------------------------------------------------------------------------------
##################
print('--##Start simulation##--')
run(duration)
print('--##End simulation##--')

#################
# Plots #########-----------------------------------------------------------------------------
#################

# Prepare data 
def bin_array(array, BIN, time_array):
    N0 = int(BIN/(time_array[1]-time_array[0]))
    N1 = int((time_array[-1]-time_array[0])/BIN)
    return array[:N0*N1].reshape((N1,N0)).mean(axis=1)

BIN=5
time_array = arange(int(TotTime/DT))*DT


LVG1=[]
LwG1=[]
LVG2=[]
LwG2=[]

LgseG1=[]
LgsiG1=[]
LgseG2=[]
LgsiG2=[]

for a in range(Nb1):
    LVG1.append(array(M2G1[a].v/mV))
    LwG1.append(array(M3G1[a].w/mamp))
    LgseG1.append(array(M4G1[a].GsynE/nS))
    LgsiG1.append(array(M5G1[a].GsynI/nS))

for a in range(Nb2):
    LVG2.append(array(M2G2[a].v/mV))
    LwG2.append(array(M3G2[a].w/mamp))
    LgseG2.append(array(M4G2[a].GsynE/nS))
    LgsiG2.append(array(M5G2[a].GsynI/nS))



if Nb1!=0:
    LfrG1=np.array(FRG1.rate/Hz)
    TimBinned,popRateG_inh=bin_array(time_array, BIN, time_array),bin_array(LfrG1, BIN, time_array)
    Lt1G1=array(M2G1.t/ms)
    Lt2G1=array(M3G1.t/ms)
    RasG_inh = array([M1G1.t/ms, [i+Nb2 for i in M1G1.i]])

if Nb2!=0:
    LfrG2=np.array(FRG2.rate/Hz)
    TimBinned,popRateG_exc=bin_array(time_array, BIN, time_array),bin_array(LfrG2, BIN, time_array)
    Lt1G2=array(M2G2.t/ms)
    Lt2G2=array(M3G2.t/ms)
    RasG_exc = array([M1G2.t/ms, M1G2.i])


#create the figure
fig=plt.figure(figsize=(12,8))
ax1=fig.add_subplot(221)
ax2=fig.add_subplot(222)
ax3=fig.add_subplot(223)
ax4=fig.add_subplot(224)
for a in range(len(LVG1)):
    ax1.plot(Lt1G1, LVG1[a], color ='r')
    ax1.plot(Lt1G1, LVG1[a], ls=(0, (2,6))) 
    ax3.plot(Lt2G1, LwG1[a],color ='r')
    ax3.plot(Lt2G1, LwG1[a],ls=(0, (2,6))) 

for a in range(len(LVG2)):
    ax1.plot(Lt1G2, LVG2[a],color ='g')
    ax1.plot(Lt1G2, LVG2[a], ls=(0, (2,6))) 
    ax3.plot(Lt2G2, LwG2[a],color ='g')
    ax3.plot(Lt2G2, LwG2[a], ls=(0, (2,6)))


if Nb1!=0:
    ax2.plot(RasG_inh[0], RasG_inh[1], '.r')
    ax4.plot(TimBinned,popRateG_inh, 'r')


if Nb2!=0:
    ax2.plot(RasG_exc[0], RasG_exc[1], '.g')
    ax4.plot(TimBinned,popRateG_exc, 'g')

ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Vm (mV)')

ax2.set_xlabel('Time (ms)')
ax2.set_ylabel('Neudel index')

ax3.set_xlabel('Time (ms)')
ax3.set_ylabel('w (nA)')

ax4.set_xlabel('Time (ms)')
ax4.set_ylabel('Firing Rate')


if Nb1!=0 and Nb2!=0:
    fig2=figure(figsize=(8,12))
    ax21=fig2.add_subplot(211)
    ax22=fig2.add_subplot(212)
    ax21.set_title('N_inh #0 vs N_Exc #0')
    ax21.plot(LVG1[0], LVG2[0], 'b')
    ax21.set_xlabel('Vm inh')
    ax21.set_ylabel('Vm exc')
    ax22.set_title('Vm vs w for neudels #0')
    ax22.plot(LVG1[0], LwG1[0], 'r')
    ax22.plot(LVG2[0], LwG2[0], 'g')
    ax22.set_xlabel('Vm')
    ax22.set_ylabel('w')

if Nb1!=0 and Nb2==0:
    fig2=figure(figsize=(8,12))
    ax21=fig2.add_subplot(211)
    ax22=fig2.add_subplot(212)
    ax21.set_title('N_inh #0 vs N_inh #1')
    ax21.plot(LVG1[0], LVG1[1], 'r')
    ax21.set_xlabel('Vm inh #0')
    ax21.set_ylabel('Vm inh #1')
    ax22.set_title('Vm vs w for neudels #0')
    ax22.plot(LVG1[0], LwG1[0], 'r')
    #ax22.plot(LVG1[1], LwG1[1], 'r')
    ax22.set_xlabel('Vm inh')
    ax22.set_ylabel('w inh')

if Nb1==0 and Nb2!=0:
    fig2=figure(figsize=(8,12))
    ax21=fig2.add_subplot(211)
    ax22=fig2.add_subplot(212)
    ax21.set_title('N_exc #0 vs N_exc #1')
    ax21.plot(LVG2[0], LVG2[1], 'g')
    ax21.set_xlabel('Vm exc #0')
    ax21.set_ylabel('Vm exc #1')
    ax22.set_title('Vm vs w for neudels #0')
    ax22.plot(LwG2[0], LVG2[0], 'g')
    #ax22.plot(LwG2[1], LVG2[1], 'g')
    ax22.set_xlabel('w inh')
    ax22.set_ylabel('Vm inh')


fig.tight_layout()

show()
