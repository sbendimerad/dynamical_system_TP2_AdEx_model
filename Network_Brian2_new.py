from brian2 import *


# Initialize lists to store results for different input rates
input_rates = [5, 20, 50] ## Change the name of the parameter you want to test and the values to be tested
mean_firing_rates_RS = []
mean_firing_rates_FS = []
CVs_RS = []
CVs_FS = []
correlations_RS = []
correlations_FS = []

for i in input_rates:
    start_scope()

    DT=0.1 # time step
    defaultclock.dt = DT*ms
    N_inh = 2000 # number of inhibitory neurons
    N_exc = 8000 # number of excitatory neurons 

    TotTime=500 #Simulation duration (ms)
    duration = TotTime*ms



    # Equations ----------------------------------------------------------------------------------
    eqs='''
    dv/dt = (-GsynE*(v-Ee)-GsynI*(v-Ei)-gl*(v-El)+ gl*Dt*exp((v-Vt)/Dt)-w + Is)/Cm : volt (unless refractory)
    dw/dt = (a*(v-El)-w)/tau_w:ampere
    dGsynI/dt = -GsynI/Tsyn : siemens
    dGsynE/dt = -GsynE/Tsyn : siemens
    Is:ampere
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

    # Populations----------------------------------------------------------------------------------


    # Population 1 - Fast Spiking

    # Population 1 - FS - Inhibitory

    G_inh = NeuronGroup(N_inh, eqs, threshold='v > -47.5*mV', reset='v = -65*mV', refractory='5*ms', method='heun')
    #init:
    G_inh.v = -65*mV
    G_inh.w = 0.0*pA
    G_inh.GsynI=0.0*nS
    G_inh.GsynE=0.0*nS
    #parameters
    G_inh.Cm = 200.*pF
    G_inh.gl = 10.*nS
    G_inh.El = -65.*mV
    G_inh.Vt = -50.*mV
    G_inh.Dt = 0.5*mV
    G_inh.tau_w = 1.0*ms
    G_inh.a = 0.0*nS
    G_inh.Is = 0.0 #[0.0 for i in range(N1)]*nA

    G_inh.Ee=0.*mV
    G_inh.Ei=-80.*mV
    G_inh.Tsyn=5.*ms

    # Population 2 - RS - Excitatory 
    b_exc = 60.*pA

    G_exc = NeuronGroup(N_exc, eqs, threshold='v > -40.0*mV', reset='v = -65*mV; w += b_exc', refractory='5*ms',  method='heun')
    G_exc.v = -65.*mV
    G_exc.w = 0.0*pA
    G_exc.GsynI=0.0*nS
    G_exc.GsynE=0.0*nS
    G_exc.Cm = 200.*pF
    G_exc.gl = 10.*nS
    G_exc.El = -70.*mV
    G_exc.Vt = -50.*mV
    G_exc.Dt = 2.*mV
    G_exc.tau_w = 1000.*ms
    G_exc.a = 0.*nS
    G_exc.Is = 0.0*nA #2.50*nA #[0.0 for i in range(N2)]*nA
    #G_exc.Is[0]=0.*nA

    G_exc.Ee=0.*mV
    G_exc.Ei=-80.*mV
    G_exc.Tsyn=5.*ms


    # external drive--------------------------------------------------------------------------
    # update_1
    ext_inp=i
    # ext_inp=4.0
    P_ed=PoissonGroup(8000, rates=ext_inp*Hz)

    # Network-----------------------------------------------------------------------------

    # quantal increment in synaptic conductances:
    # update_2
    Qi=8.0*nS
    Qe=4.0*nS
    # Qi=5.0*nS
    # Qe=1.5*nS


    # probability of connection
    prbC= 0.05

    #create synapses
    S_12 = Synapses(G_inh, G_exc, on_pre='GsynI_post+=Qi') 
    S_12.connect('i!=j', p=prbC)

    S_11 = Synapses(G_inh, G_inh, on_pre='GsynI_post+=Qi')
    S_11.connect('i!=j',p=prbC)

    S_21 = Synapses(G_exc, G_inh, on_pre='GsynE_post+=Qe')
    S_21.connect('i!=j',p=prbC)

    S_22 = Synapses(G_exc, G_exc, on_pre='GsynE_post+=Qe')
    S_22.connect('i!=j', p=prbC)

    S_ed_in = Synapses(P_ed, G_inh, on_pre='GsynE_post+=Qe')
    S_ed_in.connect(p=prbC)

    S_ed_ex = Synapses(P_ed, G_exc, on_pre='GsynE_post+=Qe')
    S_ed_ex.connect(p=prbC)


    # Recording tools -------------------------------------------------------------------------------

    M1G_inh = SpikeMonitor(G_inh)
    FRG_inh = PopulationRateMonitor(G_inh)
    M1G_exc = SpikeMonitor(G_exc)
    FRG_exc = PopulationRateMonitor(G_exc)




    # Useful trick to record global variables ------------------------------------------------------

    Gw_inh = NeuronGroup(1, 'Wtot : ampere', method='rk4')
    Gw_exc = NeuronGroup(1, 'Wtot : ampere', method='rk4')

    SwInh1=Synapses(G_inh, Gw_inh, 'Wtot_post = w_pre : ampere (summed)')
    SwInh1.connect(p=1)
    SwExc1=Synapses(G_exc, Gw_exc, 'Wtot_post = w_pre : ampere (summed)')
    SwExc1.connect(p=1)

    MWinh = StateMonitor(Gw_inh, 'Wtot', record=0)
    MWexc = StateMonitor(Gw_exc, 'Wtot', record=0)



    GV_inh = NeuronGroup(1, 'Vtot : volt', method='rk4')
    GV_exc = NeuronGroup(1, 'Vtot : volt', method='rk4')

    SvInh1=Synapses(G_inh, GV_inh, 'Vtot_post = v_pre : volt (summed)')
    SvInh1.connect(p=1)
    SvExc1=Synapses(G_exc, GV_exc, 'Vtot_post = v_pre : volt (summed)')
    SvExc1.connect(p=1)

    MVinh = StateMonitor(GV_inh, 'Vtot', record=0)
    MVexc = StateMonitor(GV_exc, 'Vtot', record=0)


    # Run simulation -------------------------------------------------------------------------------

    print('--##Start simulation##--')
    run(duration)
    print('--##End simulation##--')


    # Plots -------------------------------------------------------------------------------


    # prepare raster plot
    RasG_inh = array([M1G_inh.t/ms, [i+N_exc for i in M1G_inh.i]])
    RasG_exc = array([M1G_exc.t/ms, M1G_exc.i])



    ''' The following function is very important as it defines the time windows through which the rates will be calculated. 
    In other words, out of an ensemble of discrete spiking events, we set bins to count them, and form a mean population firing rate that is time dependent. The size of the bins most likely influence the results we can obtain, and this is one of the questions you can ask'''


    def bin_array(array, BIN, time_array):
        N0 = int(BIN/(time_array[1]-time_array[0]))
        N1 = int((time_array[-1]-time_array[0])/BIN)
        return array[:N0*N1].reshape((N1,N0)).mean(axis=1)

    BIN=5 ## Size of the time windows in ms
    time_array = arange(int(TotTime/DT))*DT



    LfrG_exc=array(FRG_exc.rate/Hz)
    TimBinned,popRateG_exc=bin_array(time_array, BIN, time_array),bin_array(LfrG_exc, BIN, time_array)

    LfrG_inh=array(FRG_inh.rate/Hz)
    TimBinned,popRateG_inh=bin_array(time_array, BIN, time_array),bin_array(LfrG_inh, BIN, time_array)

    # Compute mean firing rate for RS (excitatory) and FS (inhibitory) neurons
    mean_firing_rate_RS = np.mean(LfrG_exc)
    mean_firing_rate_FS = np.mean(LfrG_inh)



    # create the figure

    fig=figure(figsize=(8,12))
    ax1=fig.add_subplot(211)
    ax2=fig.add_subplot(212)


    ax1.plot(RasG_inh[0], RasG_inh[1], ',r')
    ax1.plot(RasG_exc[0], RasG_exc[1], ',g')

    ax1.set_xlabel('Time (ms)')
    ax1.set_ylabel('Neuron index')

    ax2.plot(TimBinned,popRateG_inh, 'r')
    ax2.plot(TimBinned,popRateG_exc, 'g')

    ax2.set_xlabel('Time (ms)')
    ax2.set_ylabel('Firing Rate (Hz)')

    name_fig='Fig_2pop_Reg_ext_'+str(ext_inp)+'.png'
    plt.savefig(name_fig)

    name_rates1='FR_2pop_Reg_ext1'+str(ext_inp)+'.npy'
    name_rates1='FR_2pop_Reg_ext2'+str(ext_inp)+'.npy'

    np.save(name_rates1,[TimBinned, popRateG_inh,popRateG_exc])
    np.save(name_rates1,[LfrG_inh,LfrG_exc])


    np.save('Wtot.npy',[np.array(MWinh.Wtot[0]/mamp),np.array(MWexc.Wtot[0]/mamp)])

    np.save('Vtot.npy',[np.array(MVinh.Vtot[0]/mV),np.array(MVexc.Vtot[0]/mV)])

    fig.tight_layout()

        # Calculate mean firing rates
    mean_firing_rates_RS.append(np.mean(FRG_exc.rate/Hz))
    mean_firing_rates_FS.append(np.mean(FRG_inh.rate/Hz))
    
    # Calculate synchronization (CV)
    ISI_RS = np.diff(M1G_exc.t/ms)
    ISI_FS = np.diff(M1G_inh.t/ms)
    CVs_RS.append(np.std(ISI_RS) / np.mean(ISI_RS) if len(ISI_RS) > 1 else 0)
    CVs_FS.append(np.std(ISI_FS) / np.mean(ISI_FS) if len(ISI_FS) > 1 else 0)
    
    # Calculate correlation
    individual_rates_RS = np.bincount(M1G_exc.i, minlength=N_exc) / duration
    individual_rates_FS = np.bincount(M1G_inh.i, minlength=N_inh) / duration
    global_rate_RS = np.mean(individual_rates_RS)
    global_rate_FS = np.mean(individual_rates_FS)
    correlations_RS.append(np.corrcoef(individual_rates_RS, [global_rate_RS]*N_exc)[0, 1])
    correlations_FS.append(np.corrcoef(individual_rates_FS, [global_rate_FS]*N_inh)[0, 1])

    # show()

# Print summarized results
print("Input rates: ", input_rates)
print("Mean firing rates (RS): ", mean_firing_rates_RS)
print("Mean firing rates (FS): ", mean_firing_rates_FS)
print("Synchronization (CV) - RS: ", CVs_RS)
print("Synchronization (CV) - FS: ", CVs_FS)
print("Correlation (RS): ", correlations_RS)
print("Correlation (FS): ", correlations_FS)
