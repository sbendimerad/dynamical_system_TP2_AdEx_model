# Model from  Di Volo et al. Neural Comp. 2019


import numpy as np
import matplotlib.pylab as plt
import mpl_toolkits.mplot3d.axes3d as axes3d
import sys
import matplotlib as mpl
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as pl
from mpmath import *
import plotly.graph_objects as go
from typing import Tuple, Iterable
from plotly.subplots import make_subplots
import plotly.express as px
import pandas as pd


def TF(P,fexc,finh,adapt):
    

    fe = fexc*(1.-gei)*pconnec*Ntot;
    fi = finh*gei*pconnec*Ntot;
    
    muGi = Qi*Ti*fi;
    muGe = Qe*Te*fe;
    muG = Gl+muGe+muGi;
    muV = (muGe*Ee+muGi*Ei+Gl*El-adapt)/muG;
    
    
    muGn = muG/Gl;
    Tm = Cm/muG;
    
    Ue =  Qe/muG*(Ee-muV);
    Ui = Qi/muG*(Ei-muV);
    
    sV = sqrt(fe*(Ue*Te)*(Ue*Te)/2./(Te+Tm)+fi*(Ui*Ti)*(Ui*Ti)/2./(Ti+Tm));
    
    
    fe= fe+1e-9;
    fi=fi+1e-9;
    Tv = ( fe*(Ue*Te)*(Ue*Te) + fi*(Qi*Ui)*(Qi*Ui)) /( fe*(Ue*Te)*(Ue*Te)/(Te+Tm) + fi*(Qi*Ui)*(Qi*Ui)/(Ti+Tm) );
    TvN = Tv*Gl/Cm;
    
    muV0=-60e-3;
    DmuV0 = 10e-3;
    sV0 =4e-3;
    DsV0= 6e-3;
    TvN0=0.5;
    DTvN0 = 1.;
    vthr=P[0]+P[1]*(muV-muV0)/DmuV0+P[2]*(sV-sV0)/DsV0+P[3]*(TvN-TvN0)/DTvN0+P[5]*((muV-muV0)/DmuV0)*((muV-muV0)/DmuV0)+P[6]*((sV-sV0)/DsV0)*((sV-sV0)/DsV0)+P[7]*((TvN-TvN0)/DTvN0)*((TvN-TvN0)/DTvN0)+P[8]*(muV-muV0)/DmuV0*(sV-sV0)/DsV0+P[9]*(muV-muV0)/DmuV0*(TvN-TvN0)/DTvN0+P[10]*(sV-sV0)/DsV0*(TvN-TvN0)/DTvN0;
    
    
    
    
    
    frout=.5/TvN*Gl/Cm*(1-erf((vthr-muV)/sqrt(2)/sV));
    
    return frout;







Gl=10*1.e-9;
Cm=200*1.e-12;
El=-65*1.e-3;



Qe=1.5*1.e-9;
Qi=5.*1.e-9;
Ti=5*1.e-3;
Te=5*1.e-3;
Ee=0;
Ei=-80*1.e-3;
pconnec=0.05;
gei=0.2;
Ntot=10000;
T=0.015


bRS=60*1.e-12;
twRS=.5;



Qee=Qe
Qei=Qe
Qii=Qi
Qie=Qi


external_input=1.5



PRS=np.load('RS-cell0_CONFIG1_fit.npy')

PFS=np.load('FS-cell_CONFIG1_fit.npy')




tfinal=10.
dt=0.001

t = np.linspace(0, tfinal, int(tfinal/dt))

f = open('time.txt', 'wb')

fecont=2;
ficont=10;
w=fecont*bRS*twRS


LSw=[]
LSfe=[]
LSfi=[]

for i in range(len(t)):
    
    fecontold=fecont
    fecont+=dt/T*(TF(PRS, fecont+external_input,ficont,w)-fecont)
    w+=dt*( -w/twRS+(bRS)*fecontold)
    ficont+=dt/T*(TF(PFS,fecontold+external_input,ficont,0.)-ficont)

    LSfe.append(float(fecont))
    LSfi.append(float(ficont))
    LSw.append(float(w))
    #f.write("%10.3E %10.3E %10.3E\n" % (t[i],fecont,ficont))

plt.plot(t, LSfe)
plt.plot(t, LSfi)
fig=plt.figure()
#plt.plot(LSfe, LSfi)
ax = fig.add_subplot(1, 1, 1, projection = '3d')
ax.plot(LSfe, LSfi, LSw)
plt.show()
#f.close()
