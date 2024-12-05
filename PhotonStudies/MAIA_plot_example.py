import numpy as np
from array import array
import matplotlib.pyplot as plt
import mplhep as hep
import uproot
import ROOT
from ROOT import TEfficiency
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

#grab whatever histos you need using uproot, etc
fig, ax = plt.subplots()
plt.errorbar(#data here)
plt.ylim(0,.1) #replace w appropriate range
plt.xlim(10.,50.)
#replace with appropriate titles
plt.xlabel("True Photon Energy [GeV]", loc='right')
plt.ylabel("Photon Energy Resolution $\sigma((E_{reco}-E_{true})/E_{true}$)", loc='top')
#MAIA formatting (including correct italicization)
hep.cms.label(exp = "Muon Collider", data = False, label = ' ', 
       rlabel='$MAIA$ Detector Concept', loc=2, italic=(1,0,0), pad=(0.0))
#place 10TeV label below 
plt.gcf().text(0.16,0.74,"($\sqrt{s}=10$ TeV)")

#mimic ROOT-style tick marks
ax.tick_params(bottom=True, top=True, left=True, right=True, which='both', direction='in')
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
ax.tick_params(axis= 'y',which='major',length = 10)
ax.tick_params(axis='y',which='minor', length=5)

#adjust to appropriate range
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_major_formatter('{x:.0f}')
ax.xaxis.set_minor_locator(MultipleLocator(1))

ax.yaxis.set_major_locator(MultipleLocator(0.02))
ax.yaxis.set_major_formatter('{x:.2f}')
ax.yaxis.set_minor_locator(MultipleLocator(0.005))

#may be helpful for resolution curve, not necessary
plt.grid(axis = 'y', linestyle = ':')
#if multiple curves
plt.legend(loc='upper right',frameon=False)



####### Specifically for EFFICIENCY plots ###########
###### To get the same type of asymmetrical errors as TEfficiency #####

pass_slice = [] #fill with the pass histo
all_slice = [] #fill with the all histo
efficiency_arr = [] #the pass_slice / all_slice
for i in range(len(efficiency_arr)-1):
    N = pass_slice[i]
    k = all_slice[i]
    #use the Clopper Pearson error calculation included in TEfficiency
    lower_error = efficiency_arr[i] - TEfficiency.ClopperPearson(k, N,  0.683,  False)
    upper_error = TEfficiency.ClopperPearson(k, N,  0.683,  True) - efficiency_arr[i] 
    up_err_arr.append(upper_error)
    low_err_arr.append(lower_error)
asymm_errs = [low_err_arr, up_err_arr]
E_errs = []
E_arr = []
for i in range(len(E_slice)-1):
    E_arr.append((E_slice[i+1]+E_slice[i])/2)
    E_errs.append((E_slice[i+1]-E_slice[i])/2)


