from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TH1D,TH1, TFile, TLorentzVector, TMath, TTree, TVector3, TEfficiency, TGraphAsymmErrors
from math import *
from optparse import OptionParser
from array import array
import os
import fnmatch
import uproot
import mplhep as hep
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

#histofiles = ['ntuple_photons_BIBjets_DR1_HF1_0_50.root', 'ntuple_photons_BIBjets_noDR_HF1_0_50.root']#,'ntuple_photons_BIBjets_DR5_HF5_0_50.root','ntuple_photons_BIBjets_DR5_HF4_0_50.root']#, 
        #'ntuple_photons_BIBjets_DR1_0_50.root', 'ntuple_photons_BIBjets_DR05_0_50.root']
histofiles = ['ntuple_photons_pTmatch_12-09_0_50.root', 'ntuple_photons_pTmatch_12-09_50_250.root', 'ntuple_photons_pTmatch_12-09_250_1000_BATCH1.root', 'ntuple_photons_pTmatch_12-09_250_1000_BATCH2.root']
histofiles_B = ['ntuple_photons_BIBpTmatch_0_50.root', 'ntuple_photons_BIBpTmatch_50_250.root', 'ntuple_photons_BIBpTmatch_250_1000.root']
cutoffs = [ 0.1]# 0.1, 0.05]
filenum = 0
fig, ax = plt.subplots()


up_err_arr = []
low_err_arr = []
up_err_arr_B = []
low_err_arr_B = []
for filename in histofiles:
    file = uproot.open(filename)
    pass_histo = file['E_pass;1']
    all_histo = file['E_all;1']

    pass_histo_np = pass_histo.to_numpy()
    all_histo_np = all_histo.to_numpy()
    pass_histo_arr = pass_histo_np[0]
    all_histo_arr = all_histo_np[0]
    if filenum == 0:
        pass_slice_0 = pass_histo_arr[0:]
        all_slice_0 = all_histo_arr[0:]
    elif filenum == 1:
        pass_slice_1 = pass_histo_arr[0:]
        all_slice_1 = all_histo_arr[0:]
    elif filenum == 2:
        pass_slice_2 = pass_histo_arr[0:]
        all_slice_2 = all_histo_arr[0:]
    elif filenum == 3:
        pass_slice_3 = pass_histo_arr[0:]
        all_slice_3 = all_histo_arr[0:]
    all_histo_arr = all_histo_np[0]
    #efficiency_arr = pass_slice/all_slice
    E_arr = pass_histo_np[1]
    E_slice = E_arr[1:]
    filenum = filenum+1


filenum = 0
for filename in histofiles_B:
    file = uproot.open(filename)
    pass_histo = file['E_pass;1']
    all_histo = file['E_all;1']

    pass_histo_np = pass_histo.to_numpy()
    all_histo_np = all_histo.to_numpy()
    pass_histo_arr = pass_histo_np[0]
    all_histo_arr = all_histo_np[0]
    if filenum == 0:
        pass_slice_0_B = pass_histo_arr[0:]
        all_slice_0_B = all_histo_arr[0:]
    elif filenum == 1:
        pass_slice_1_B = pass_histo_arr[0:]
        all_slice_1_B = all_histo_arr[0:]
    elif filenum == 2:
        pass_slice_2_B = pass_histo_arr[0:]
        all_slice_2_B = all_histo_arr[0:]
    all_histo_arr = all_histo_np[0]
    #efficiency_arr = pass_slice/all_slice


    #upErrors=array('d')
    #lowErrors=array('d')
    #binEff=array('d')
    #photonEff = TGraphAsymmErrors(pass_histo, all_histo)
    E_arr_B = pass_histo_np[1]
    E_slice_B = E_arr[1:]
    #plot using MAIA conventions
    #color = plt.cm.plasma(filenum/len(cutoffs))
    filenum = filenum + 1
pass_slice = pass_slice_0+pass_slice_1+pass_slice_2+pass_slice_3
all_slice = all_slice_0+all_slice_1+all_slice_2+all_slice_3
efficiency_arr = pass_slice/all_slice

pass_slice_B = pass_slice_0_B + pass_slice_1_B + pass_slice_2_B
all_slice_B = all_slice_0_B + all_slice_1_B + all_slice_2_B
efficiency_arr_B = pass_slice_B / all_slice_B

for i in range(len(efficiency_arr)-1):
    N = pass_slice[i]
    k = all_slice[i]
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

for i in range(len(efficiency_arr_B)-1):
    N = pass_slice_B[i]
    k = all_slice_B[i]
    lower_error = efficiency_arr_B[i] - TEfficiency.ClopperPearson(k, N,  0.683,  False)
    upper_error = TEfficiency.ClopperPearson(k, N,  0.683,  True) - efficiency_arr_B[i]
    up_err_arr_B.append(upper_error)
    low_err_arr_B.append(lower_error)
asymm_errs = [low_err_arr[6:], up_err_arr[6:]]
asymm_errs_B = [low_err_arr_B[6:], up_err_arr_B[6:]]
#print(len(up_err_arr))
E_errs = []
E_arr = []
for i in range(len(E_slice)-1):
    E_arr.append((E_slice[i+1]+E_slice[i])/2)
    E_errs.append((E_slice[i+1]-E_slice[i])/2)

print(len(E_arr[6:]))
print(len(efficiency_arr[6:len(efficiency_arr)-1]))
print(len(E_errs[6:]))
print(len(asymm_errs[0]))


plt.errorbar(E_arr[6:], efficiency_arr[6:len(efficiency_arr)-1], xerr = E_errs[6:], yerr = asymm_errs,fmt=' ',
        color='blue', label = "Without BIB Overlay")#label="BIB, jets, HF < 0.1, dR <= "+str(cutoffs[filenum]))
plt.errorbar(E_arr[6:], efficiency_arr_B[6:len(efficiency_arr)-1], xerr = E_errs[6:], yerr = asymm_errs_B, fmt = ' ', color = 'red', label="With BIB Overlay") 
plt.axhline(1.,ls='--', lw=0.5)
ax.set_ylim(0.,1.4)
ax.set_xlim(E_slice[6],E_slice[len(E_slice)-1])
ax.set_xlabel("True Photon E [GeV]", loc='right')
ax.set_ylabel("Photon Reconstruction Efficiency", loc='top')

hep.cms.label(exp = "Muon Collider", data = False, label = "no BIB (v0.4 Lattice)", 
       rlabel='$MAIA$ Detector Concept', loc=2, italic=(1,0,0), pad=(0.0))
plt.gcf().text(0.16,0.74,"$\sqrt{s}=10$ TeV")
plt.gcf().text(0.16, 0.7, "30 GeV < True Photon E < 1 TeV")
#plt.text(670,0.135, "MAIA", fontstyle = 'italic')
#plt.text(755, 0.135, "Detector Concept")
#plt.text(670, 0.13, "BIB Overlay in [-0.5, 15] ns", fontsize = 'small')
ax.legend(frameon=False)

ax.tick_params(bottom=True, top=True, left=True, right=True, which='both', direction='in')
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
ax.tick_params(axis= 'y',which='major',length = 10)
ax.tick_params(axis='y',which='minor', length=5)

ax.xaxis.set_major_locator(MultipleLocator(200))
ax.xaxis.set_major_formatter('{x:.0f}')
ax.xaxis.set_minor_locator(MultipleLocator(50))

ax.yaxis.set_major_locator(MultipleLocator(0.2))
ax.yaxis.set_major_formatter('{x:.1f}')
ax.yaxis.set_minor_locator(MultipleLocator(0.04))
plt.savefig("gammaEff_PFOs_maia_E.pdf")
plt.close()
print("gammaEff_PFOs_maia_E.pdf created")
