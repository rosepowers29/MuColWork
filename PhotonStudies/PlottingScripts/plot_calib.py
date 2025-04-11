import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from optparse import OptionParser
from array import array
import os
import logging
import itertools
import fnmatch
import csv
import pandas as pd
import mplhep as hep
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

EBins = array('d', (10.,50., 100., 150., 200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800.))#, 850., 900., 950., 1000.))
EBins_lowE = array('d', (0., 50., 100., 150., 200., 250.))#, 30., 40., 50.))
#ThetaBins = array('d', (20.*TMath.Pi()/180, 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180., 90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., 160*TMath.Pi()/180))
#ThetaBins = np.linspace(0.175,2.96,30)
endcap_bins = np.linspace(0.225, 0.625, 5)
endcap_bins_2 = np.linspace(2.517, 2.91, 5)
tb_bins = np.linspace(0.625,2.517,20)
ebins_1 = np.array([0.175])
ebins_2 = np.array([2.96])
#ThetaBins_neutrons = np.concatenate((ebins_1, tb_bins, ebins_2))
ThetaBins = np.concatenate((endcap_bins[:2], tb_bins, endcap_bins_2[1:]))


with open('responseMap_photons_newregions_noBIB.csv', 'r') as csvToRead:
    calibmap=csv.reader(csvToRead)
    corr_matrix = list(calibmap)
csvToRead.close()

#convert the matrix back into floats

corr_map = []
for thetaSlice in corr_matrix:
    corrfloats = []
    for corr in thetaSlice:
        if float(corr)>0.:
            corrfloats.append(float(corr))
        else:
            corrfloats.append(-1.)
    corr_map.append(corrfloats)

fig, ax = plt.subplots()


#plot with pcolormesh
cmap = mpl.cm.get_cmap("plasma").copy()
cmap.set_under("white")
#endcap_mesh1 = ax.pcolormesh(EBins_endcap, endcap_bins, corr_map[0:2], cmap = cmap,vmin = 0., vmax = 3.5)
#endcap_mesh2 = ax.pcolormesh(EBins_endcap, endcap_bins_2, corr_map[len(corr_map)-3: len(corr_map)-1], cmap = cmap,
#        vmin = 0., vmax = 3.5)

mesh = ax.pcolormesh(EBins, ThetaBins, corr_map,cmap = cmap, 
        vmin = np.min(corr_map), vmax = np.max(corr_map))
#norm=mpl.colors.LogNorm(vmin=np.min(corr_map), vmax=np.max(corr_map)))
cbar = fig.colorbar(mesh)
cbar.set_label(label="$E_{true}/E_{reco}$",loc='top' )

#plt.xlabel("Summed ECal Hit E [GeV]", loc='right')
plt.xlabel("Reconstrcuted Photon E [GeV]", loc='right')
plt.ylabel("Reconstructed $\\theta$ [rad]", loc='top')
#use mplhep for labels
hep.cms.label(exp = "Muon Collider", data = False, 
       rlabel='$MAIA$ Detector Concept', loc=0, italic=(1,0,0), pad=(0.0))
#plt.gcf().text(0.1,0.9,"($\sqrt{s}=10$ TeV)")

ax.tick_params(bottom=True, top=True, left=True, right=True, which='both', direction='in')
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
ax.tick_params(axis= 'y',which='major',length = 10)
ax.tick_params(axis='y',which='minor', length=5)

ax.xaxis.set_major_locator(MultipleLocator(200))
ax.xaxis.set_major_formatter('{x:.0f}')
ax.xaxis.set_minor_locator(MultipleLocator(40))

ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.yaxis.set_major_formatter('{x:.1f}')
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
#hep.cms.label(exp = "Muon Collider", data = False, rlabel = "MAIA Detector Concept", loc=0, italic=(1,0,0))
plt.savefig('responseMap_photons_newregions_v08v2.pdf')
print("Created file 'responseMap_photons_newregions_v08v2.pdf'")
plt.close()
