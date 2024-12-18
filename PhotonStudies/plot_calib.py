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

EBins = array('d', (0., 50., 100., 150., 200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800.))#, 850., 900., 950., 1000.))
EBins_lowE = array('d', (0., 50., 100., 150., 200., 250.))#, 30., 40., 50.))
#ThetaBins = array('d', (20.*TMath.Pi()/180, 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180., 90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., 160*TMath.Pi()/180))
ThetaBins = np.linspace(0.175,2.96,30)


with open('responseMap_BIBPFOs_reco.csv', 'r') as csvToRead:
    calibmap=csv.reader(csvToRead)
    corr_matrix = list(calibmap)


#plot the CalibMaps in both cases
#convert the matrix back into floats

corr_map = []
for thetaSlice in corr_matrix:
    corrfloats = []
    for corr in thetaSlice:
        corrfloats.append(float(corr))
    corr_map.append(corrfloats)

fig, ax = plt.subplots()


#plot with pcolormesh
mesh = ax.pcolormesh(EBins, ThetaBins, corr_map,cmap = plt.cm.plasma, vmin = np.min(corr_map), vmax = np.max(corr_map))#norm=mpl.colors.LogNorm(vmin=np.min(corr_map), vmax=np.max(corr_map)))
cbar = fig.colorbar(mesh)
cbar.set_label(label="$E_{true}/E_{reco}$",loc='top' )

plt.xlabel("Reconstructed E [GeV]", loc='right')
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
plt.savefig('responseMap_recoPFOs_BIB.pdf')
print("Created file 'responseMap_recoPFOs_BIB.pdf'")
plt.close()
