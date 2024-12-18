from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TLorentzVector, TTree, TVector3
from ROOT import TH1D, TH2D, TFile, TTree, TColor, TCanvas, TLegend, TLatex, TLine, TMath, TGraphErrors, TF1, TProfile, TProfile2D
from ROOT import kBlack, kBlue, kRed, kYellow, kGreen, kGray, kOrange, kViolet
from ROOT import gStyle, gPad
from ROOT import gROOT
from ROOT import TStyle
from math import *
import numpy as np
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

parser= OptionParser()
parser.add_option("-m", "--inFile_0_50_10T", dest='inFile_0_50_10T',
                  default="histos_photon.root", help="Name of the ROOT file")

parser.add_option("-n", "--inFile_50_250_10T", dest='inFile_50_250_10T',
                  default="histos_photon.root", help="Name of the ROOT file")

parser.add_option("-o", "--inFile_250_1000_10T", dest='inFile_250_1000_10T',
                 default="histos_photon.root", help="Name of the ROOT file")
#parser.add_option("-p", "--inFile_250_1000_10T_2", dest = 'inFile_250_1000_10T_2',
#                 default = "histos_photon.root", help = "Name of the ROOT file")
(options, args) = parser.parse_args()

#load files -- the slices
fFile_0_50 = TFile(options.inFile_0_50_10T, "READ")
fFile_50_250 = TFile(options.inFile_50_250_10T, "READ")
fFile_250_1000 = TFile(options.inFile_250_1000_10T, "READ")
#fFile_250_1000_2 = TFile(options.inFile_250_1000_10T_2, "READ")


# Energy binning (EPJC style), angular binning
EBins = array('d', (0., 50., 100., 150., 200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800.))#, 850., 900., 950., 1000.))
EBins_lowE = array('d', (10.,15., 20., 25.,30.,35., 40.,45., 50.))
#ThetaBins = array('d', (20.*TMath.Pi()/180, 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180., 90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., 160*TMath.Pi()/180))
ThetaBins = np.linspace(0.175,2.96,30)



#get files
files = [fFile_0_50, fFile_50_250, fFile_250_1000]#, fFile_250_1000_2]


#####################################
# BINNED CALIBRATION                #
# double for loop, take the         #
# Etrue/Erec ratio for each theta/E #
# slice and apply to the 2D bins.   #
#####################################

corr_matrix = [] #for each theta bin, a list of avgs for each E bin
#uncomment the below to regenerate the calibration matrix

for tbin in range(0, len(ThetaBins)-1):
    rmax = 0
    ThMin = ThetaBins[tbin]
    ThMax = ThetaBins[tbin+1]
    Elist = []
    for Ebin in range(0, len(EBins)-1):
        EMin = EBins[Ebin]
        EMax = EBins[Ebin+1]
        ratios = []
        for file in files:
            tree = file.Get("photon_tree")
            #tree = file.Get("jet_tree")
            for entry in tree:
                E_truth = entry.E_truth
                theta = entry.theta_truth
                E_reco = entry.E
                theta_reco = entry.theta
                #get rid of weird theta values
                if theta_reco < 0:
                    continue
                # get rid of negative energy values
                if E_reco < 0:
                    continue
                #get rid of lower tail, file-specific
                #if file == fFile_50_250 and E_reco < 35:
                #    continue
                #if file == fFile_250_1000 and E_reco < 100:
                #    continue
                #finally, dump the ratio into a list
                #use RECO values for Binning!!
                if EMin <= E_reco and E_reco < EMax:
                    if ThMin <= theta_reco and theta_reco < ThMax:
                        ratios.append(E_truth/E_reco)
        if len(ratios) > 0:
            avg_ratio = np.average(ratios)
        else:
            print("No entries in bin: ", Ebin, " Emin: ", EMin) 
            avg_ratio = 0
        Elist.append(avg_ratio)
    print("Theta slice ", tbin, " done")
    print("Starting theta slice ",tbin+1)
    corr_matrix.append(Elist)


#write calibration map to a csv file
with open('responseMap_BIBPFOs_reco.csv', 'w', newline = '') as csvfile:
    mapwriter = csv.writer(csvfile)
    mapwriter.writerows(corr_matrix)
