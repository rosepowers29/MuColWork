from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TLorentzVector, TTree, TVector3
from ROOT import TH1F, TH1D, TH2D, TFile, TTree, TColor, TCanvas, TLegend, TLatex, TLine, TMath, TGraphErrors, TF1, TProfile, TProfile2D
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
from scipy.stats import crystalball
import scipy as scp
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


parser= OptionParser()
parser.add_option("-m", "--inFile_0_50_10T", dest='inFile_0_50_10T',
                  default="histos_photon.root", help="Name of the ROOT file")

parser.add_option("-n", "--inFile_50_250_10T", dest='inFile_50_250_10T',
                  default="histos_photon.root", help="Name of the ROOT file")

parser.add_option("-o", "--inFile_250_1000_10T", dest='inFile_250_1000_10T',
                 default="histos_photon.root", help="Name of the ROOT file")

#parser.add_option("-p", "--inFile_250_1000_10T_2", dest = 'inFile_250_1000_10T_2',
#                 default = "histos_photon.root", help="Name of the ROOT file")

parser.add_option("-r", "--region", dest = 'region', default = 'e', help = "theta region")
(options, args) = parser.parse_args()

#load files -- the slices
fFile_0_50 = TFile(options.inFile_0_50_10T, "READ")
fFile_50_250 = TFile(options.inFile_50_250_10T, "READ")
fFile_250_1000 = TFile(options.inFile_250_1000_10T, "READ")
#fFile_250_1000_2 = TFile(options.inFile_250_1000_10T_2, "READ")

#grab the region
region = options.region

# Energy binning (EPJC style), angular binning
EBins = array('d', (30.,50., 100., 150., 200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800.))#, 850., 900., 950., 1000.))
EBins_calib = array('d', (10.,50., 100., 150., 200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800.))#, 850., 900., 950., 1000.))
EBins_lowE = array('d', (10.,15., 20., 25.,30.,35., 40.,45., 50.))
#ThetaBins = array('d', (20.*TMath.Pi()/180, 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180., 90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., 160*TMath.Pi()/180))
EBins_endcap = array('d', (0.,100.,200.,300.,400.,500.,600.,700.,800.))
#tb_bins = np.linspace(0.6899,2.45,7)
#endcap_bins = np.linspace(0.175, 0.6899, 3)
#endcap_bins_2 = np.linspace(2.45, 2.96, 3)
#ThetaBins_neutrons = np.concatenate((endcap_bins[:2], tb_bins, endcap_bins_2[1:]))

#ThetaBins = np.linspace(0.175,2.96,30)
NPhotonBins = np.linspace(0,50,50)


endcap_bins = np.linspace(0.225, 0.625, 5)
endcap_bins_2 = np.linspace(2.517, 2.91, 5)
tb_bins = np.linspace(0.625,2.517,20)
ThetaBins = np.concatenate((endcap_bins[:2], tb_bins, endcap_bins_2[1:]))

#declare arrays
e_arr = array('d')
th_arr = array('d')
th_err_arr=array('d')
sigma_arr = array('d')
sigma_arr_th = array('d')
e_err_arr = array('d')
sigma_err_arr = array('d')
sigma_err_arr_th = array('d')
mu_arr = array('d')
mu_err_arr = array('d')
mu_arr_th = array('d')
mu_arr_th_err = array('d')


#get files
files = [fFile_0_50, fFile_50_250, fFile_250_1000]#, fFile_250_1000_2]


with open('../ResponseMaps/responseMap_neutrons_noBIB.csv', 'r') as csvToRead:
    calibmap=csv.reader(csvToRead)
    corr_matrix = list(calibmap)
csvToRead.close()

#convert the matrix back into floats

corr_map = []
for thetaSlice in corr_matrix:
    corrfloats = []
    for corr in thetaSlice:
        corrfloats.append(float(corr))
    corr_map.append(corrfloats)



################### RESOLUTION STUDY #################
cx = TCanvas("", "", 800, 600)
gStyle.SetOptStat(1)
for Ebin in range(0, len(EBins)-1):
    #ratios=[]
    if region == 'th':
        break
    EMin = EBins[Ebin]
    EMax = EBins[Ebin+1]
    proj_name = "E reso, "+str(EBins[Ebin])+"<E<"+str(EBins[Ebin+1])
    file_name = "Ereso"+str(Ebin)
    if EMin < 50:
        lim = 1.
        bins = 75
    elif 50<= EMin <= 200:
        lim = 2.
        bins = 70
    elif 200< EMin < 500.:
        lim = 1.5
        bins = 50
    else:
        lim = 1.
        bins = 50

    h_my_proj_2 = TH1D(proj_name, proj_name, bins, -lim, lim)

    for file in files:
        #tree = file.Get("jet_tree")
        #tree = file.Get("photon_tree")
        tree = file.Get("neutron_tree")
        for entry in tree:
            E_truth = entry.E_truth #entry.photon_E
            if E_truth < EMin or E_truth > EMax:
                continue
            theta = entry.theta_truth #entry.photon_theta_angle
            E_reco = entry.E #entry.jet_energy
            E_corr = entry.E #entry.jet_energy
            theta_reco = entry.theta #entry.jet_theta
            phi_truth = entry.phi_truth
            phi_reco = entry.phi
            #NPhoton = entry.N_photon
            #if(theta_reco > 0): #matched
                #h_theta_multiplicity.Fill(theta, NPhoton)
                #h_E_multiplicity.Fill(E_truth, NPhoton)

            dR_reco_true = np.sqrt((phi_reco-phi_truth)**2+(theta_reco-theta)**2) 
            #restrict to endcap region (NO SOLENOID CONTACT)
    
            if region == 'e':
                if theta > 0.625 and theta <2.52 or theta < 0.225 or theta > 2.91:
                    continue
            #restrict to central-barrel region (minimal solenoid contact)
            elif region == 'b':
                if theta < 0.99 or theta > 2.15:
                    continue
            #restrict to transition region (maximal solenoid contact)
            elif region == 't':
                if theta < 0.625 or theta > 2.225 or 0.99 <theta< 2.15:
                    continue
            #restrict to barrel and trandition region
            elif region == 'tb':
                if theta < 0.6899 or theta > 2.45:
                    continue
            #all-inclusive
            elif region == 'a':
                if theta < 0.18 or theta > 2.96:
                   continue
            #get rid of weird theta values
            if theta_reco < 0:
                continue
            # get rid of negative energy values, restrict range
            if E_reco < 10.:
                continue
            
            for tbin in range(0, len(ThetaBins)-1):
                ThMin = ThetaBins[tbin]
                ThMax = ThetaBins[tbin+1]
                foundCorr = False
                for Ecorrbin in range(0, len(EBins_calib)-1):
                    ECMin = EBins_calib[Ecorrbin]
                    ECMax = EBins_calib[Ecorrbin+1]
                    #get corr factor
                    corr = float(corr_matrix[tbin][Ecorrbin])
                    #if entry is in the 2d bin, correct energy
                    if ECMin <= E_reco and E_reco < ECMax:
                        foundCorr = True
                        if ThMin <= theta_reco and theta_reco < ThMax:
                            if corr != 0:# and EMin != 100.:
                                E_corr = E_reco * corr
                            else:
                                print(corr)
                    
                        
            neg_lim = -10.
            pos_lim = 10.
            #print((E_corr - E_truth)/E_truth)
            if (E_corr - E_truth)/E_truth > neg_lim and (E_corr - E_truth)/E_truth < pos_lim:
                h_my_proj_2.Fill((E_corr - E_truth)/E_truth)
    #now start fitting the gaussians
    
    lim = 0.0
    if EMin<50000:
        gaussFit1 = TF1("gaussfit", "gaus", -1., 1.)
        gaussFit1.SetLineColor(kRed)
        gaussFit1.SetParameter(1, 0.)
        gaussFit1.SetParameter(2, h_my_proj_2.GetRMS())
        h_my_proj_2.Fit(gaussFit1, "E")
        gStyle.SetOptFit(0o111);
        h_my_proj_2.Draw("HIST")
        gaussFit1.Draw("LSAME")
        cx.Update()
        #cx.SaveAs("slices_corr_full"+file_name+".root")
        sigma = gaussFit1.GetParameter(2)
        sigma_err = gaussFit1.GetParError(2)
        bincenter = (EMax-EMin)/2
        e_arr.append(EMin+bincenter)
        e_err_arr.append(bincenter)
        sigma_arr.append(sigma)
        mu_arr.append(gaussFit1.GetParameter(1))
        mu_err_arr.append(gaussFit1.GetParError(1))
        print("SIGMA=",sigma)
        #print(h1_reso_E.GetBinCenter(bin+1))
        sigma_err_arr.append(sigma_err)
# theta scan
if region == 'th':

    for Tbin in range(0, len(ThetaBins)-1):
        ThMin = ThetaBins[Tbin]
        ThMax = ThetaBins[Tbin+1]
        proj_name = "E reso, "+str(np.round(ThetaBins[Tbin],2))+"<theta<"+str(np.round(ThetaBins[Tbin+1],2))
        file_name = "Ereso"+str(Tbin)
        lim = 0.5
        if ThMin > 0.577 and ThMin < 2.96:
            lim = 2.
        #elif ThMin > 2.0 and ThMin < 2.58:
        #    lim = 5.
        else:
            lim = 5.
        h_my_proj_2 = TH1D(proj_name, proj_name, 50,-lim,lim)

        for file in files:
            #tree = file.Get("jet_tree")
            #tree = file.Get("photon_tree")
            tree = file.Get("neutron_tree")
            for entry in tree:
                E_truth = entry.E_truth #entry.photon_E
                theta = entry.theta_truth #entry.photon_theta_angle
                E_reco = entry.E #entry.jet_energy
                E_corr = entry.E #entry.jet_energy
                theta_reco = entry.theta #entry.jet_theta
                #restrict to barrel region
                #if theta < 0.678 or theta > 2.46:
                #    continue
                #get rid of weird theta values
                if theta_reco < 0:
                    continue
                # get rid of negative energy values, restrict range
                #only include the energy plateau
                if E_truth < 600.:
                    continue
                if E_reco < 5.:
                    continue
                for Ebin in range(0, len(EBins)-1):
                    EMin = EBins[Ebin]
                    EMax = EBins[Ebin+1]
                    #get corr factor from profile
                    corr = float(corr_matrix[Tbin][Ebin])
                    #if entry is in the 2d bin, correct energy
                    if EMin <= E_reco and E_reco < EMax:
                        if ThMin <= theta_reco and theta_reco < ThMax:
                            #if corr > 0:
                            E_corr = E_reco * corr
                            #isolate the peak of the distribution so we can fit a Gaussian
                            neg_lim = -4.
                            pos_lim = 4.
                            #print((E_corr - E_truth)/E_truth)
                            if (E_corr - E_truth)/E_truth > neg_lim and (E_corr - E_truth)/E_truth < pos_lim:
                                h_my_proj_2.Fill((E_corr - E_truth)/E_truth)
                                h_response_profile_E.Fill(E_truth, (E_truth)/E_corr)
                                h_response_profile_th.Fill(theta, (E_truth)/E_corr)
           

        #now start fitting the gaussians
        lim = 0.0
        if EMin<10000:
            gaussFit1 = TF1("gaussfit", "gaus", -1., 1.)
            gaussFit1.SetLineColor(kRed)
            gaussFit1.SetParameter(1, 0.)
            gaussFit1.SetParameter(2, h_my_proj_2.GetRMS())
            h_my_proj_2.Fit(gaussFit1, "E")
            gStyle.SetOptFit(0o111);
            h_my_proj_2.Draw("HIST")
            gaussFit1.Draw("LSAME")
            cx.Update()
            cx.SaveAs("slices_corr_full"+file_name+".root")
            sigma = gaussFit1.GetParameter(2)
            sigma_err = gaussFit1.GetParError(2)
            bincenter = (ThMax-ThMin)/2
            e_arr.append(h_reso_theta.GetBinCenter(Tbin+1))
            e_err_arr.append(bincenter)
            sigma_arr.append(sigma)
            mu_arr.append(gaussFit1.GetParameter(1))
            mu_err_arr.append(gaussFit1.GetParError(1))
            print("SIGMA=",sigma)
            #print(h1_reso_E.GetBinCenter(bin+1))
            sigma_err_arr.append(sigma_err)



#save arrays to a csv for easy access
res_info_noBIB = pd.DataFrame({'E': e_arr, 'sigma': sigma_arr, 'E_err': e_err_arr, 'sigma_err': sigma_err_arr})
res_info_noBIB.to_csv("../ResArrays/noBIB_v2_resarrays_neutrons_"+region+".csv")
#change from noBIB to BIB and vice versa depending on case

