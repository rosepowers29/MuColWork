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

#this is not really useful
######## ECAL dimensions #######
zlim = 230.7 #cm
R_inner_lim = 150.0 #cm
R_outer_lim = 185.7 #cm
X_0 = 8.897 #cm -- rad length for AL
###############################

####### ROOT gStyle ##########
gStyle.SetPadTopMargin(0.09)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadBottomMargin(0.16)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetOptStat(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
################################

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


# Declare response profiles and reso plots
h_response_profile_E = TProfile('E_profile', 'E_profile', len(EBins)-1, EBins, 150,800)
h_response_profile_th = TProfile('Th_profile', 'Th_profile', len(ThetaBins)-1, ThetaBins, 0, 50)
h_2d_response = TProfile2D('T-E_profile', 'T-E_profile', len(ThetaBins)-1, ThetaBins, len(EBins_lowE)-1, EBins_lowE)
h_2d_resp_corr = TProfile2D('TEcorr', 'TEcorr', len(ThetaBins)-1, ThetaBins, len(EBins)-1, EBins)
h_reso_E = TH1D('reso_E', 'reso_E', len(EBins)-1, EBins)
h_reso_theta = TH1D('reso_theta', 'reso_theta', len(ThetaBins)-1, ThetaBins)

h_E_v_theta = TH2D ('E_v_theta', 'E_v_theta', len(ThetaBins)-1, ThetaBins, len(EBins)-1, EBins)
h_theta_multiplicity = TH2D('N_vs_theta', 'N_vs_theta', len(ThetaBins)-1, ThetaBins, len(NPhotonBins)-1, NPhotonBins) 
h_E_multiplicity = TH2D('N_vs_E', 'N_vs_E', len(EBins)-1, EBins, len(NPhotonBins)-1, NPhotonBins)  

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

#none of the commented chunk below is used anymore
'''
#Define piecewise function for radlengths as fn of theta
def NX_0(theta):
    if 0.577 < theta and theta < 0.678:
        return (zlim/cos(theta)-R_inner_lim/sin(theta))/X_0
    elif 0.678 < theta and theta < TMath.Pi()-0.678:
        return ((R_outer_lim-R_inner_lim)/sin(theta))/X_0
    elif TMath.Pi()-0.678 < theta and theta < TMath.Pi()-0.577:
        return (zlim/cos(TMath.Pi()-theta)-R_inner_lim/sin(TMath.Pi()-theta))/X_0
    else:
        return 0

def corr_ratio(x,A):
    t = NX_0(x)
    return (2**t)*A
'''

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
'''
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
            #tree = file.Get("photon_tree")
            tree = file.Get("jet_tree")
            for entry in tree:
                E_truth = entry.photon_theta
                theta = entry.photon_theta_angle
                E_reco = entry.jet_energy
                theta_reco = entry.jet_theta
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
with open('calibMap_PFOs', 'w', newline = '') as csvfile:
    mapwriter = csv.writer(csvfile)
    mapwriter.writerows(corr_matrix)

#after initial calibrtion, read in the csv file and just pull values from there
# have block below uncommented if reading in file: corr_matrix will be defined differently
'''
with open('responseMap_photons_newregions_noBIB.csv', 'r') as csvToRead:
    calibmap=csv.reader(csvToRead)
    corr_matrix = list(calibmap)
csvToRead.close()

#plot the CalibMaps in both cases
#convert the matrix back into floats

corr_map = []
for thetaSlice in corr_matrix:
    corrfloats = []
    for corr in thetaSlice:
        corrfloats.append(float(corr))
    corr_map.append(corrfloats)



'''
#plot with pcolormesh
plt.pcolormesh(EBins, ThetaBins, corr_map,colormap = plt.cm.plasma, vmin=np.min(corr_map), vmax=np.max(corr_map))
plt.colorbar("$E_{jet}/E_{true}$", loc='top')
plt.xlabel("True E [GeV]",loc='right')
plt.ylabel("True Theta [rad]",loc='top')
#use mplhep for labels
hep.cms.label(exp = "Muon Collider", data = False, label = "(2023 Lattice)", 
       rlabel='$MAIA$ Detector Concept', loc=0, italic=(1,0,0), pad=(0.0))
plt.gcf().text(0.1,0.9,"($\sqrt{s}=10$ TeV)")
ax.tick_params(which='both', direction='in')
#hep.cms.label(exp = "Muon Collider", data = False, rlabel = "MAIA Detector Concept", loc=0, italic=(1,0,0))
plt.savefig('calibmap_PFOS.pdf')
plt.close()
print("calibmap done")


################### RESOLUTION STUDY #################
'''
profile=[]
stdevs = []
cx = TCanvas("", "", 800, 600)
gStyle.SetOptStat(1)
ratios=[]
bad_dRs = []
peak_dRs = []
bad_nPFO = []
peak_nPFO = []
badthetas = []
goodthetas = []
goodresponse = []
badresponse = []
goodphis=[]
badphis=[]
goodEreco=[]
badEreco=[]
for Ebin in range(0, len(EBins)-1):
    #ratios=[]
    if region == 'th':
        break
    EMin = EBins[Ebin]
    EMax = EBins[Ebin+1]
    proj_name = "E reso, "+str(EBins[Ebin])+"<E<"+str(EBins[Ebin+1])
    file_name = "Ereso"+str(Ebin)
    if EMin < 100:
        h_my_proj_2 = TH1D(proj_name, proj_name, 70, -10., 10.)
    elif 100<= EMin < 200:
        h_my_proj_2 = TH1D(proj_name, proj_name, 70, -5., 5.)
    elif 200<= EMin < 500.:
        h_my_proj_2 = TH1D(proj_name, proj_name, 50, -1.,1.)
    else:
        h_my_proj_2 = TH1D(proj_name, proj_name, 50, -1. ,1.)


    for file in files:
        #tree = file.Get("jet_tree")
        tree = file.Get("photon_tree")
        #tree = file.Get("neutron_tree")
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

            #cut out nozzles
            #if theta < 0.18 or theta > 2.96:
            #    continue
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
                        #isolate the peak of the distribution so we can fit a Gaussian
                        #BIB limits
                        
            if 30. == EMin:
                neg_lim = -20.
                pos_lim = 20.
            #elif 300. <= EMin < 1000.:
            #    neg_lim = -1.
            #    pos_lim = 1.
            else:
                neg_lim = -10.
                pos_lim = 10.
                        
            #neg_lim = -10.
            #pos_lim = 10.
            #print((E_corr - E_truth)/E_truth)
            if (E_corr - E_truth)/E_truth > neg_lim and (E_corr - E_truth)/E_truth < pos_lim:
                h_my_proj_2.Fill((E_corr - E_truth)/E_truth)
                if EMin == 30:
                    if -0.6 < (E_corr-E_truth)/E_truth < -0.2:
                        bad_dRs.append(dR_reco_true)
                        badthetas.append(theta_reco)
                        badresponse.append(E_truth/E_reco)
                        badphis.append(phi_reco)
                        badEreco.append(E_reco)
                    else:
                        goodthetas.append(theta_reco)
                        goodresponse.append(E_truth/E_reco)
                        peak_dRs.append(dR_reco_true)
                        goodphis.append(phi_reco)
                        goodEreco.append(E_reco)
                if EMin <= 150:
                    h_response_profile_E.Fill(E_truth, (E_truth)/E_reco)
                    ratios.append(E_truth/E_reco)

    if EMin == 30:
        plt.hist(bad_dRs, bins=40,range=(0,.04),color="red",histtype='step',density=True, label="Neg Peak dR")
        plt.hist(peak_dRs, bins=40,range=(0,.04),color="green", histtype='step',density=True,label="0 Peak dR")
        plt.xlabel("dR(truth, reco)")
        plt.legend()
        plt.savefig("dRs_30_"+region+".pdf")
        plt.close()

        plt.hist(badthetas, bins=50,range=(0.225,2.91),color="red",histtype='step',density=True, label="Neg Peak $\\theta_{reco}$")
        plt.hist(goodthetas, bins=50,range=(0.225,2.91),color="green", histtype='step',density=True,label="0 Peak $\\theta_{reco}$")
        plt.xlabel("$\\theta_{reco}$")
        plt.legend()
        plt.savefig("thetas_30_"+region+".pdf")
        plt.close()

        plt.hist(badresponse, bins=30,range=(1.,2.4),color="red",histtype='step', density=True,label="Neg Peak Response")
        plt.hist(goodresponse, bins=30,range=(1.,2.4),color="green", histtype='step',density=True,label="0 Peak Response")
        plt.xlabel("$E_{true}/E_{reco}$")
        plt.legend()
        plt.savefig("responses_30_"+region+".pdf")
        plt.close()
        
        plt.hist(badphis, bins=30,range=(0.,2*np.pi),color="red",histtype='step', density=True,label="Neg Peak Response")
        plt.hist(goodphis, bins=30,range=(0.,2*np.pi),color="green", histtype='step',density=True,label="0 Peak Response")
        plt.xlabel("$\phi_{reco}")
        plt.legend()
        plt.savefig("phis_30_"+region+".pdf")
        plt.close()
        
        plt.hist(badEreco, bins=50,range=(0.,100.),color="red",histtype='step', density=True,label="Neg Peak Response")
        plt.hist(goodEreco, bins=50,range=(0.,100.),color="green", histtype='step',density=True,label="0 Peak Response")
        plt.xlabel("$E_{reco}$")
        plt.legend()
        plt.savefig("Ereco_30_"+region+".pdf")
        plt.close()
        
        print("done")
    '''
    dat, bins, patches = plt.hist(ratios,bins=15,range=(0.8,3),color="green",label="Endcap")
    #popt, pcov = scp.optimize.curve_fit(crystalball, bins, dat, p0=[0.5,1])
    plt.scatter(bins,crystalball.pdf(bins,0.5, 1.5))
    plt.show()


    plt.xlabel("E_true/E_reco")
    plt.legend()
    plt.savefig("distr_"+region+".pdf")
    plt.close()
    profile.append(np.average(ratios))
    stdevs.append(np.std(ratios))
    #now start fitting the gaussians
    '''
    lim = 0.0
    if EMin<50000:
        if EMin<400:
            gaussFit1 = TF1("gaussfit", "gaus", -1., 1.)
        #else:
        #    gaussFit1 = TF1("gaussfit", "gaus", -1., -0.8)
        #elif EMin>700:
        #    gaussFit1 = TF1("gaussfit", "gaus", -.02, 0.01)   
        elif EMin == 300:
            gaussFit1 = TF1("gaussfit", "gaus", -1., 0.02)
        else:
            gaussFit1 = TF1("gaussfit", "gaus", -0.1, 0.1)
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
        bincenter = (EMax-EMin)/2
        e_arr.append(EMin+bincenter)
        e_err_arr.append(bincenter)
        sigma_arr.append(sigma)
        mu_arr.append(gaussFit1.GetParameter(1))
        mu_err_arr.append(gaussFit1.GetParError(1))
        print("SIGMA=",sigma)
        #print(h1_reso_E.GetBinCenter(bin+1))
        sigma_err_arr.append(sigma_err)


#print(len(profile))
#print(len(EBins[1:]))
#plt.errorbar(EBins[1:], profile, yerr = stdevs, ls = '', marker = "o", c='black')
#plt.xlabel("True E [GeV]")
#plt.ylabel("E_true/E_reco")
#plt.ylim(0,2)
#plt.savefig("barrel_profile.pdf")


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
            if EMin<400:
                gaussFit1 = TF1("gaussfit", "gaus", -1., 1.)
            #else:
            #   gaussFit1 = TF1("gaussfit", "gaus", -1., -0.8)
            #elif EMin>700:
            #    gaussFit1 = TF1("gaussfit", "gaus", -.02, 0.02)   
            else:
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




#declare the graphs
#gr_reso_E = TGraphErrors(len(e_arr), e_arr, sigma_arr, e_err_arr, sigma_err_arr)
#gr_reso_Th = TGraphErrors(len(th_arr), th_arr, sigma_arr_th,th_err_arr, sigma_err_arr_th)
#gr_response_E = TGraphErrors(len(e_arr), e_arr, mu_arr, e_err_arr, mu_err_arr)
#gr_response_Th = TGraphErrors(len(th_arr), th_arr, mu_arr_th, th_err_arr, mu_arr_th_err)



#print("E:", e_arr)
#print("sigma:", sigma_arr)
#print("sig error:", sigma_err_arr)
#print("E error:", e_err_arr)

#save arrays to a csv for easy access
res_info_noBIB = pd.DataFrame({'E': e_arr, 'sigma': sigma_arr, 'E_err': e_err_arr, 'sigma_err': sigma_err_arr})
res_info_noBIB.to_csv("ccBIB_v2_resarrays_photons_"+region+".csv")
#change from noBIB to BIB and vice versa depending on case
'''
plt.scatter(e_arr, sigma_arr, c='red')
plt.xlabel("Theta [rad]")
plt.ylabel("Resolution")
plt.savefig("PFOcorrected_theta_reso.pdf")


cE1=TCanvas("", "", 800, 600)
gr_reso_E.SetTitle(" ")
gr_reso_E.GetYaxis().SetTitle("Jet #sigma_{E}/E")
gr_reso_E.GetYaxis().SetTitleOffset(1.4)
gr_reso_E.GetXaxis().SetTitleOffset(1.2)
gr_reso_E.GetXaxis().SetTitle("True jet energy [GeV]")
gr_reso_E.SetLineColor(kBlue)
gr_reso_E.SetLineWidth(2)
gr_reso_E.Draw("AP")
cE1.Update()


#ROOT fit

resoFit = TF1(    "resofit", "sqrt([0]*[0]/x+[1]*[1]/(x*x)+[2]*[2])", 30., 1000., 3)
resoFit.SetLineColor(kBlue)
resoFit.SetParName(0, "Stochastic")
resoFit.SetParName(1, "Noise")
resoFit.SetParName(2, "Constant")
resoFit.SetParameter(0, 0.5)
resoFit.FixParameter(1, 0.)
resoFit.SetParameter(2, 0.)
gStyle.SetOptFit(0)
gr_reso_E.Fit(resoFit)
#resoFit10.SetLineColor(kBlue)
#resoFit.Write("ResoFit")
gr_reso_E.Draw("AP")
resoFit.Draw("LSAME")
#gr_reso_E.SaveAs("Ereso_gr_th.root")
#resoFit.SaveAs("Ereso_fit_th.root")
cE1.Update()
cE1.SaveAs("Ereso_endcap_matrixcorr.root")


cN1 = TCanvas("", "", 800, 600)
h_theta_multiplicity.SetTitle("")
h_theta_multiplicity.SetYTitle("NPhotons")
h_theta_multiplicity.SetXTitle("True Photon #theta [rad]")
h_theta_multiplicity.Draw("COLZ")
cN1.SetLogz()
cN1.SetRightMargin(0.18)
cN1.Update()
cN1.SaveAs("NPhoton_vs_th.root")

cN2 = TCanvas("", "", 800, 600)
h_E_multiplicity.SetTitle("")
h_E_multiplicity.SetYTitle("NPhotons")
h_E_multiplicity.SetXTitle("True Photon E [GeV]")
h_E_multiplicity.Draw("COLZ")
cN2.SetRightMargin(0.18)
cN2.Update()
cN2.SaveAs("NPhoton_vs_E.root")
'''
#MAIA format python plots (!!)
#not actually used anymore -- just the skeleton of how we do it
#redefine the fit function
'''
def py_resofit(stoch, noise, const, E):
    return np.sqrt((stoch**2)/E+(noise/E)**2+const**2)

Epoints = np.linspace(2.,1000,1000)
resopoints=[]
#grab the fit values from the ROOT fit
for energy in Epoints:
    resopoints.append(py_resofit(resoFit.GetParameter(0), resoFit.GetParameter(1), resoFit.GetParameter(2), energy))
print("E:", e_arr)
print("sigma:", sigma_arr)
print("sig error:", sigma_err_arr)
print("E error:", e_err_arr)
print("stoch:", resoFit.GetParameter(0), ", noise:", resoFit.GetParameter(1), ", const:", resoFit.GetParameter(2))

plt.errorbar(e_arr, sigma_arr, yerr = sigma_err_arr, xerr = e_err_arr, fmt=".", color='red')
plt.plot(Epoints, resopoints,'--', c='blue', label="Fit")
plt.ylim(0,0.14)
plt.xlabel("True Photon Energy [GeV]")
plt.ylabel("Photon Jet Energy Resolution $\sigma/E$")
hep.cms.label(exp = "Muon Collider", data = False, rlabel = "MAIA Detector Concept", loc=0, italic=(1,0,0))
plt.show()
plt.savefig("Ereso_EE_maia_PFOS_BIB.pdf")
plt.close()



cE2=TCanvas("", "", 800, 600)
#gr_response_E.SetTitle(" ")
#gr_response_E.GetYaxis().SetTitle("Photon #mu_{E}/E")
#gr_response_E.GetYaxis().SetTitleOffset(1.4)
#gr_response_E.GetXaxis().SetTitleOffset(1.2)
#gr_response_E.GetXaxis().SetTitle("True jet energy [GeV]")
#gr_response_E.SetLineColor(kOrange)
#gr_response_E.SetLineWidth(2)
#gr_response_E.Draw("AP")
h_response_profile_E.SetTitle("")
h_response_profile_E.SetAxisRange(0.,5.,"Y")
cE2.Update()
h_response_profile_E.GetYaxis().SetTitle("Neutron E_{true}/E_{reco}")
h_response_profile_E.GetYaxis().SetTitleOffset(1.4)
h_response_profile_E.GetXaxis().SetTitleOffset(1.2)
h_response_profile_E.GetXaxis().SetTitle("True neutron energy [GeV]")
h_response_profile_E.SetLineColor(kRed)
h_response_profile_E.SetLineWidth(2)
h_response_profile_E.Draw("HIST")
cE2.Update()
cE2.SaveAs("Eresp_neutrons.root")



cTh1=TCanvas("", "", 800, 600)
gr_reso_Th.SetTitle(" ")
gr_reso_Th.GetYaxis().SetTitle("Photon #sigma_{E}/E")
gr_reso_Th.GetYaxis().SetTitleOffset(1.4)
gr_reso_Th.GetXaxis().SetTitleOffset(1.2)
gr_reso_Th.GetXaxis().SetTitle("Photon #theta [rad]")
gr_reso_Th.SetLineColor(kBlue)
gr_reso_Th.SetLineWidth(2)
gr_reso_Th.Draw("AP")
cTh1.Update()
cTh1.SaveAs("ThReso_corr.root")


cTh2=TCanvas("", "", 800, 600)
#gr_response_Th.SetTitle(" ")
#gr_response_Th.GetYaxis().SetTitle("Photon #mu_{E}/E")
#gr_response_Th.GetYaxis().SetTitleOffset(1.4)
#gr_response_Th.GetXaxis().SetTitleOffset(1.2)
#gr_response_Th.GetXaxis().SetTitle("Photon #theta [rad]")
#gr_response_Th.SetLineColor(kOrange)
#gr_response_Th.SetLineWidth(2)
#gr_response_Th.Draw("AP")
h_response_profile_th.SetTitle("")
h_response_profile_th.GetYaxis().SetTitle("Jet E_{true}/E_{reco}")
h_response_profile_th.GetYaxis().SetTitleOffset(1.4)
h_response_profile_th.GetXaxis().SetTitleOffset(1.2)
h_response_profile_th.GetXaxis().SetTitle("True jet theta [rad]")
h_response_profile_th.SetLineColor(kRed)
h_response_profile_th.SetLineWidth(2)
h_response_profile_th.Draw("HIST PE")
cTh2.Update()
cTh2.SaveAs("ThResp_corr_BIB.root")
'''
