import os
import logging
from ROOT import TH1D, TH2D, TFile, TTree, TColor, TCanvas, TLegend, TLatex, TLine, TMath, TGraphErrors, TF1, TProfile
from ROOT import kBlack, kBlue, kRed, kYellow, kGreen, kGray, kOrange, kViolet
from ROOT import gStyle, gPad
from ROOT import gROOT
from ROOT import TStyle
from optparse import OptionParser
import itertools
from math import *
from array import array


def check_output_directory(output_path):
    # checks if output directory exists; if not, mkdir
    if not os.path.exists(str(output_path)):
        os.makedirs(output_path)


# Options
parser = OptionParser()
#parser.add_option("-i", "--inFile_0_50",  dest='inFile_0_50',
                  #default="histos_photon.root", help="Name of the ROOT file")
#parser.add_option("-j", "--inFile_50_250",  dest='inFile_50_250',
                  #default="histos_photon.root", help="Name of the ROOT file")
#parser.add_option("-k", "--inFile_0_50_BIB",  dest='inFile_0_50_BIB',
                  #default="histos_photon.root", help="Name of the ROOT file")
#parser.add_option("-l", "--inFile_50_250_BIB",  dest='inFile_50_250_BIB',
                  #default="histos_photon.root", help="Name of the ROOT file")
parser.add_option("-m", "--inFile_0_50_10T", dest='inFile_0_50_10T',
                  default="histos_photon.root", help="Name of the ROOT file")
parser.add_option("-n", "--inFile_50_250_10T", dest='inFile_50_250_10T',
                  default="histos_photon.root", help="Name of the ROOT file")
parser.add_option("-o", "--inFile_250_1000_10T", dest='inFile_250_1000_10T',
                  default="histos_photon.root", help="Name of the ROOT file")
#parser.add_option("-p", "--inFile_250_1000",  dest='inFile_250_1000',
                  #default="histos_photon.root", help="Name of the ROOT file")
#parser.add_option("-q", "--inFile_1000_5000",  dest='inFile_1000_5000',
                  #default="histos_photon.root", help="Name of the ROOT file")

(options, args) = parser.parse_args()

# Load files
#fFile_0_50 = TFile(options.inFile_0_50, "READ")
#fFile_50_250 = TFile(options.inFile_50_250, "READ")
#fFile_0_50_BIB = TFile(options.inFile_0_50_BIB, "READ")
#fFile_50_250_BIB = TFile(options.inFile_50_250, "READ")
fFile_0_50_10T = TFile(options.inFile_0_50_10T, "READ")
fFile_50_250_10T = TFile(options.inFile_50_250_10T, "READ")
fFile_250_1000_10T = TFile(options.inFile_250_1000_10T, "READ")
#fFile_1000_5000 = TFile(options.inFile_1000_5000, "READ")


gStyle.SetPadTopMargin(0.09)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadBottomMargin(0.16)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetOptStat(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

arrBins_theta = array('d', (0., 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., TMath.Pi()))
arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 50.,
                        100., 250., 500., 1000., 1500., 2500., 5000.))
#arrBins_Eres = array('d', (0.,1.,2., 5.,7.5, 10.,12.5, 15., 20., 30., 40., 50., 60., 70., 80., 100., 125., 150., 175., 200., 250., 300., 350., 400.,
                           #500., 600., 700., 800., 900., 1000.))
arrBins_Eres = array('d', (0.,50.,100.,150.,200.,250.,300.,350.,400.,450.,500.,550.,600.,650.,700.,750.,800.,850.,
                           900.,950.,1000.))
#arrBins_Eres = array('d', (0., 5.,  10., 15., 20., 25., 50., 100., 250., 500., 1000., 1200., 2500., 5000.))

arrBins_response = array('d', (0., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 17., 20.,
                               25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90, 95., 100.,
                               110., 120., 130., 140., 150., 200., 250., 300., 350., 400., 450., 500., 550.,
                               600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400., 1500., 2000.,
                               2500., 3000., 4000., 5000.))

h_response_3 = TH2D('photon_response', 'photon_response',
                  1000, 0, 500, 1000, 0, 500)
h_response_10 = TH2D('photon_response_10', 'photon_response_10',
                  2000, 0, 1000, 100, -10, 10)
h_response_10_BIB = TH2D('photon_response_10bib', 'photon_response_10bib',
                  1000, 0, 500, 1000, 0, 500)
h_response_profile_3 = TProfile('photon_response_profile_3', 'photon_response_profile_3',
                              len(arrBins_response)-1, arrBins_response, 0, 5000)
h_response_profile_10 = TProfile('photon_response_profile_10', 'photon_response_profile_10',
                                 len(arrBins_Eres)-1, arrBins_Eres, 0, 5000)
h_response_profile_10_BIB = TProfile('photon_response_profile_10bib', 'photon_response_profile_10bib',
                              len(arrBins_response)-1, arrBins_response, 0, 5000)

h_deltaresponse_3 = TH2D('delta_response_3', 'delta_response_3',
                       100, 0, 500, 100, -1., 1.)
h_deltaresponse_10 = TH2D('delta_response_10', 'delta_response_10',
                       100, 0, 500, 100, -1., 1.)
h_deltaresponse_10_BIB = TH2D('delta_response_10bib', 'delta_response_10bib',
                       100, 0, 500, 100, -1., 1.)

#files_3T = [fFile_0_50, fFile_50_250]#, fFile_250_1000, fFile_1000_5000]
files_10T = [fFile_0_50_10T, fFile_50_250_10T, fFile_250_1000_10T]
#files_10T_BIB = [fFile_0_50_BIB]#, fFile_50_250_BIB]

h_reso_theta = TH1D('reso_theta', 'reso_theta',
                            len(arrBins_theta)-1, arrBins_theta)
h_reso_theta_20_50 = TH1D('reso_theta_20_50', 'reso_theta_20_50',
                          len(arrBins_theta)-1, arrBins_theta)
h_reso_theta_50_250 = TH1D('reso_theta_50_250', 'reso_theta_50_250',
                           len(arrBins_theta)-1, arrBins_theta)
h_reso_theta_250_1000 = TH1D('reso_theta_250_1000', 'reso_theta_250_1000',
                             len(arrBins_theta)-1, arrBins_theta)
h_reso_theta_1000_5000 = TH1D('reso_theta_1000_5000', 'reso_theta_1000_5000',
                              len(arrBins_theta)-1, arrBins_theta)
h1_reso_E = TH1D('reso_E', 'reso_E',
                 len(arrBins_Eres)-1, arrBins_Eres)
h2_reso_E = TH1D('reso_E_10', 'reso_E_10', 
                 len(arrBins_Eres)-1, arrBins_Eres)
h3_reso_E = TH1D('reso_E_BIB', 'reso_E_BIB',
                 len(arrBins_Eres)-1, arrBins_Eres)

e_arr = array('d')
e_arr_10 = array('d')
e_arr_bib = array('d')
sigma_arr = array('d')
sigma_arr_10 = array('d')
sigma_arr_bib = array('d')
mu_arr_10 = array('d')
e_err_arr = array('d')
e_err_arr_10 = array('d')
e_err_arr_bib = array('d')
sigma_err_arr = array('d')
sigma_err_arr_10 = array('d')
sigma_err_arr_bib = array('d')
mu_err_arr_10 = array('d')

cx = TCanvas("", "", 800, 600)
gStyle.SetOptStat(1)
'''
for bin in range(0, len(arrBins_Eres)-1):
    minE = arrBins_Eres[bin]
    maxE = arrBins_Eres[bin+1]

    proj_name = "ph_E"+str(arrBins_Eres[bin])+"_py"
    print(proj_name)
    binrange = 0
    if arrBins_Eres[bin]< 10.0:
        binrange = 1.5
    elif arrBins_Eres[bin]<70.0:
        binrange = 1.
    else:
        binrange = 0.5
    h_my_proj = TH1D(proj_name, proj_name, 150, -binrange, binrange)

    for file in files_3T:
        tree = file.Get("photon_tree")
        for entry in tree:
            h_response_3.Fill(entry.E_truth, entry.E)
            h_response_profile_3.Fill(entry.E_truth, entry.E)
            h_deltaresponse_3.Fill(
                entry.E_truth, (entry.E - entry.E_truth)/entry.E_truth)
            
        
        for entry in tree:
            if entry.E > 0:
                #if entry.theta < 0.61 or entry.theta > 2.53: #entry.theta < 0.612:
                if (entry.theta < 0.96 and entry.theta >0.61) or (entry.theta > 2.18 and entry.theta < 2.53):
                    if (entry.E_truth > minE) and (entry.E_truth < maxE):
                        if (entry.E - entry.E_truth)/entry.E_truth > -0.5:
                            h_my_proj.Fill((entry.E - entry.E_truth)/entry.E_truth)
                        
    lim = 0.0
    if minE<250:
        if minE<50:
            lim = 0.5
        else:
            lim = 0.25
        gaussFit = TF1("gaussfit", "gaus", -1.0*lim, lim)
        gaussFit.SetLineColor(kRed)
        gaussFit.SetParameter(1, 0.)
        gaussFit.SetParameter(2, h_my_proj.GetRMS())
        h_my_proj.Fit(gaussFit, "E")
        gStyle.SetOptFit(0o111);
        h_my_proj.Draw("HIST")
        gaussFit.Draw("LSAME")
        cx.SaveAs("slices_ph/"+proj_name+".pdf")
        sigma = gaussFit.GetParameter(2)
        sigma_err = gaussFit.GetParError(2)
        if bin > 0:
            e_arr.append(h1_reso_E.GetBinCenter(bin+1))
            e_err_arr.append((maxE-minE)/2)
            sigma_arr.append(sigma)
            print("SIGMA=",sigma)
            #print(h1_reso_E.GetBinCenter(bin+1))
            #print(sigma)
            sigma_err_arr.append(sigma_err)
        else:
            gaussFit = TF1("gaussfit", "gaus", -0.15, 0.15)
            gaussFit.SetLineColor(kRed)
            gaussFit.SetParameter(1, 0.)
            gaussFit.SetParameter(2, h_my_proj.GetRMS())
            h_my_proj.Fit(gaussFit, "E")
            h_my_proj.Draw("HIST")
            gaussFit.Draw("LSAME")
            cx.SaveAs("slices_ph/"+proj_name+".pdf")
            sigma = gaussFit.GetParameter(2)
            sigma_err = gaussFit.GetParError(2)
            if bin > 0:
                e_arr.append(h1_reso_E.GetBinCenter(bin+1))
                e_err_arr.append(maxE-minE)
                sigma_arr.append(sigma)
                sigma_err_arr.append(sigma_err)

'''
'''
for bin in range(0, len(arrBins_Eres)-1):
    minE = arrBins_Eres[bin]
    maxE = arrBins_Eres[bin+1]

    proj_name = "ph_E"+str(arrBins_Eres[bin])+"_py"
    print(proj_name)
    h_my_proj = TH1D(proj_name, proj_name, 150, -5.0, 8.0)

    for file in files_10T_BIB:
        tree = file.Get("photon_tree")
        for entry in tree:
            h_response_10_BIB.Fill(entry.E_truth, entry.E)
            h_response_profile_10_BIB.Fill(entry.E_truth, entry.E)
            h_deltaresponse_10_BIB.Fill(
                entry.E_truth, (entry.E - entry.E_truth)/entry.E_truth)
            
        
        for entry in tree:
            if entry.E > 0 and entry.theta < TMath.Pi()/4:
                if (entry.E_truth > minE) and (entry.E_truth < maxE):
                    if (entry.E - entry.E_truth)/entry.E_truth > -1.5:
                        h_my_proj.Fill((entry.E - entry.E_truth)/entry.E_truth)
                        

    if minE<50 and minE > 5:
        gaussFit = TF1("gaussfit", "gaus", -1.75, 1.75)
        gaussFit.SetLineColor(kRed)
        gaussFit.SetParameter(1, 0.)
        gaussFit.SetParameter(2, h_my_proj.GetRMS())
        h_my_proj.Fit(gaussFit, "E")
        gStyle.SetOptFit(0o111);
        h_my_proj.Draw("HIST")
        gaussFit.Draw("LSAME")
        cx.SaveAs("slices_ph/"+proj_name+".pdf")
        sigma = gaussFit.GetParameter(2)
        sigma_err = gaussFit.GetParError(2)
        if bin > 0:
            e_arr_bib.append(h3_reso_E.GetBinCenter(bin+1))
            e_err_arr_bib.append((maxE-minE)/2)
            sigma_arr_bib.append(sigma)
            print("SIGMA=",sigma)
            #print(h1_reso_E.GetBinCenter(bin+1))
            #print(sigma)
            sigma_err_arr_bib.append(sigma_err)
        else:
            gaussFit = TF1("gaussfit", "gaus", -0.15, 0.15)
            gaussFit.SetLineColor(kRed)
            gaussFit.SetParameter(1, 0.)
            gaussFit.SetParameter(2, h_my_proj.GetRMS())
            h_my_proj.Fit(gaussFit, "E")
            h_my_proj.Draw("HIST")
            gaussFit.Draw("LSAME")
            cx.SaveAs("slices_ph/"+proj_name+".pdf")
            sigma = gaussFit.GetParameter(2)
            sigma_err = gaussFit.GetParError(2)
            if bin > 0:
                e_arr_bib.append(h3_reso_E.GetBinCenter(bin+1))
                e_err_arr_bib.append(maxE-minE)
                sigma_arr_bib.append(sigma)
                sigma_err_arr_bib.append(sigma_err)



'''       
for bin in range(0, len(arrBins_Eres)-1):
    minE = arrBins_Eres[bin]
    maxE = arrBins_Eres[bin+1]

    proj_name = "ph_E"+str(arrBins_Eres[bin])+"_py"
    print(proj_name)
    binrange = 0
    if arrBins_Eres[bin]< 50.0:
        binrange = 1.5
    elif arrBins_Eres[bin]<100.0:
        binrange = 1.
    elif arrBins_Eres[bin]<500.0:
        binrange = 0.5
    else:
        binrange = 0.25
    h_my_proj = TH1D(proj_name, proj_name, 150, -binrange, binrange)

    for file in files_10T:
        tree = file.Get("photon_tree")
        for entry in tree:
            if entry.theta >1.01 and entry.theta < 2.13:
                h_response_10.Fill(entry.E_truth, entry.E/entry.E_truth)
                #h_response_10.Fill(entry.E_truth, entry.E)
                h_response_profile_10.Fill(entry.E_truth, (entry.E-entry.E_truth)/entry.E_truth)
                h_deltaresponse_10.Fill(
                    entry.E_truth, (entry.E - entry.E_truth)/entry.E_truth)
            
     
        for entry in tree:
            if entry.E > 0:
                #if(entry.theta > 0.65 and  entry.theta < 1.01) or (entry.theta > 2.13 and entry.theta <2.53): #if entry.theta < 0.612:
                if entry.theta >1.01 and entry.theta < 2.13:
                    if (entry.E_truth > minE) and (entry.E_truth < maxE):
                        if (entry.E - entry.E_truth)/entry.E_truth > -0.9:
                            h_my_proj.Fill((entry.E - entry.E_truth)/entry.E_truth)
                        
    lim = 0.0
    if minE<50000:
        if minE<50:
            lim = 0.5
        else:
            lim = 0.25
        gaussFit = TF1("gaussfit", "gaus", -1.0*lim, lim)
        gaussFit.SetLineColor(kBlue)
        gaussFit.SetParameter(1, 0.)
        gaussFit.SetParameter(2, h_my_proj.GetRMS())
        h_my_proj.Fit(gaussFit, "E")
        gStyle.SetOptFit(0o111);
        h_my_proj.Draw("HIST")
        #gaussFit.Draw("LSAME")
        cx.SaveAs("slices_ph/"+proj_name+".pdf")
        sigma = gaussFit.GetParameter(2)
        sigma_err = gaussFit.GetParError(2)
        if bin > 0:
            bincenter = (maxE-minE)/2
            e_arr_10.append(h2_reso_E.GetBinCenter(bin+1))
            e_err_arr_10.append(bincenter)
            sigma_arr_10.append(sigma)
            mu_arr_10.append(gaussFit.GetParameter(1))
            mu_err_arr_10.append(gaussFit.GetParError(1))
            print("SIGMA=",sigma)
            #print(h1_reso_E.GetBinCenter(bin+1))
            sigma_err_arr_10.append(sigma_err)
        else:
            gaussFit = TF1("gaussfit", "gaus", -0.15, 0.15)
            gaussFit.SetLineColor(kBlue)
            gaussFit.SetParameter(1, 0.)
            gaussFit.SetParameter(2, h_my_proj.GetRMS())
            h_my_proj.Fit(gaussFit, "E")
            h_my_proj.Draw("HIST")
            gaussFit.Draw("LSAME")
            cx.SaveAs("slices_ph/"+proj_name+".pdf")
            sigma = gaussFit.GetParameter(2)
            sigma_err = gaussFit.GetParError(2)
            if bin > 0:
                e_arr_10.append(h1_reso_E.GetBinCenter(bin+1))
                e_err_arr_10.append(maxE-minE)
                sigma_arr_10.append(sigma)
                sigma_err_arr_10.append(sigma_err)
        
        

#print(sigma_arr_10)
#h_reso_E_3 = TGraphErrors(len(e_arr), e_arr, sigma_arr, e_err_arr, sigma_err_arr)
h_reso_E_10 = TGraphErrors(len(e_arr_10), e_arr_10, sigma_arr_10, e_err_arr_10, sigma_err_arr_10)
#h_response_E_10 = TGraphErrors(len(e_arr_10), e_arr_10, mu_arr_10, e_err_arr_10, mu_err_arr_10)
#h_reso_E_10_bib = TGraphErrors(len(e_arr_bib), e_arr_bib, sigma_arr_bib, e_err_arr_bib, sigma_err_arr_bib)

c2_2 = TCanvas("", "", 800, 600)
#h_response_10.Draw("HIST")
#c2_2.Update()
#c2_2.SaveAs("2dhist.root")

'''
h_reso_E_3.SetTitle(" ")
#h_reso_E.SetMaximum(1.2)
#h_reso_E.SetMinimum(0.005)
h_reso_E_3.GetYaxis().SetTitle("Photon #sigma_{E}/E")
h_reso_E_3.GetYaxis().SetTitleOffset(1.4)
h_reso_E_3.GetXaxis().SetTitleOffset(1.2)
h_reso_E_3.GetXaxis().SetTitle("True photon energy [GeV]")
'''

h_reso_E_10.SetTitle(" ")
#h_reso_E.SetMaximum(1.2)
#h_reso_E.SetMinimum(0.005)
h_reso_E_10.GetYaxis().SetTitle("Photon #sigma_{E}/E")
h_reso_E_10.GetYaxis().SetTitleOffset(1.4)
h_reso_E_10.GetXaxis().SetTitleOffset(1.2)
h_reso_E_10.GetXaxis().SetTitle("True photon energy [GeV]")

'''
h_reso_E_10_bib.SetTitle(" ")
#h_reso_E.SetMaximum(1.2)
#h_reso_E.SetMinimum(0.005)
h_reso_E_10_bib.GetYaxis().SetTitle("Photon #sigma_{E}/E")
h_reso_E_10_bib.GetYaxis().SetTitleOffset(1.4)
h_reso_E_10_bib.GetXaxis().SetTitleOffset(1.2)
h_reso_E_10_bib.GetXaxis().SetTitle("True photon energy [GeV]")

'''
'''
resoFit3 = TF1(
    "resofit", "sqrt([0]*[0]/x+[1]*[1]/(x*x)+[2]*[2])", 0., 250., 3)
resoFit3.SetLineColor(kRed)
resoFit3.SetParName(0, "Stochastic (3TeV)")
resoFit3.SetParName(1, "Noise (3TeV)")
resoFit3.SetParName(2, "Constant (3TeV)")
resoFit3.SetParameter(0, 0.1)
resoFit3.FixParameter(1, 0.)
resoFit3.SetParameter(2, 0.01)
gStyle.SetOptFit(0)
h_reso_E_3.Fit(resoFit3)
#resoFit3.SetLineColor(kRed)
'''
resoFit10 = TF1(
    "resofit", "sqrt([0]*[0]/x+[1]*[1]/(x*x)+[2]*[2])", 0., 1000., 3)
resoFit10.SetLineColor(kBlue)
resoFit10.SetParName(0, "Stochastic (NCS, 10TeV)")
resoFit10.SetParName(1, "Noise (NCS, 10TeV)")
resoFit10.SetParName(2, "Constant (NCS, 10TeV)")
resoFit10.SetParameter(0, 0.5)
resoFit10.FixParameter(1, 0.)
resoFit10.SetParameter(2, 0.)
gStyle.SetOptFit(0)
h_reso_E_10.Fit(resoFit10)
#resoFit10.SetLineColor(kBlue)

'''
resoFit10BIB = TF1(
    "resofit", "sqrt([0]*[0]/x+[1]*[1]/(x*x)+[2]*[2])", 0., 50., 3)
resoFit10BIB.SetLineColor(kRed)
resoFit10BIB.SetParName(0, "Stochastic (NCS, 10TeV, BIB)")
resoFit10BIB.SetParName(1, "Noise (NCS, 10TeV, BIB)")
resoFit10BIB.SetParName(2, "Constant (NCS, 10TeV, BIB)")
resoFit10BIB.SetParameter(0, 0.5)
resoFit10BIB.SetParameter(1, 0.)
resoFit10BIB.SetParameter(2, 0.)
h_reso_E_10_bib.Fit(resoFit10BIB)

'''
'''
h_reso_E_3.SetLineColor(kRed)
h_reso_E_3.SetLineWidth(2)
h_reso_E_3.Draw("AP")
'''

h_reso_E_10.SetLineColor(kBlue)
h_reso_E_10.SetLineWidth(2)
h_reso_E_10.Draw("AP")


'''
h_reso_E_10_bib.SetLineColor(kRed)
h_reso_E_10_bib.SetLineWidth(2)
h_reso_E_10_bib.Draw("AP")
'''

'''
resoFit3.Draw("LSAME")
h_reso_E_3.Draw("PSAME")
resoFit10.Draw("LSAME")
h_reso_E_10.Draw("PSAME")
#resoFit10BIB.Draw("LSAME")
#h_reso_E_10_bib.Draw("PSAME")
'''

#h_reso_E_10.GetXaxis().SetRangeUser(0.,1000.)
#h_reso_E_10.GetYaxis().SetRangeUser(0.01, 0.4)

c2_2.Update()
#c2_2.SetLogy()
#c2_2.SetLogx()

#h_reso_E_10.GetYaxis().SetRangeUser(0.01,0.4)
#c2_2.Update()

leg = TLegend(.65, 0.72, .95, 0.92)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.035)
#leg.AddEntry(h_reso_E_3, "#sqrt{s}=3 TeV", "L")
#leg.AddEntry(h_reso_E_10, "#sqrt{s}=10 TeV", "L")
#leg.AddEntry(h_reso_E_10_bib, "10 TeV, No cell sub., ?% BIB", "L")
leg.Draw()

t3 = TLatex()
t3.SetTextFont(42)
t3.SetTextColor(1)
t3.SetTextSize(0.025)
t3.SetTextAlign(12)
t3.SetNDC()

t4 = TLatex()
t4.SetTextFont(42)
t4.SetTextColor(1)
t4.SetTextSize(0.025)
t4.SetTextAlign(12)
t4.SetNDC()

t5 = TLatex()
t5.SetTextFont(42)
t5.SetTextColor(1)
t5.SetTextSize(0.025)
t5.SetTextAlign(12)
t5.SetNDC()

t6 = TLatex()
t6.SetTextFont(42)
t6.SetTextColor(1)
t6.SetTextSize(0.025)
t6.SetTextAlign(12)
t6.SetNDC()


#t3.DrawLatex(0.15,0.65, 'Stochastic (3TeV)  '+str(round(resoFit3.GetParameter(0),3))+' #pm '+str(round(resoFit3.GetParError(0),4))+'\n'+'Noise (3TeV)  '+str(round(resoFit3.GetParameter(1),3))+' #pm '+str(round(resoFit3.GetParError(1),4))+'\n'+'Constant (3TeV)  '+str(round(resoFit3.GetParameter(2),3))+' #pm '+str(round(resoFit3.GetParError(2),4)))

#t4.DrawLatex(0.15,0.35, 'Stochastic (10TeV BIB)  '+str(round(resoFit10BIB.GetParameter(0),3))+' #pm '+str(round(resoFit10BIB.GetParError(0),4))+'\n'+'Noise (10TeV BIB)  '+str(round(resoFit10BIB.GetParameter(1),3))+' #pm '+str(round(resoFit10BIB.GetParError(1),4))+'\n'+'Constant (10TeV BIB)  '+str(round(resoFit10BIB.GetParameter(2),3))+' #pm '+str(round(resoFit10BIB.GetParError(2),4)))


t3.DrawLatex(0.15, 0.94, 'No BIB, Barrel Region (1.01<#theta<2.13)')
t4.DrawLatex(0.82, 0.94, '#sqrt{s} = 10 TeV')

c2_2.SaveAs("reso_vs_E_correctedgeo_10T.root")
c2_2.SaveAs("reso_vs_E_10T_BARREL_v2_EPJC_bin.pdf")

'''
c1_1=TCanvas("", "", 800, 600)
#h_response_E_10.SetLineColor(kOrange)
#h_response_E_10.Draw("AP")
h_response_profile_10.GetYaxis().SetTitle("Photon #mu_{E}/E")
h_response_profile_10.GetYaxis().SetTitleOffset(1.4)
h_response_profile_10.GetXaxis().SetTitleOffset(1.2)
h_response_profile_10.GetXaxis().SetTitle("True photon energy [GeV]")
#h_response_profile_10.Draw("AEP")
#h_response_E_10.Draw("HIST")
h_response_profile_10.Draw("HIST PE")
h_response_profile_10.GetYaxis().SetRangeUser(-5.,5.)
c1_1.Update()
#h_response_E_10.Draw("PSAME")


leg_1 = TLegend(.65, 0.72, .95, 0.92)
leg_1.SetBorderSize(0)
leg_1.SetFillColor(0)
leg_1.SetFillStyle(0)
leg_1.SetTextFont(42)
leg_1.SetTextSize(0.035)
leg_1.AddEntry(h_response_profile_10, "Response Profile", "PE")
#leg_1.AddEntry(h_response_E_10, "Response Gaussian Means", "PE")
leg_1.Draw()
c1_1.Update()
#c1_1.SaveAs("response_profile.root")

#fFile_0_50.Close()
#fFile_50_250.Close()
#fFile_0_50_10T.Close()
#fFile_0_50_BIB.Close()
#fFile_50_250_10T.Close()
#fFile_250_1000.Close()
#fFile_1000_5000.Close()

'''
