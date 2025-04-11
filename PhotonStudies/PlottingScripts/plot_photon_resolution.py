import os
import logging
from ROOT import TH1D, TH2D, TFile, TTree, TColor, TCanvas, TLegend, TLatex, TLine, TMath, TEfficiency, TF1
from ROOT import kBlack, kBlue, kRed, kYellow, kGreen, kGray, kOrange
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
parser.add_option("-i", "--inFile",   dest='inFile',
                  default="ntup_tracks.root", help="Name of the ROOT file")
parser.add_option("-o", "--outFolder",   dest='outFolder',
                  default="/data", help="Name of the output folder")
(options, args) = parser.parse_args()

# Load files
fFile = TFile(options.inFile, "READ")

gStyle.SetPadTopMargin(0.09)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadBottomMargin(0.16)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetOptStat(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

tree = fFile.Get("photon_tree")

arrBins_theta = array('d', (0., 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., TMath.Pi()))
arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 50.,
                  100., 250., 500., 1000., 2500., 5000.))

h_resolution = TH2D('photon_resolution_theta', 'photon_resolution_theta', len(
    arrBins_theta)-1, arrBins_theta, 100, -400, 200)


c0 = TCanvas("", "", 800, 600)
h_frame = TH1D('framec3', 'framec3', len(arrBins_E)-1, arrBins_E)

# Efficiency vs E
h_reco_E = fFile.Get("matched_E")
h_truth_E = fFile.Get("truth_E")

eff_rec = TEfficiency(h_reco_E, h_truth_E)
eff_rec.SetLineWidth(2)
eff_rec.SetLineColor(kRed)

h_frame.SetLineColor(kBlack)
h_frame.SetLineWidth(2)
h_frame.SetTitle("")
h_frame.GetYaxis().SetTitle("Photon reconstruction efficiency")
h_frame.GetYaxis().SetTitleOffset(1.2)
h_frame.GetXaxis().SetTitleOffset(1.2)
h_frame.GetXaxis().SetRangeUser(0., 500)
h_frame.GetXaxis().SetTitle("Truth photon E [GeV]")
h_frame.SetMinimum(0.)
h_frame.SetMaximum(1.05)

h_frame.Draw()
eff_rec.Draw("E0SAME")
t2_3 = TLatex()
t2_3.SetTextFont(42)
t2_3.SetTextColor(1)
t2_3.SetTextSize(0.035)
t2_3.SetTextAlign(12)
t2_3.SetNDC()
t2_3.DrawLatex(.62, 0.31, 'Photon particle gun')
t2_3.DrawLatex(.62, 0.25, 'Cell threshold 2 MeV')

t3 = TLatex()
t3.SetTextFont(42)
t3.SetTextColor(1)
t3.SetTextSize(0.035)
t3.SetTextAlign(12)
t3.SetNDC()
t3.DrawLatex(0.15, 0.94, 'Background hits overlay in [-0.5, 10] ns range')

t4 = TLatex()
t4.SetTextFont(42)
t4.SetTextColor(1)
t4.SetTextSize(0.035)
t4.SetTextAlign(12)
t4.SetNDC()
t4.DrawLatex(0.82, 0.94, '#sqrt{s} = 10 TeV')

c0.SaveAs("photon_efficiency_vs_E.pdf")


# Efficiency vs theta
h_frame2 = TH1D('framec', 'framec', len(arrBins_theta)-1, arrBins_theta)

c = TCanvas("", "", 800, 600)

h_reco_theta = fFile.Get("matched_theta")
h_truth_theta = fFile.Get("truth_theta")

eff_rec = TEfficiency(h_reco_theta, h_truth_theta)
eff_rec.SetLineColor(kRed)
eff_rec.SetLineWidth(2)

h_frame2.SetLineColor(kBlack)
h_frame2.SetLineWidth(2)
h_frame2.SetTitle("")
h_frame2.GetYaxis().SetTitle("Photon reconstruction efficiency")
h_frame2.GetYaxis().SetTitleOffset(1.2)
h_frame2.GetXaxis().SetTitleOffset(1.2)
# h_frame2.GetXaxis().SetRangeUser(-0.24, 0.5)
h_frame2.GetXaxis().SetTitle("Truth photon #theta [rad]")
h_frame2.SetMinimum(0.8)
h_frame2.SetMaximum(1.05)

h_frame2.Draw()
eff_rec.Draw("E0SAME")

t2_3.DrawLatex(.62, 0.31, 'Photon particle gun')
t2_3.DrawLatex(.62, 0.25, 'Cell threshold 2 MeV')
t3.DrawLatex(0.15, 0.94, 'Background hits overlay in [-0.5, 10] ns range')
t4.DrawLatex(0.82, 0.94, '#sqrt{s} = 10 TeV')

c.SaveAs("photon_efficiency_vs_theta.pdf")

for entry in tree:
    h_resolution.Fill(entry.theta_truth, entry.E-entry.E_truth)

h_reso_theta = TH1D('reso_theta', 'reso_theta',
                    len(arrBins_theta)-1, arrBins_theta)
for bin in range(0, len(arrBins_theta)-1):
    h_my_proj = h_resolution.ProjectionY("_py", bin, bin+1)
    gaussFit = TF1("gaussfit", "gaus")
    h_my_proj.Fit(gaussFit, "E")
    sigma = gaussFit.GetParameter(2)
    sigma_err = gaussFit.GetParError(2)
    h_reso_theta.SetBinContent(bin+1, sigma/1000.)
    h_reso_theta.SetBinError(bin+1, sigma_err/1000.)

c2 = TCanvas("", "", 800, 600)
h_reso_theta.SetTitle(" ")
h_reso_theta.SetLineColor(kBlack)
h_reso_theta.SetLineWidth(2)
h_reso_theta.SetMaximum(0.035)
h_reso_theta.SetMinimum(0.015)
h_reso_theta.GetYaxis().SetTitle("Photon energy resolution #sigma_{E}/E")
h_reso_theta.GetYaxis().SetTitleOffset(1.4)
h_reso_theta.GetXaxis().SetTitleOffset(1.2)
# h_reso_theta.GetXaxis().SetRangeUser(-0.24, 0.5)
h_reso_theta.GetXaxis().SetTitle("Truth photon #theta [rad]")
h_reso_theta.Draw("E0")

leg = TLegend(.6, .75, .9, .85)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.035)
leg.AddEntry(h_reso_theta, "E_{#gamma} = 1 TeV", "L")
leg.Draw()

t3.DrawLatex(0.15, 0.94, 'Background hits overlay in [-0.5, 10] ns range')
t4.DrawLatex(0.82, 0.94, '#sqrt{s} = 10 TeV')

c2.SaveAs("reso_vs_theta.pdf")

fFile.Close()
