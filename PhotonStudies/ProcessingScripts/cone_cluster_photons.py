from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TH1D, TH2D, TH1, TFile, TLorentzVector, TMath, TTree, TVector3, TEfficiency
from math import *
from optparse import OptionParser
from array import array
import os
import fnmatch
import uproot
import mplhep as hep
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

#########################
parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile ntup_photons.root',
                  type=str, default='ntup_photons.root')
(options, args) = parser.parse_args()

arrBins_theta = np.linspace(0.175, TMath.Pi()-0.175, 30)#array('d', (0., 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            #90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., TMath.Pi()))
arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 50., 100., 250., 500., 1000.))#, 2500., 5000.))
arrBins_E_lowrange = array('d', (0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 75., 100., 125., 150., 200., 250., 300., 350., 400., 450., 500., 600., 700., 800., 900., 1000.))
arrBins_fakes = np.linspace(0,1000,250)
arrBins_dR = np.linspace(0,5.,100)

h_E_pass = TH1D('E_pass', 'E_pass', len(arrBins_E_lowrange)-1, arrBins_E_lowrange)
h_E_all = TH1D('E_all', 'E_all', len(arrBins_E_lowrange)-1, arrBins_E_lowrange)

h_theta_pass = TH1D('theta_pass', 'theta_pass', len(arrBins_theta)-1, arrBins_theta)
h_theta_all = TH1D('theta_all', 'theta_all', len(arrBins_theta)-1, arrBins_theta)

# declare histograms
h_truth_E = TH1D('truth_E', 'truth_E', len(arrBins_E)-1, arrBins_E)
h_truth_theta = TH1D('truth_theta', 'truth_theta', len(arrBins_theta)-1, arrBins_theta)

h_truth_pT = TH1D('truth_pT', 'truth_pT', len(arrBins_E)-1, arrBins_E) #GeV: true generator pT

h_reco_photon_E = TH1D('reco_photon_E', 'reco_photon_E', len(arrBins_E_lowrange)-1, arrBins_E_lowrange)

h_hit_dR = TH1D('dR', 'dR', len(arrBins_dR)-1, arrBins_dR)
h_min_hit_dR = TH1D('mindR', 'mindR', len(arrBins_dR)-1, arrBins_dR)
h_hit_r = TH1D('hit_r', 'hit_r', len(arrBins_dR)-1, arrBins_dR)
histos_list = [h_truth_E, h_truth_theta, h_truth_pT, h_reco_photon_E, h_hit_dR,
        h_min_hit_dR, h_hit_r,
        h_E_pass, h_E_all,
        h_theta_pass, h_theta_all]
# Histo list for writing to outputs

for histo in histos_list:
    histo.SetDirectory(0)

####################################
photon_tree = TTree("photon_tree", "photon_tree")
E = array('d', [0])
phi = array('d', [0])
theta = array('d', [0])
E_truth = array('d', [0])
phi_truth = array('d', [0])
theta_truth = array('d', [0])
photon_tree.Branch("E",  E,  'var/D')
photon_tree.Branch("phi", phi, 'var/D')
photon_tree.Branch("theta", theta, 'var/D')
photon_tree.Branch("E_truth",  E_truth,  'var/D')
photon_tree.Branch("theta_truth", theta_truth, 'var/D')
photon_tree.Branch("phi_truth", phi_truth, 'var/D')


to_process = []

if os.path.isdir(options.inFile):
    for r, d, f in os.walk(options.inFile):
        for file in f:
            to_process.append(os.path.join(r, file))
else:
    to_process.append(options.inFile)

#create the mesh for the dR contour plot
dthetas = np.linspace(0.,np.pi,800)
dphis = np.linspace(0., 2*np.pi, 800)
dRs = []

filenum=0
nevts_proc = 0
for file in to_process:
    # create a reader and open an LCIO file
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    try:
        reader.open(file)
    except Exception:
        #let it skip the bad files without breaking
        continue
    filenum=filenum+1
    #if filenum > 1:
    #    break
    # loop over all events in the file
    for ievt, event in enumerate(reader):
        #if filenum % 10 == 0:
        print(" ")
        print("File "+str(filenum))
        print("Processing event " + str(ievt))

        print(event.getCollectionNames())

        #get the truth particle, always the first in the collection
        mcpCollection = event.getCollection('MCParticle')
        truePhoton = mcpCollection[0]

        #get energy, momentum and fill truth histos
        trueE = truePhoton.getEnergy()
        h_truth_E.Fill(trueE)
        dp3 = truePhoton.getMomentum()
        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], trueE)
        h_truth_theta.Fill(tlv.Theta())

        #truth entries in photon tree
        E_truth[0] = trueE
        phi_truth[0] = abs(tlv.Phi())
        theta_truth[0] = tlv.Theta()
        

        #now get the ECal hits

        cone_photon_E = 0.
        hitphis = []
        hitthetas = []
        ecal_coll = ['HcalBarrelCollectionConed','HcalEndcapCollectionConed']
        #ecal_relcoll = ['EcalBarrelRelationsSimSel','EcalEndcapRelationsSimSel']
        min_hit_dR = 10.
        for icoll, coll in enumerate(ecal_coll):

            try:
                ECALhitCollection = event.getCollection(coll)
                #relationCollection = event.getCollection(ecal_relcoll[icoll])
                #relation = UTIL.LCRelationNavigator(relationCollection)

                encoding = ECALhitCollection.getParameters(
                ).getStringVal(EVENT.LCIO.CellIDEncoding)
                decoder = UTIL.BitField64(encoding)
                for hit in ECALhitCollection:
                    cellID = int(hit.getCellID0())
                    decoder.setValue(cellID)
                    layer = decoder["layer"].value()
                    #get the hit position and time
                    #determine if it is within dR = 0.05 from the truth particle
                    hit_time = hit.getTime()
                    hit_pos = hit.getPosition()
                    hit_vec = TLorentzVector(hit_pos[0], hit_pos[1], hit_pos[2], hit_time)
                    
                    h_hit_r.Fill(np.sqrt(hit_pos[0]**2+ hit_pos[1]**2+ hit_pos[2]**2))

                    
                    
                    hit_dR = tlv.DeltaR(hit_vec)
                    if hit_dR < min_hit_dR:
                        min_hit_dR = hit_dR
                    h_hit_dR.Fill(hit_dR)
                    if hit_dR < 0.05:
                        cone_photon_E += hit.getEnergy()
                        hitthetas.append(abs(hit_vec.Theta()))
                        hitphis.append(abs(hit_vec.Phi()))
            except: 
                print("No", coll, "found") 
        '''        
        
        dRs.clear()
        for itheta in dthetas:
            temp = []
            dtheta = theta_truth[0] - itheta
            for iphi in dphis:
                dphi = phi_truth[0] - iphi
                dR = np.sqrt(dtheta**2+dphi**2)
                if dR > 0.3:
                    temp.append(-.01)
                else:
                    temp.append(dR)
            dRs.append(temp)


         
        dRMax =.3
        fig, ax = plt.subplots()
        contr = ax.contourf(dphis, dthetas, dRs, cmap=plt.cm.Greys, vmin = 0.0, vmax = dRMax,levels=20,)
        cbar = fig.colorbar(contr)
        cbar.set_label('dR')
        ax.scatter(hitphis, hitthetas, color='springgreen', s = 1, label = "ECalHits")
        ax.scatter(phi_truth, theta_truth, color = 'magenta', s = 8, label = "Truth Photon")
        ax.set_xlabel("$\phi$", loc='right')
        ax.set_ylabel("$\\theta$", loc='top')
        ax.set_xlim(phi_truth[0]-np.pi/20.,phi_truth[0]+np.pi/20.)
        ax.set_ylim(theta_truth[0]-np.pi/20., theta_truth[0]+np.pi/20.)
        plt.grid(axis = 'both', linestyle = ":")
        ax.legend()
        fig.savefig("evtdisplay"+str(ievt)+".pdf")
        plt.close(fig)    
        '''     
        h_min_hit_dR.Fill(min_hit_dR)
        if cone_photon_E == 0.:
            h_reco_photon_E.Fill(-1.)
            E[0] = -1.
            theta[0] = -1.
            phi[0] = -1.
        else:
            h_reco_photon_E.Fill(cone_photon_E)
            E[0] = cone_photon_E
            theta[0] = np.average(hitthetas)
            phi[0] = np.average(hitphis)
            if theta_truth[0] > 0.175 and theta_truth[0] < 2.96:
                h_E_pass.Fill(E[0])
            if E_truth[0] > 100:
                h_theta_pass.Fill(theta[0])
        if theta_truth[0] > 0.175 and theta_truth[0] < 2.96:
            h_E_all.Fill(E[0])
        if E_truth[0] > 100:
            h_theta_all.Fill(theta[0])
        photon_tree.Fill()
    reader.close()
    
# write histograms
output_file = TFile(options.outFile, 'RECREATE')
for histo in histos_list:
    histo.Write()
photon_tree.Write()

output_file.Close()
