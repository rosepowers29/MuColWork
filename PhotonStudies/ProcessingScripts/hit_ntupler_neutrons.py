from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TList, TH1D, TH2D, TH1, TFile, TLorentzVector, TMath, TTree, TVector3, TEfficiency
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
arrBins_dR = np.linspace(0,5.,100)



####################################
neutron_tree = TTree("neutron_tree", "neutron_tree")
E = array('d', [0])
phi = array('d', [0])
theta = array('d', [0])
eta = array('d', [0])
pdg_id = array('d', [0])
eta_truth = array('d', [0])
E_truth = array('d', [0])
phi_truth = array('d', [0])
theta_truth = array('d', [0])
pdg_truth = array('d', [0])
pT = array('d', [0])
pT_truth = array('d', [0])
N_neutron = array('d', [0])

truth_ProdVtx_x = array('d', [0])
truth_ProdVtx_y = array('d', [0])
truth_ProdVtx_z = array('d', [0])
ProdVtx_x = array('d', [0])
ProdVtx_y = array('d', [0])
ProdVtx_z = array('d', [0])
DecayVtx_x = array('d', [0])
DecayVtx_y = array('d', [0])
DecayVtx_z = array('d', [0])


Hhit_x = array('d')
Hhit_y = array('d', [0])
Hhit_z = array('d', [0])
Hhit_phi =array('d', [0])
Hhit_eta = array('d', [0])
Hhit_E = array('d', [0])
Hhit_depth = array('d', [0])
Hhit_time = array('d', [0])
##########################################################
neutron_tree.Branch("E",  E,  'var/D')
neutron_tree.Branch("phi", phi, 'var/D')
#neutron_tree.Branch("theta", theta, 'var/D')
neutron_tree.Branch("pT_truth", pT_truth, 'var/D')
neutron_tree.Branch("pT", pT, 'var/D')
neutron_tree.Branch("eta", eta, 'var/D')
neutron_tree.Branch("pdg_id", pdg_id, 'var/D')
neutron_tree.Branch("E_truth",  E_truth,  'var/D')
neutron_tree.Branch("phi_truth", phi_truth, 'var/D')
#neutron_tree.Branch("theta_truth", theta_truth, 'var/D')
neutron_tree.Branch("eta_truth", eta_truth, 'var/D')
neutron_tree.Branch("pdg_truth", pdg_truth, 'var/D')
neutron_tree.Branch("truth_ProdVtx_x", truth_ProdVtx_x, 'var/D')
neutron_tree.Branch("truth_ProdVtx_y", truth_ProdVtx_y, 'var/D')
neutron_tree.Branch("truth_ProdVtx_z", truth_ProdVtx_z, 'var/D')
neutron_tree.Branch("ProdVtx_x", ProdVtx_x, 'var/D')
neutron_tree.Branch("ProdVtx_y", ProdVtx_y, 'var/D')
neutron_tree.Branch("ProdVtx_z", ProdVtx_z, 'var/D')
neutron_tree.Branch("DecayVtx_x", DecayVtx_x, 'var/D')
neutron_tree.Branch("DecayVtx_y", DecayVtx_y, 'var/D')
neutron_tree.Branch("DecayVtx_z", DecayVtx_z, 'var/D')
#hit info
neutron_tree.Branch("Hhit_x.", Hhit_x, 'var/D',99)
neutron_tree.Branch("Hhit_y.", Hhit_y, 'var/D',99)
neutron_tree.Branch("Hhit_z.", Hhit_z, 'var/D',99)
neutron_tree.Branch("Hhit_phi.", Hhit_phi, 'var/D',99)
neutron_tree.Branch("Hhit_eta.", Hhit_eta, 'var/D',99)
neutron_tree.Branch("Hhit_E.", Hhit_E, 'var/D',99)
neutron_tree.Branch("Hhit_time.", Hhit_time, 'var/D',99)
neutron_tree.Branch("Hhit_depth.", Hhit_depth, 'var/D',99)

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
        trueNeutron = mcpCollection[0]

        #get energy, momentum and fill truth histos
        trueE = trueNeutron.getEnergy()
        dp3 = trueNeutron.getMomentum()
        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], trueE)

        #truth entries in neutron tree
        E_truth[0] = trueE
        phi_truth[0] = abs(tlv.Phi())
        eta_truth[0] = tlv.Eta()
                

        #now get the HCal hits

        cone_neutron_E = 0.
        hitphis = []
        hitthetas = []
        hcal_coll = ['HcalBarrelCollectionConed','HcalEndcapCollectionConed']
        #ecal_relcoll = ['EcalBarrelRelationsSimSel','EcalEndcapRelationsSimSel']
        min_hit_dR = 10.
        try:
            HCALhitCollection = event.getCollection('HcalBarrelCollectionConed')
        except:
            print("No Hcal Barrel collection")
        for icoll, coll in enumerate(hcal_coll):

            try:
                HCALhitCollection = event.getCollection(coll)

                encoding = HCALhitCollection.getParameters(
                ).getStringVal(EVENT.LCIO.CellIDEncoding)
                decoder = UTIL.BitField64(encoding)
                for hit in HCALhitCollection:
                    cellID = int(hit.getCellID0())
                    decoder.setValue(cellID)
                    layer = decoder["layer"].value()
                    #get the hit position and time
                    #determine if it is within dR = 0.05 from the truth particle
                    hit_time = hit.getTime()
                    hit_pos = hit.getPosition()
                    hit_vec = TLorentzVector(hit_pos[0], hit_pos[1], hit_pos[2], hit_time)
                    Hhit_x.append(hit_pos[0])

                    
                    
                    hit_dR = tlv.DeltaR(hit_vec)
                    if hit_dR < min_hit_dR:
                        min_hit_dR = hit_dR
                    if hit_dR < 0.05:
                        cone_neutron_E += hit.getEnergy()
                        hitthetas.append(abs(hit_vec.Theta()))
                        hitphis.append(abs(hit_vec.Phi()))
            
            except: 
                print("No", coll, "found") 
        if cone_neutron_E == 0.:
            E[0] = -1.
            theta[0] = -1.
            phi[0] = -1.
        else:
            E[0] = cone_neutron_E
            theta[0] = np.average(hitthetas)
            phi[0] = np.average(hitphis)
        neutron_tree.Fill()
    reader.close()
    
# write histograms
output_file = TFile(options.outFile, 'RECREATE')
neutron_tree.Write()

output_file.Close()
