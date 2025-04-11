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

#Global Flag
do_calo_tree = False

# declare histograms
h_ind_PFO_E = TH1D('ind_PFO_E', 'ind_PFO_E', len(arrBins_E)-1, arrBins_E)
h_truth_E = TH1D('truth_E', 'truth_E', len(arrBins_E)-1, arrBins_E)
h_truth_theta = TH1D('truth_theta', 'truth_theta', len(arrBins_theta)-1, arrBins_theta)
h_matched_theta = TH1D('matched_theta', 'matched_theta', len(arrBins_theta)-1, arrBins_theta)
h_matched_E = TH1D('matched_E', 'matched_E', len(arrBins_E)-1, arrBins_E)

h_truth_pT = TH1D('truth_pT', 'truth_pT', len(arrBins_E)-1, arrBins_E) #GeV: true generator pT

h_NFakes = TH1D('NFakes', 'NFakes', len(arrBins_fakes)-1, arrBins_fakes)

h_Npfo = TH1D('Npfo', "Npfo", 1000, 0, 1000)
h_EMpfo_E = TH1D('EMpfo_E', 'EMpfo_E', len(arrBins_E)-1, arrBins_E)
h_HADpfo_E = TH1D('HADpfo_E', 'HADpfo_E', len(arrBins_E)-1, arrBins_E)
h_matchedEMpfo_E = TH1D('matchedEMpfo_E', 'matchedEMpfo_E', len(arrBins_E)-1, arrBins_E)
h_matchedEMpfo_theta = TH1D('matchedEMpfo_theta', 'matchedEMpfo_theta',
                            len(arrBins_theta)-1, arrBins_theta)
h_matchedHADpfo_E = TH1D('matchedHADpfo_E', 'matchedHADpfo_E', len(arrBins_E)-1, arrBins_E)
h_matchedHADpfo_theta = TH1D('matchedHADpfo_theta', 'matchedHADpfo_theta',
                             len(arrBins_theta)-1, arrBins_theta)
h_pfo_type = TH1D('pfo_type', "pfo_type", 3000, 0, 3000)
h_mc_type = TH1D('mc_type', "mc_type", 3000, 0, 3000)
h_reco_photon_E = TH1D('reco_photon_E', 'reco_photon_E', len(arrBins_E_lowrange)-1, arrBins_E_lowrange)
h_jet_const_type = TH1D('jet_const_type', "jet_const_type", 3000, 0, 3000) #pfo object PDG number
h_pfo_pair_dR = TH1D('pfo_pair_dR', 'pfo_pair_dR', 100, 0, TMath.Pi())


h_deltaEM_E = TH1D('deltaEM_E', 'deltaEM_E', 250, -1000, 1000)
h_deltaHAD_E = TH1D('deltaHAD_E', 'deltaHAD_E', 250, -1000, 1000)
h_delta_E_sumE = TH1D('delta_E_sumE', 'delta_E_sumE', 250, -1000, 1000)

# Low-level sim hit distributions
h_ECAL_simhit_E = TH1D('ECAL_simhit_E', 'ECAL_simhit_E', 100, 0, 20)  # GeV
h_ECAL_simhit_layer = TH1D(
    'ECAL_simhit_layer', 'ECAL_simhit_layer', 100, 0, 100)
h_ECAL_simhit_layer_ele = TH1D(
    'ECAL_simhit_layer_ele', 'ECAL_simhit_layer_ele', 100, 0, 100)
h_ECAL_simhit_layer_gamma = TH1D(
    'ECAL_simhit_layer_gamma', 'ECAL_simhit_layer_gamma', 100, 0, 100)
h_ECAL_simhit_layer_other = TH1D(
    'ECAL_simhit_layer_other', 'ECAL_simhit_layer_other', 100, 0, 100)

h_HCAL_simhit_E = TH1D('HCAL_simhit_E', 'HCAL_simhit_E', 100, 0, 20)  # GeV
h_HCAL_simhit_layer = TH1D(
    'HCAL_simhit_layer', 'HCAL_simhit_layer', 100, 0, 100)

# Low-level digitised hit distributions
h_ECAL_nhits = TH1D('ECAL_nhits', 'ECAL_nhits', 1000, 0, 2000)
h_ECAL_hit_time = TH1D('ECAL_hit_time', 'ECAL_hit_time', 100, -10, 10)  # ns
h_ECAL_hit_E = TH1D('ECAL_hit_E', 'ECAL_hit_E', 100, 0, 20)  # GeV
h_ECAL_hit_R = TH1D('ECAL_hit_R', 'ECAL_hit_R', 100, 1700, 4000)  # m
h_ECAL_hit_layer = TH1D('ECAL_hit_layer', 'ECAL_hit_layer', 100, 0, 100)

h_HCAL_hit_time = TH1D('HCAL_hit_time', 'HCAL_hit_time', 100, -10, 10)  # ns
h_HCAL_hit_E = TH1D('HCAL_hit_E', 'HCAL_hit_E', 100, 0, 20)  # GeV
h_HCAL_hit_R = TH1D('HCAL_hit_R', 'HCAL_hit_R', 100, 1700, 4000)  # m
h_HCAL_hit_layer = TH1D('HCAL_hit_layer', 'HCAL_hit_layer', 100, 0, 100)

# Aggregated energy info
h_sumE = TH1D('sumE', 'sumE', 120, 0, 6000)  # GeV
h_ECAL_sumE = TH1D('ECAL_sumE', 'ECAL_sumE', 120, 0, 6000)  # GeV
h_HCAL_sumE = TH1D('HCAL_sumE', 'HCAL_sumE', 120, 0, 6000)  # GeV
h_EMfrac = TH1D('EMfrac', 'EMfrac', 100, 0, 1)  # GeV
h_EMfrac_PFO = TH1D('EMfrac_PFO', 'EMfrac_PFO', 100, 0, 1)  # GeV

### Jet Histos
h_jet_energy = TH1D('jet_energy', 'jet_energy', len(arrBins_E)-1, arrBins_E)
h_jet_dR = TH1D('jet_dR', 'jet_dR', 100, 0, TMath.Pi())
h_min_jet_dR = TH1D('min_jet_dR', 'min_jet_dR', 100, 0, TMath.Pi())
h_jet_multiplicity = TH1D('jet_multiplicity', 'jet_multiplicity', 100, 0, 150)
h_dR_vs_E = TH2D ('dR_vs_E', 'dR_vs_E', 100, 0, TMath.Pi(), len(arrBins_E)-1, arrBins_E)

h_E_pass = TH1D('E_pass', 'E_pass', len(arrBins_E_lowrange)-1, arrBins_E_lowrange)
h_E_all = TH1D('E_all', 'E_all', len(arrBins_E_lowrange)-1, arrBins_E_lowrange)

h_theta_pass = TH1D('theta_pass', 'theta_pass', len(arrBins_theta)-1, arrBins_theta)
h_theta_all = TH1D('theta_all', 'theta_all', len(arrBins_theta)-1, arrBins_theta)


h_jet_dR_pairs = TH1D('jet_dR_pairs', 'jet_dR_pairs', 25, 0, 15)
h_jet_hadFrac = TH1D('jet_had_frac', 'jet_had_frac', 50,0,1)

h_dR_vs_hadFrac = TH2D('dR_vs_hadFrac', 'dR_vs_hadFrac', 100, 0, TMath.Pi(), 50, 0, 1)

# Histo list for writing to outputs
histos_list = [h_ind_PFO_E, h_NFakes, h_ind_PFO_E, h_truth_E, h_truth_theta, h_matched_theta,
               h_EMpfo_E, h_matched_E,
               h_HADpfo_E,
               h_matchedEMpfo_E, h_matchedEMpfo_theta,
               h_matchedHADpfo_E, h_matchedHADpfo_theta,
               h_deltaEM_E, h_deltaHAD_E, h_delta_E_sumE,
               h_reco_photon_E,
               h_Npfo, h_pfo_type, h_mc_type,
               h_pfo_pair_dR,
               h_jet_const_type, h_truth_pT,
               h_ECAL_hit_time, h_ECAL_hit_E, h_ECAL_hit_R,
               h_ECAL_nhits,
               h_HCAL_hit_time, h_HCAL_hit_E, h_HCAL_hit_R,
               h_ECAL_simhit_E, h_HCAL_simhit_E,
               h_ECAL_sumE, h_HCAL_sumE,
               h_EMfrac, h_EMfrac_PFO,
               h_ECAL_hit_layer, h_HCAL_hit_layer,
               h_ECAL_simhit_layer, h_ECAL_simhit_layer_ele, h_ECAL_simhit_layer_gamma, h_ECAL_simhit_layer_other,
               h_HCAL_simhit_layer, 
               h_jet_energy, h_jet_dR, h_min_jet_dR, h_jet_multiplicity,
               h_E_pass, h_E_all,
               h_theta_pass, h_theta_all,
               h_jet_dR_pairs,
               h_dR_vs_E, h_jet_hadFrac, h_dR_vs_hadFrac
               #h_NFakes
               ]

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
pT = array('d', [0])
pT_truth = array('d', [0])
N_photon = array('d', [0])
photon_tree.Branch("E",  E,  'var/D')
photon_tree.Branch("phi", phi, 'var/D')
photon_tree.Branch("theta", theta, 'var/D')
photon_tree.Branch("E_truth",  E_truth,  'var/D')
photon_tree.Branch("phi_truth", phi_truth, 'var/D')
photon_tree.Branch("theta_truth", theta_truth, 'var/D')
photon_tree.Branch("N_photon", N_photon, 'var/D')


####################################
calo_tree = TTree("calo_tree", "calo_tree")
E_cell = array('d', [0])
#theta_cell = array('d', [0])
E_sim = array('d', [0])
#theta_sim = array('d', [0])
isECAL = array('i', [0])
calo_tree.Branch("E_cell",  E_cell,  'var/D')
#calo_tree.Branch("theta_cell", theta_cell, 'var/D')
calo_tree.Branch("E_sim",  E_sim,  'var/D')
#calo_tree.Branch("theta_sim", theta_sim, 'var/D')
calo_tree.Branch("isECAL",  isECAL,  'var/I')

##############jet tree
jet_tree = TTree("jet_tree", "jet_tree")
photon_E = array('d', [0]) #initial neutron energy
photon_pT = array('d', [0]) #initial neutron pT
photon_phi_angle = array('d', [0]) #initial neutron phi
photon_theta_angle = array('d', [0]) #initial neutron theta
#neutron_index = array('d', [0]) #initial neutron indexing number
jet_energy = array('d', [0]) #reconstructed jet energy
jet_pT = array('d', [0]) #reconstructed jet pT
jet_phi = array('d', [0]) #reconstructed jet phi
jet_theta = array('d', [0]) #reconstructed jet theta
jet_dR = array('d', [0]) #angular distance between jet and neutron

jet_tree.Branch("photon_E", photon_E, 'var/D')
jet_tree.Branch("photon_pT", photon_pT, 'var/D')
jet_tree.Branch("photon_phi_angle", photon_phi_angle, 'var/D')
jet_tree.Branch("photon_theta_angle", photon_theta_angle, 'var/D')
#jet_tree.Branch("neutron_index", neutron_index, 'var/I')
jet_tree.Branch("jet_energy", jet_energy, 'var/D')
jet_tree.Branch("jet_pT", jet_pT, 'var/D')
jet_tree.Branch("jet_phi", jet_phi, 'var/D')
jet_tree.Branch("jet_theta", jet_theta, 'var/D')
jet_tree.Branch("jet_dR", jet_dR, 'var/D')

##############PFO tree (PID agnostic)
PFO_tree = TTree("PFO_tree", "PFO_tree")
MCP_E = array('d', [0]) #initial  energy
MCP_pT = array('d', [0]) #initial  pT
MCP_phi_angle = array('d', [0]) #initial phi
MCP_theta_angle = array('d', [0]) #initial  theta
#neutron_index = array('d', [0]) #initial neutron indexing number
PFO_energy = array('d', [0]) #reconstructed  energy
PFO_pT = array('d', [0]) #reconstructed pT
PFO_phi = array('d', [0]) #reconstructed phi
PFO_theta = array('d', [0]) #reconstructed theta
PFO_dR = array('d', [0]) #angular distance between PFO and MCP
PFO_hadFrac = array('d', [0]) #hadronic energy fraction of PFO (from clusters)

#PFO_tree.Branch("MCP_E", MCP_E, 'var/D')
#PFO_tree.Branch("MCP_pT", MCP_pT, 'var/D')
#PFO_tree.Branch("MCP_phi_angle", MCP_phi_angle, 'var/D')
#PFO_tree.Branch("MCP_theta_angle", MCP_theta_angle, 'var/D')
#jet_tree.Branch("neutron_index", neutron_index, 'var/I')
PFO_tree.Branch("PFO_energy", PFO_energy, 'var/D')
PFO_tree.Branch("PFO_pT", PFO_pT, 'var/D')
PFO_tree.Branch("PFO_phi", PFO_phi, 'var/D')
PFO_tree.Branch("PFO_theta", PFO_theta, 'var/D')
PFO_tree.Branch("PFO_dR", PFO_dR, 'var/D')
PFO_tree.Branch("PFO_hadFrac", PFO_hadFrac, 'var/D')


to_process = []

if os.path.isdir(options.inFile):
    for r, d, f in os.walk(options.inFile):
        for file in f:
            to_process.append(os.path.join(r, file))
else:
    to_process.append(options.inFile)

filenum=0
nevts_proc = 0
for file in to_process:
    # create a reader and open an LCIO file
    #if filenum > 1000:
    #    break
    #if nevts_proc > 100:
    #    break
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    try:
        reader.open(file)
    except Exception:
        #let it skip the bad files without breaking
        continue
    filenum=filenum+1
    # loop over all events in the file
    for ievt, event in enumerate(reader):
        if ievt % 1 == 0:
            print(" ")
            print("File "+str(filenum))
            print("Processing event " + str(ievt))

        print(event.getCollectionNames())

        # Fill the truth-level histos, the first particle is always the gun
        mcpCollection = event.getCollection('MCParticle')
        #for mcp in mcpCollection:
        #    h_mc_type.Fill(abs(mcp.getPDG()))
        h_truth_E.Fill(mcpCollection[0].getEnergy())
        dp3 = mcpCollection[0].getMomentum()
        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcpCollection[0].getEnergy())
        h_truth_theta.Fill(tlv.Theta())

        E_truth[0] = mcpCollection[0].getEnergy()
        phi_truth[0] = tlv.Phi()
        theta_truth[0] = tlv.Theta()
        MCP_E[0]=E_truth[0]
        MCP_phi_angle[0]=phi_truth[0]
        MCP_theta_angle[0]=theta_truth[0]
        MCP_pT[0]=tlv.Perp()

        #h_E_all.Fill(E_truth[0])

        '''
        for mcp in mcpCollection:
            if mcp.getGeneratorStatus() == 1 and len(mcp.getParents()) == 0:
                dp3 = mcp.getMomentum()
                tlv = TLorentzVector()
                tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcp.getEnergy())
                h_truth_E.Fill(mcp.getEnergy())
                h_truth_theta.Fill(tlv.Theta())
        '''
        
        # Fill the reco-level histos
        pfoCollection = event.getCollection('PandoraPFOs')
        h_Npfo.Fill(len(pfoCollection))
        #restrict to high-multiplicity events
        #if len(pfoCollection) < 100:
        #    print("skipping")
        #    continue
        nevts_proc += 1
        #print(nevts_proc)
        # Match true pfo with closest reco PFO in deltaR
        # Or with highest pT (better for BIB)
        matchedEM_E = -1.
        matchedEM_theta = -1.
        matchedHAD_E = -1.
        matchedHAD_theta = -1.
        allEM_E = 0.
        allHAD_E = 0.

        minDREM = 0.1
        minDRHAD = 0.5
        maxHadFrac = 0.1
        maxpT = 0.0

        E[0] = -1
        phi[0] = -4
        theta[0] = -1
        
        EcalE=0.
        HcalE=0.
        HcalFrac = 0.
        PFO_hadFrac[0]=0.
        #noncentral = True
        #initialize NFakes, we will subtract from it as we exclude on energy basis
        NFakes = len(pfoCollection)-1
        NPhoton = 0
        n_ind=0
        for pfo in pfoCollection:
            h_pfo_type.Fill(abs(pfo.getType()))
            dp3 = pfo.getMomentum()
            tlv_pfo = TLorentzVector()
            tlv_pfo.SetPxPyPzE(dp3[0], dp3[1], dp3[2], pfo.getEnergy())
            n_pair_ind = -1
            if filenum == 1 and ievt == 2:
                h_ind_PFO_E.Fill(pfo.getEnergy())
            '''
            for pairPFO in pfoCollection:
                n_pair_ind = n_ind+1
                if n_pair_ind < n_ind:
                    continue
                else:
                    pairdp3 = pairPFO.getMomentum()
                    tlv_pairPFO = TLorentzVector()
                    tlv_pairPFO.SetPxPyPzE(pairdp3[0], pairdp3[1], pairdp3[2], pairPFO.getEnergy())
                    pairdR = tlv_pairPFO.DeltaR(tlv_pfo)
                    h_pfo_pair_dR.Fill(pairdR)
            n_ind = n_ind + 1
            '''
            pT_pfo = tlv_pfo.Perp()
            
            if abs(pfo.getType()) == 22:
                allEM_E = allEM_E + pfo.getEnergy()
                NPhoton = NPhoton+1
                h_reco_photon_E.Fill(pfo.getEnergy())
            elif abs(pfo.getType()) == 2112:
                allHAD_E = allHAD_E + pfo.getEnergy()
                #don't include neutrons in fakes
                NFakes = NFakes-1
                
            #fill N_photon with the number of non-neutron PFO    
            N_photon[0]=NPhoton
            dR = tlv_pfo.DeltaR(tlv)
            if pfo.getEnergy() < 20.:
                #subtract away from NFakes
                NFakes = NFakes - 1
                continue
            if mcpCollection[0].getEnergy() < 10.:
                continue
            
            if pfo.getEnergy()/mcpCollection[0].getEnergy() > 0.: #change threshold here if needed
                if pT_pfo > maxpT and abs(pfo.getType()) == 22:
                    if dR<minDREM:
                        maxpT = pT_pfo
                        minDREM = dR
                        matchedEM_E = pfo.getEnergy()
                        matchedEM_theta = tlv_pfo.Theta()
                        E[0] = pfo.getEnergy()
                        phi[0] = tlv_pfo.Phi()
                        theta[0] = tlv_pfo.Theta()
                        #if abs(theta[0]) < TMath.Pi()/6:
                        #   print("central")
                if dR < minDRHAD and abs(pfo.getType()) == 2112:
                    minDRHAD = dR
                    matchedHAD_E = pfo.getEnergy()
                    matchedHAD_theta = tlv_pfo.Theta()
            

        print(allEM_E, allHAD_E, matchedEM_E, matchedHAD_E)
        #reader.close()

        h_EMpfo_E.Fill(allEM_E)
        h_HADpfo_E.Fill(allHAD_E)
        
        #h_E_all.Fill(mcpCollection[0].getEnergy())

        if matchedEM_E > 0:
            h_matchedEMpfo_E.Fill(matchedEM_E)
            h_matchedEMpfo_theta.Fill(matchedEM_theta)
            h_deltaEM_E.Fill(matchedEM_E-mcpCollection[0].getEnergy())
            #cut out edges of acceptance
            if tlv.Theta() > 0.175 and tlv.Theta() < 2.967:
                h_E_pass.Fill(mcpCollection[0].getEnergy())
            if mcpCollection[0].getEnergy() > 100:
                h_theta_pass.Fill(tlv.Theta())
            
        else:
            print("No match")
        if matchedHAD_E > 0:
            h_matchedHADpfo_E.Fill(matchedHAD_E)
            h_matchedHADpfo_theta.Fill(matchedHAD_theta)
            h_deltaHAD_E.Fill(matchedHAD_E-mcpCollection[0].getEnergy())

        if allHAD_E+allEM_E > 0:
            h_EMfrac_PFO.Fill(allEM_E/(allHAD_E+allEM_E))
        print(NFakes) 
        h_NFakes.Fill(NFakes)
        photon_tree.Fill()

        PFO_energy[0]=E[0]
        PFO_phi[0]=phi[0]
        PFO_theta[0]=theta[0]
        if tlv.Theta() > 0.175 and tlv.Theta() < 2.967:
            h_E_all.Fill(mcpCollection[0].getEnergy())
        if mcpCollection[0].getEnergy() > 100:
            h_theta_all.Fill(tlv.Theta())
        PFO_tree.Fill()
        reader.close()
        '''
        try:
            jetCollection = event.getCollection('JetOut')
        except Exception:
            continue

        minDRjet = 9999999.
        maxHadFracjet = 0.1
        jet_hadFrac = 1.0
        matchedJetE = -1.

        E[0] = -1
        pT[0] = -1
        phi[0] = -4
        theta[0] = -1
        #final_minDRjet[0] = -TMath.Pi()

        jet_energy[0] = -1
        jet_dR[0] = -1
        
        print(len(jetCollection), "jets")
        #jet_multiplicity[0] = len(jetCollection)
        h_jet_multiplicity.Fill(len(jetCollection))
        n_ind = -1
        for jet in jetCollection:
            n_ind = n_ind+1
            dp3 = jet.getMomentum()
            tlv_pfo = TLorentzVector()
            tlv_pfo.SetPxPyPzE(dp3[0], dp3[1], dp3[2], jet.getEnergy())

            dR = tlv_pfo.DeltaR(tlv)
            #clusters=jet.getClusters()
            #print(clusters.size())

        '''
        '''
            n_pair_ind = -1
            for pairJet in jetCollection:
                n_pair_ind = n_pair_ind+1
                if n_pair_ind < n_ind:
                    continue
                else:
                    pairdp3 = pairJet.getMomentum()
                    tlv_pairjet = TLorentzVector()
                    tlv_pairjet.SetPxPyPzE(pairdp3[0], pairdp3[1], pairdp3[2], pairJet.getEnergy())
                    pairdR = tlv_pairjet.DeltaR(tlv_pfo)
                    h_jet_dR_pairs.Fill(pairdR)
        '''
        '''    
            
            h_jet_energy.Fill(jet.getEnergy())
            h_jet_dR.Fill(dR)
            h_dR_vs_E.Fill(dR, jet.getEnergy())
            jet_theta[0] = tlv_pfo.Theta()
            jet_phi[0] = tlv_pfo.Phi()
            photon_E[0] = E_truth[0]
            photon_pT[0] = pT_truth[0]
            photon_phi_angle[0] = phi_truth[0]
            photon_theta_angle[0] = theta_truth[0]
            #neutron_index[0] = n_ind
            #jet_energy[0] = jet.getEnergy()
            #jet_pT[0] = tlv_pfo.Perp()
            #jet_phi[0] = tlv_pfo.Phi()
            #jet_theta[0] = tlv_pfo.Theta()
            #jet_dR[0] = dR
            
        ''' 
        ''' 
                if jet_n > 0:
                #print("a jet is found!")
                if jet_energy[0] == -1:
                    jet_energy[0] = jet.getEnergy()
                else:
                    jet_energy.append(jet.getEnergy())
                jet_dR.append(dR)
        '''
        '''   
            if jet_energy[0] == -1:
                jet_energy[0] = jet.getEnergy()
                #print("1st jet found: jet_energy = {:3e}".format(jet_energy[0]))
            #else:
                #jet_energy.append(jet.getEnergy())
                #print("2nd jet found: jet_energy = {:3e}".format(jet_energy[-1]))

            if jet_dR[0] == -1:
                jet_dR[0] = dR
            else:
                jet_dR.append(dR)
        
        
            #jet_hadFrac=1.0
            #EcalE = 0.0
            #HcalE = 0.0
            matchedJetE = -1.0

            if dR < minDRjet:
                #minDRjet = dR
                #print(dR)
                #get the hadronic content of each of the constituent PF's
                #if we are above the threshold, discard the jet
                #below is not necessary for no-BIB samples
        '''
        '''
                PFs = jet.getParticles()
                #print("PFs:",len(PFs))
                EcalE = 0.
                HcalE = 0.
                for particle in PFs:
                    clusters = particle.getClusters()
                    #print(len(clusters))
                    for cluster in clusters:
                        subEnergies = cluster.getSubdetectorEnergies()
                        #[0] is ecal energy, [1] is hcal
                        EcalE = EcalE + subEnergies[0]
                        HcalE = HcalE + subEnergies[1]
                #now get the hadronic fraction
                if (HcalE+EcalE) > 0:
                    jet_hadFrac = HcalE / (HcalE + EcalE)
                else: 
                    jet_hadFrac = 1.0
                
                
                if jet_hadFrac > maxHadFracjet:
                    #print("hadFrac:",jet_hadFrac)
                    continue
                else:
                    maxHadFracjet = jet_hadFrac
        '''
        '''
                minDRjet = dR
                #h_jet_hadFrac.Fill(jet_hadFrac)
                matchedJetE=jet.getEnergy()
                pT[0] = tlv.Perp()
                phi[0] = tlv_pfo.Phi()
                theta[0] = tlv_pfo.Theta()
        
        jet_energy[0] = matchedJetE
        jet_tree.Fill()
        '''
        '''
        #final_minDRjet[0] = minDRjet
        h_min_jet_dR.Fill(minDRjet)
        #print("HadFrac_Jet = ", jet_hadFrac)
        h_E_all.Fill(mcpCollection[0].getEnergy())
        if minDRjet<1000 and jet_hadFrac < 0.1:
            h_matched_theta.Fill(tlv.Theta())
            h_matched_E.Fill(tlv.Energy())
            h_E_pass.Fill(mcpCollection[0].getEnergy())
        else:
            print("dR:",minDRjet, "HF:", jet_hadFrac)
        
        clusterCollection = event.getCollection('PandoraClusters')
        minDRclus = 99999.
        for cluster in clusterCollection:

            px = cluster.getEnergy()*sin(cluster.getITheta())*cos(cluster.getIPhi())
            py = cluster.getEnergy()*sin(cluster.getITheta())*sin(cluster.getIPhi())
            pz = cluster.getEnergy()*cos(cluster.getITheta())

            tlv_clus = TLorentzVector()
            tlv_clus.SetPxPyPzE(px, py, pz, cluster.getEnergy())

            dR = tlv_clus.DeltaR(tlv)

            if dR < minDRclus:
                minDRclus = dR
                
                E[0] = cluster.getEnergy()
                phi[0] = cluster.getIPhi()
                theta[0] = cluster.getITheta()
                


        #if noncentral:
        photon_tree.Fill()

        for jet in jetCollection:
            dp3 = jet.getMomentum()
            tlv_pfo = TLorentzVector()
            tlv_pfo.SetPxPyPzE(dp3[0], dp3[1], dp3[2], jet.getEnergy())

            dR = tlv_pfo.DeltaR(tlv)

            if dR < minDRjet:
                minDRjet = dR
        
        if minDRjet<0.4:
            h_matched_theta.Fill(tlv.Theta())
            h_matched_E.Fill(tlv.Energy())
        
        # Fill the hit-level histos and aggregated energy
    
        ECAL_sumE = 0.
        ECAL_sumhits = 0.
        
        try:
            ECALhits_endcap = event.getCollection(EcalEndcapCollectionSel)
            ECAL_sumhits += len(ECALhits_endcap)
        except:
            print("No EcalEndcapCollectionSel in event")

        try:
            ECALhits_barrel = event.getCollection(EcalBarrelCollectionSel)
            ECAL_sumhits += len(ECALhits_barrel)
        except:
            print("No EcalBarrelCollectionSel in event")
        
        ecal_coll = ['EcalBarrelCollectionSel','EcalEndcapCollectionSel']
        #ecal_relcoll = ['EcalBarrelRelationsSimSel','EcalEndcapRelationsSimSel']
        
        for icoll, coll in enumerate(ecal_coll):

            try:
                ECALhitCollection = event.getCollection(coll)
         
                #relationCollection = event.getCollection(ecal_relcoll[icoll])
                #relation = UTIL.LCRelationNavigator(relationCollection)

                encoding = ECALhitCollection.getParameters(
                ).getStringVal(EVENT.LCIO.CellIDEncoding)
                decoder = UTIL.BitField64(encoding)
                
                ECAL_sumhits += len(ECALhitCollection)

                for hit in ECALhitCollection:
                    cellID = int(hit.getCellID0())
                    decoder.setValue(cellID)
                    layer = decoder["layer"].value()
                    h_ECAL_hit_time.Fill(hit.getTime())
                    h_ECAL_hit_E.Fill(hit.getEnergy())
                    h_ECAL_hit_layer.Fill(layer, hit.getEnergy())

                    ECAL_sumE = ECAL_sumE + hit.getEnergy()
                    pos = hit.getPosition()
                    h_ECAL_hit_R.Fill(sqrt(pos[0]*pos[0]+pos[1]*pos[1]))

                    E_cell[0] = hit.getEnergy()
                    #sim_hit = relation.getRelatedToObjects(hit)[0]
                    #E_sim[0] = sim_hit.getEnergy()
                    isECAL[0] = 1

                    calo_tree.Fill()

                h_ECAL_sumE.Fill(ECAL_sumE)
            except:
                print("No", coll, "found")
        


        h_ECAL_nhits.Fill(ECAL_sumhits)
        HCAL_sumE = 0.
        hcal_coll = ['HcalBarrelCollectionConed','HcalEndcapCollectionConed']
        hcal_relcoll = ['HcalBarrelRelationsSimSel','HcalEndcapRelationsSimSel']
        
        for icoll, coll in enumerate(hcal_coll):

            try:
                HCALhitCollection = event.getCollection(coll)

                #relationCollection = event.getCollection(hcal_relcoll[icoll])
                #relation = UTIL.LCRelationNavigator(relationCollection)

                encoding = HCALhitCollection.getParameters(
                ).getStringVal(EVENT.LCIO.CellIDEncoding)
                decoder = UTIL.BitField64(encoding)

                for hit in HCALhitCollection:
                    cellID = int(hit.getCellID0())
                    decoder.setValue(cellID)
                    layer = decoder["layer"].value()

                    h_HCAL_hit_time.Fill(hit.getTime())
                    h_HCAL_hit_E.Fill(hit.getEnergy())
                    h_HCAL_hit_layer.Fill(layer, hit.getEnergy())
                    HCAL_sumE = HCAL_sumE + hit.getEnergy()
                    pos = hit.getPosition()
                    h_HCAL_hit_R.Fill(sqrt(pos[0]*pos[0]+pos[1]*pos[1]))

                    E_cell[0] = hit.getEnergy()
                    sim_hit = relation.getRelatedToObjects(hit)[0]
                    E_sim[0] = sim_hit.getEnergy()
                    isECAL[0] = 0

                    calo_tree.Fill()

                h_HCAL_sumE.Fill(HCAL_sumE)
            except:
                print("No", coll, "found")

        print(ECAL_sumE, HCAL_sumE)

        h_sumE.Fill(ECAL_sumE+HCAL_sumE)

        if ECAL_sumE+HCAL_sumE > 0:
            h_EMfrac.Fill(ECAL_sumE/(ECAL_sumE+HCAL_sumE))
        else:
            h_EMfrac.Fill(0)
        h_delta_E_sumE.Fill(ECAL_sumE+HCAL_sumE-mcpCollection[0].getEnergy())
    '''
    reader.close()
    
# write histograms
output_file = TFile(options.outFile, 'RECREATE')
for histo in histos_list:
    histo.Write()
photon_tree.Write()
#jet_tree.Write()
#PFO_tree.Write()
if do_calo_tree:
    calo_tree.Write()

'''
histofile = uproot.open(options.outFile)
print(histofile.classnames())
pass_histo = histofile['E_pass;1']
all_histo = histofile['E_all;1']

#TEfficiency
upErrors=array('d')
lowErrors=array('d')
binEff=array('d')
photonEff = TEfficiency(pass_histo, all_histo)
for Ebin in range(len(arrBins_E)-1):
    upErrors[Ebin]=photonEff.GetEfficiencyErrorUp()
    lowErrors[Ebin]=photonEff.GetEfficiencyErrorLow()
    binEff[Ebin]=photonEff.GetEfficiency()
#plot using MAIA conventions
plt.errorbar(arrBins_E, binEff, xerr = (arrBins_E[1]-arrBins_E[0])/2, yerr = [lowErrors, upErrors], fmt='.', color='red',label="No HadFrac cutoff, dR <= 0.5")
plt.xlabel("True MC Photon Energy [GeV]")
plt.ylabel("Photon Reco Efficiency")
hep.cms.label(exp = "Muon Collider", data = False, rlabel = "MAIA Detector Concept", loc=0, italic=(1,0,0))
plt.legend()
plt.show()
plt.savefig("gammaEff_maia.pdf")
plt.close()
'''
output_file.Close()
