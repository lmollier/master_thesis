###FOR TTBAR BACKGROUND
####FOR EMU + TAU
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from helpers import delta_r, delta_phi, inv_mass_3p, cos_opening_angle, inv_mass, delta_eta
from samples import signal_samples
from helpers import files_from_dir, files_from_dirs
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import warnings
from scipy.stats import chisquare
from scipy.stats import ks_2samp
import os
from ROOT import TH1D



plt.style.use([hep.style.ROOT])
curr_os = os.name
if(curr_os == "nt"):
    local_dir = "C:\\Users\\lucas\\Desktop\\PDM\\data\\backgrounds\\tt2l2nu"
if(curr_os == "posix"):
    local_dir = "/media/sf_PDM/data/backgrounds/tt2l2nu"

warnings.filterwarnings("ignore")
lumi = 60000. # pb-1
xsec = 87.315 #pb
scale = lumi * xsec
samples = files_from_dir(local_dir),
i_file = 1
tot_events = 0
file = samples[0][0]
events = NanoEventsFactory.from_root(file, schemaclass=NanoAODSchema, entry_stop=10000).events()
events['SelElectron'] = events.Electron[(events.Electron.pt > 24.) & (events.Electron.mvaFall17V2Iso_WP90 > 0.5)]
events['SelMuon'] = events.Muon[(events.Muon.pt > 24.) & (events.Muon.mediumPromptId) & (events.Muon.pfRelIso03_all < 0.2) & (np.abs(events.Muon.dxy) < 0.005)]
muon1, electron1 = ak.unzip(ak.cartesian([events.SelMuon, events.SelElectron], nested=True))
match1 = ak.any(muon1.jetIdx == electron1.jetIdx, axis=-1, mask_identity=False)        
events['SelMuon'] = events.SelMuon[(~(match1))]    
events['SelTau'] = events.Tau[(events.Tau.pt > 20.) & (abs(events.Tau.eta) < 2.3)& (events.Tau.idDeepTau2017v2p1VSmu > 0.5) & (events.Tau.idDeepTau2017v2p1VSe > 0.5) & (events.Tau.idDeepTau2017v2p1VSjet >=8)]
tau2, electron2 = ak.unzip(ak.cartesian([events.SelTau, events.SelElectron], nested=True))
match2 = ak.any(tau2.jetIdx == electron2.jetIdx, axis=-1, mask_identity=False)
tau3, muon3 = ak.unzip(ak.cartesian([events.SelTau, events.SelMuon], nested=True))
match3 = ak.any(tau3.jetIdx == muon3.jetIdx, axis=-1, mask_identity=False)
events['SelTau'] = events.SelTau[((~(match2) & ~(match3)) )]  

events = events[ak.num(events.SelMuon) == 1]
events = events[ak.num(events.SelElectron) == 1]
jets = events.Jet[events.Jet.pt > 25.]
bjets = jets[jets.btagDeepFlavB > 0.2770]
one_bjet = ak.num(bjets) == 1
no_bjet = ak.num(bjets) == 0

CR_loose = events[one_bjet]
CR_medium = events[one_bjet]
CR_medium['SelTau'] = CR_medium.SelTau[CR_medium.SelTau.idDeepTau2017v2p1VSjet >=63]

CR_loose = CR_loose[ak.num(CR_loose.SelTau) == 1]
CR_medium = CR_medium[ak.num(CR_medium.SelTau) == 1]

CR_loose_kin = ak.flatten(CR_loose.SelTau.pt)
CR_medium_kin = ak.flatten(CR_medium.SelTau.pt)

bins_pt_tau = []
width = 3
for i in range(20,100):
    if((i) %width == 0):
        bins_pt_tau.append(i)
for i in range(100,220):
    if(i%(width*4) == 0):
        bins_pt_tau.append(i)  

hpx    = TH1D( 'hpx', 'This is the px distribution', 2, [1,2])


