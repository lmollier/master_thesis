###FOR TTBAR BACKGROUND
####FOR Etau + TAU
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from coffea import hist
from helpers import delta_r, delta_phi, inv_mass_3p, cos_opening_angle, inv_mass, delta_eta
from samples import signal_samples
from helpers import files_from_dir, files_from_dirs
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from ROOT import TH1D, TFile
import warnings
from scipy.stats import chisquare
from scipy.stats import ks_2samp
import os
plt.style.use([hep.style.CMS])
curr_os = os.name
local_dir = "/eos/user/l/lmollier/PDM/data/bck/nanoAOD/tt2l2nu"
print("##########################################\n##########################################")
print(local_dir)
nb_of_e = 0
nb_of_tau = 2
nb_of_mu = 1
workingpoint = 15
path_root_file_histo = "/afs/cern.ch/user/l/lmollier/PDM/results_histograms/histos_dileptons_fake_rate.root"
f = TFile(path_root_file_histo ,"RECREATE")
warnings.filterwarnings("ignore")
lumi = 60000. # pb-1
xsec = 87.315 #pb
scale = lumi * xsec
samples = files_from_dir(local_dir),
i_file = 1
tot_events = 0
file = samples[0][0]
print("PF and stuff like that")
events = NanoEventsFactory.from_root(file, schemaclass=NanoAODSchema).events()
print("events are charged!")
events['SelElectron'] = events.Electron[(events.Electron.pt > 24.) & (events.Electron.mvaFall17V2Iso_WP90 > 0.5)]
events['SelMuon'] = events.Muon[(events.Muon.pt > 24.) & (events.Muon.mediumPromptId) & (events.Muon.pfRelIso03_all < 0.2) & (np.abs(events.Muon.dxy) < 0.005)]
muon1, electron1 = ak.unzip(ak.cartesian([events.SelMuon, events.SelElectron], nested=True))
match1 = ak.any(muon1.jetIdx == electron1.jetIdx, axis=-1, mask_identity=False)        
events['SelMuon'] = events.SelMuon[(~(match1))]    
events['SelTau'] = events.Tau[(events.Tau.pt > 20.) & (abs(events.Tau.eta) < 2.3)& (events.Tau.idDeepTau2017v2p1VSmu > 0.5) & (events.Tau.idDeepTau2017v2p1VSe > 0.5) & (events.Tau.idDeepTau2017v2p1VSjet >=7)]
tau2, electron2 = ak.unzip(ak.cartesian([events.SelTau, events.SelElectron], nested=True))
match2 = ak.any(tau2.jetIdx == electron2.jetIdx, axis=-1, mask_identity=False)
tau3, muon3 = ak.unzip(ak.cartesian([events.SelTau, events.SelMuon], nested=True))
match3 = ak.any(tau3.jetIdx == muon3.jetIdx, axis=-1, mask_identity=False)
events['SelTau'] = events.SelTau[((~(match2) & ~(match3)) )]  
print("selection of leptons done!")
events = events[ak.num(events.SelMuon) +ak.num(events.SelElectron) == nb_of_mu + nb_of_e]
#events = events[ak.num(events.SelElectron) == nb_of_e]
events = events[ak.num(events.SelTau) == nb_of_tau]
print(str(nb_of_e + nb_of_mu)+" electron or muon, and "+str(nb_of_tau)+" taus required!")

events['LeadTau'] = events.SelTau[:,0]
events['SubTau'] = events.SelTau[:,1]
temp_2_passed_1 = events.LeadTau.idDeepTau2017v2p1VSjet >=workingpoint 
temp_2_passed_2 =  events.SubTau.idDeepTau2017v2p1VSjet >=workingpoint
temp_2_passed = temp_2_passed_1 & temp_2_passed_2
events_2_passed = events[temp_2_passed]
temp_2_failed_1 = ~(events.LeadTau.idDeepTau2017v2p1VSjet >=workingpoint)
temp_2_failed_2 = ~(events.SubTau.idDeepTau2017v2p1VSjet >=workingpoint)
temp_2_failed = temp_2_failed_1 & temp_2_failed_2
events_2_failed = events[temp_2_failed]
temp_1_failed_1_passed_1 = events.LeadTau.idDeepTau2017v2p1VSjet >=workingpoint
temp_1_failed_1_passed_2 = ~events.SubTau.idDeepTau2017v2p1VSjet >=workingpoint
temp_1_failed_1_passed_3 = ~events.LeadTau.idDeepTau2017v2p1VSjet >=workingpoint
temp_1_failed_1_passed_4 = events.SubTau.idDeepTau2017v2p1VSjet >=workingpoint
temp_1_failed_1_passed = (temp_1_failed_1_passed_1 & temp_1_failed_1_passed_2) | (temp_1_failed_1_passed_3 & temp_1_failed_1_passed_4) 
events_1_failed_1_passed = events[temp_1_failed_1_passed]
print("***************************")
print(len(events_1_failed_1_passed))
print(len(events_2_failed))
print(len(events_2_passed))
print("***************************")

print("2 passed, 2 failed, 1 failed-1 passed are created!")
#####
#binning for pt tau observable (analysis based on it)
#####
bins_tau_pt = hist.Bin("pT", "tau pT [GeV]",50, 20, 220)
#####
#compute the pt binning of 2 passed etc...
######

h_2P = TH1D("2 taus passed","2 taus passed", 50, 20, 220)

temp = events_2_passed.LeadTau.pt
temp_weight = events_2_passed.genWeight

for i in range(0, len(temp)):
    h_2P.Fill(temp[i], temp_weight[i])

h_1P1F = TH1D("1 tau passed 1 failed","1 tau passed 1 failed", 50, 20, 220)
temp = events_1_failed_1_passed.LeadTau.pt
temp_weight = events_1_failed_1_passed.genWeight
for i in range(0, len(temp)):
    h_1P1F.Fill(temp[i], temp_weight[i])
h_2F = TH1D("2 taus failed","2 taus failed", 50, 20, 220)
temp = events_2_failed.LeadTau.pt
temp_weight = events_2_failed.genWeight
for i in range(0, len(temp)):
    h_2F.Fill(temp[i], temp_weight[i])   
print("histogram of 2P, 1PF etc... are made")

h_2F.Write()
h_2P.Write()
h_1P1F.Write()
#compute the prompt rate
# which is the rate at which prompt taus pass the selection
#first choose the prompt taus and then apply the selection, compute the ratio
print("#######################################################")
print("Prompt rate computation")

events_pr = NanoEventsFactory.from_root(file, schemaclass=NanoAODSchema).events()
print("Events are charged!")
events_pr['SelElectron'] = events_pr.Electron[(events_pr.Electron.pt > 24.) & (events_pr.Electron.mvaFall17V2Iso_WP90 > 0.5)]
events_pr['SelMuon'] = events_pr.Muon[(events_pr.Muon.pt > 24.) & (events_pr.Muon.mediumPromptId) & (events_pr.Muon.pfRelIso03_all < 0.2) & (np.abs(events_pr.Muon.dxy) < 0.005)]
muon1, electron1 = ak.unzip(ak.cartesian([events_pr.SelMuon, events_pr.SelElectron], nested=True))
match1 = ak.any(muon1.jetIdx == electron1.jetIdx, axis=-1, mask_identity=False)        
events_pr['SelMuon'] = events_pr.SelMuon[(~(match1))]  
events_pr['SelTau'] = events_pr.Tau[(events_pr.Tau.pt > 20.) & (abs(events_pr.Tau.eta) < 2.3)& (events_pr.Tau.idDeepTau2017v2p1VSmu > 0.5) & (events_pr.Tau.idDeepTau2017v2p1VSe > 0.5) & (events_pr.Tau.idDeepTau2017v2p1VSjet >=7)]
tau2, electron2 = ak.unzip(ak.cartesian([events_pr.SelTau, events_pr.SelElectron], nested=True))
match2 = ak.any(tau2.jetIdx == electron2.jetIdx, axis=-1, mask_identity=False)
tau3, muon3 = ak.unzip(ak.cartesian([events_pr.SelTau, events_pr.SelMuon], nested=True))
match3 = ak.any(tau3.jetIdx == muon3.jetIdx, axis=-1, mask_identity=False)
events_pr['SelTau'] = events_pr.SelTau[((~(match2) & ~(match3)) )]  
print("Selection of leptons is made!")
prompt_before_sel = ak.flatten(events_pr.SelTau.pt)

h_prompt_bef = TH1D("prompt taus before selection","prompt taus before selection", 50, 20, 220)
for i in range(0, len(prompt_before_sel)):
    h_prompt_bef.Fill(prompt_before_sel[i])
print("histogram of denum of prompt rate is made")

#after selection
events_pr['SelTau'] = events_pr.SelTau[events_pr.SelTau.idDeepTau2017v2p1VSjet >=15]
prompt_after_sel = ak.flatten(events_pr.SelTau.pt)
h_prompt_aft = TH1D("prompt taus after selection","prompt taus after selection", 50, 20, 220)
for i in range(0, len(prompt_after_sel)):
    h_prompt_aft.Fill(prompt_after_sel[i])
print("histogram of num of prompt rate is made")

h_prompt_rate = h_prompt_aft/h_prompt_bef
h_prompt_rate.SetName("histogram of prompt rate")
h_prompt_rate.SetTitle("histogram of prompt rate")

print("histogram of prompt rate is made")
h_prompt_rate.Write()

###compute the fake rate
### with pt tau
### in the CR region (one bjet) and before and after ID tau selection
print("#######################################################")

print("fake rate computation")
nb_of_tau = 1
nb_of_e = 1
nb_of_mu = 1
events_fr = NanoEventsFactory.from_root(file, schemaclass=NanoAODSchema).events()
print("events are charged!")
events_fr['SelElectron'] = events_fr.Electron[(events_fr.Electron.pt > 24.) & (events_fr.Electron.mvaFall17V2Iso_WP90 > 0.5)]
events_fr['SelMuon'] = events_fr.Muon[(events_fr.Muon.pt > 24.) & (events_fr.Muon.mediumPromptId) & (events_fr.Muon.pfRelIso03_all < 0.2) & (np.abs(events_fr.Muon.dxy) < 0.005)]
muon1, electron1 = ak.unzip(ak.cartesian([events_fr.SelMuon, events_fr.SelElectron], nested=True))
match1 = ak.any(muon1.jetIdx == electron1.jetIdx, axis=-1, mask_identity=False)        
events_fr['SelMuon'] = events_fr.SelMuon[(~(match1))]    
events_fr['SelTau'] = events_fr.Tau[(events_fr.Tau.pt > 20.) & (abs(events_fr.Tau.eta) < 2.3)& (events_fr.Tau.idDeepTau2017v2p1VSmu > 0.5) & (events_fr.Tau.idDeepTau2017v2p1VSe > 0.5) & (events_fr.Tau.idDeepTau2017v2p1VSjet >=7)]
tau2, electron2 = ak.unzip(ak.cartesian([events_fr.SelTau, events_fr.SelElectron], nested=True))
match2 = ak.any(tau2.jetIdx == electron2.jetIdx, axis=-1, mask_identity=False)
tau3, muon3 = ak.unzip(ak.cartesian([events_fr.SelTau, events_fr.SelMuon], nested=True))
match3 = ak.any(tau3.jetIdx == muon3.jetIdx, axis=-1, mask_identity=False)
events_fr['SelTau'] = events_fr.SelTau[((~(match2) & ~(match3)) )]  
print("Selection of leptons is made!")
events_fr = events_fr[ak.num(events_fr.SelMuon) == nb_of_mu]
events_fr = events_fr[ak.num(events_fr.SelElectron) == nb_of_e]
jets = events_fr.Jet[events_fr.Jet.pt > 25.]
bjets = jets[jets.btagDeepFlavB > 0.2770]
one_bjet = ak.num(bjets) == 1
no_bjet = ak.num(bjets) == 0
print(str(nb_of_e)+" electron "+str(nb_of_mu)+" muon, and "+str(nb_of_tau)+" taus required!")

CR_loose = events_fr[one_bjet]
CR_medium = events_fr[one_bjet]
print("events with one bjet are selected")
CR_medium['SelTau'] = CR_medium.SelTau[CR_medium.SelTau.idDeepTau2017v2p1VSjet >=workingpoint]
print("events with one bjet and medium passed are selected")

CR_loose = CR_loose[ak.num(CR_loose.SelTau) == nb_of_tau]
CR_medium = CR_medium[ak.num(CR_medium.SelTau) == nb_of_tau]
CR_loose_kin = ak.flatten(CR_loose.SelTau.pt)
CR_medium_kin = ak.flatten(CR_medium.SelTau.pt)


h_fakerate_bef = TH1D("fake rate taus before selection","fake rate taus before selection", 50, 20, 220)
for i in range(0, len(CR_loose_kin)):
    h_fakerate_bef.Fill(CR_loose_kin[i])
print("histogram of fake rate loose is made")

h_fakerate_aft = TH1D("fake rate taus after selection","fake rate taus after selection", 50, 20, 220)
for i in range(0, len(CR_medium_kin)):
    h_fakerate_aft.Fill(CR_medium_kin[i])
print("histogram of fake rate loose  is made")




h_fake_rate = h_fakerate_aft/h_fakerate_bef
h_fake_rate.SetName("histogram of ratio of fake rate")
h_fake_rate.SetTitle("histogram of ratio of fake rate")

print("histogram of ratio of fake rate is made")
h_fake_rate.Write()


######
#COMPUTATION OF THE OVERALL DILEPTONS FAKERATE
######
print("computation of the dilepton fake rate via the formula: ")
print("N signal = p^2 * Npp")
print("Npp = (1/(p-f)^2) * ( (1-f)^2 * 2P - f*(1-f)*1P1F + f^2* 2F)")
h_1 = TH1D("histogram filled with one","histogram filled with one", 50, 20, 220)
for i in range(20, 220, 4):
    h_1.Fill(i)
h_1.Write()

temp1 = h_1/((h_prompt_rate-h_fake_rate)*(h_prompt_rate-h_fake_rate))
temp2 = (h_1-h_fake_rate)*(h_1-h_fake_rate)
temp3 = h_fake_rate*(h_1-h_fake_rate)
temp4 = h_fake_rate*h_fake_rate
N_pp = temp1 * ( temp2 * h_2P - temp3 * h_1P1F + temp4 * h_2F)

temp5 = h_prompt_rate*h_prompt_rate
N_sig = temp5*N_pp
N_sig.SetName("Number of signal")
N_sig.SetTitle("Number of signal")

N_sig.Write()
f.Close()
print("all the histograms are stored in: " + str(path_root_file_histo))



