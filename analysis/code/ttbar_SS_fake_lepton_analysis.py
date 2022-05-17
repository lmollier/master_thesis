import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from helpers import delta_r, delta_phi, inv_mass_3p, cos_opening_angle, inv_mass, delta_eta
from samples import signal_samples
from helpers import files_from_dir, files_from_dirs
import numpy as np
import matplotlib.pyplot as plt

bck_dir = '/eos/user/l/lmollier/PDM/data/bck_v3/MC_background/'
local_dir_bck_TT_To_2L2Nu = bck_dir+'TT_To_2L2Nu'



def jet_veto_selection(events, dr_cut=0.5):
    jets = events.Jet[events.Jet.pt > 30.]
    l1 = events.SelLepton[:,0]
    l2 = events.SelLepton[:,1]
    l3 = events.SelLepton[:,2]
    dr1 = jets.delta_r(l1)
    jets = jets[dr1 > 0.5]
    dr2 = jets.delta_r(l2)
    jets = jets[dr2 > 0.5]
    dr3 = jets.delta_r(l3)
    jets = jets[dr3 > 0.5]
    events  = events[ak.num(jets) <= 2]
    return events

def bjet_veto(events):
    jets = events.Jet[events.Jet.pt > 25.]
    bjets = jets[jets.btagDeepFlavB > 0.2770]
    l1 = events.SelLepton[:,0]
    l2 = events.SelLepton[:,1]
    l3 = events.SelLepton[:,2]
    dr1 = bjets.delta_r(l1)
    bjets = bjets[dr1 > 0.5]
    dr2 = bjets.delta_r(l2)
    bjets = bjets[dr2 > 0.5]
    dr3 = bjets.delta_r(l3)
    bjets = bjets[dr3 > 0.5]
    cut = (ak.num(bjets) == 0)
    events = events[cut]
    return events




samples = files_from_dir(local_dir_bck_TT_To_2L2Nu),
i_file = 1
tot_fake_e = 0
tot_fake_mu = 0
tot_fake_tau = 0
tot_not_viable = 0
tot_events_selected = 0
for file in samples[0]:
    NanoAODSchema.warn_missing_crossrefs = False
    events = NanoEventsFactory.from_root(file, schemaclass=NanoAODSchema).events()
    print(str(i_file)+"/"+str(len(samples[0])))
    i_file+=1
    events['weight'] = events.genWeight
    events['SelElectron'] = events.Electron[(events.Electron.pt > 24.) & (events.Electron.mvaFall17V2Iso_WP90 > 0.5)]
    events['SelMuon'] = events.Muon[(events.Muon.pt > 24.) & (events.Muon.mediumPromptId) & (events.Muon.pfRelIso03_all < 0.2) & (np.abs(events.Muon.dxy) < 0.005)]
    muon1, electron1 = ak.unzip(ak.cartesian([events.SelMuon, events.SelElectron], nested=True))
    match1 = ak.any(muon1.jetIdx == electron1.jetIdx, axis=-1, mask_identity=False)
    notinjet1 = ak.any(muon1.jetIdx == -1) 
    events['SelMuon'] = events.SelMuon[(~(match1))]    
    events['OtherMuon'] = events.Muon[~( (events.Muon.pt > 24.) & (events.Muon.mediumPromptId) & (events.Muon.pfRelIso03_all < 0.2) & (np.abs(events.Muon.dxy) < 0.005))]
    events['SelTau'] = events.Tau[(events.Tau.pt > 20.) & (abs(events.Tau.eta) < 2.3)& (events.Tau.idDeepTau2017v2p1VSmu > 0.5) & (events.Tau.idDeepTau2017v2p1VSe > 0.5) & (events.Tau.idDeepTau2017v2p1VSjet >=8)]
    tau2, electron2 = ak.unzip(ak.cartesian([events.SelTau, events.SelElectron], nested=True))
    match2 = ak.any(tau2.jetIdx == electron2.jetIdx, axis=-1, mask_identity=False)
    tau3, muon3 = ak.unzip(ak.cartesian([events.SelTau, events.SelMuon], nested=True))
    match3 = ak.any(tau3.jetIdx == muon3.jetIdx, axis=-1, mask_identity=False)
    notinjet2 = ak.any(tau2.jetIdx == int(-1))
    notinjet3 = ak.any(tau3.jetIdx == int(-1))
    events['SelTau'] = events.SelTau[((~(match2) & ~(match3)) )]  
    events = events[ak.num(events.SelElectron) + ak.num(events.SelMuon) + ak.num(events.SelTau) == 3]   
    events['SelLepton'] = ak.concatenate([events.SelElectron, events.SelMuon, events.SelTau],axis = -1)
    events = jet_veto_selection(events)
    events = bjet_veto(events)
    events = events[ak.num(events.SelElectron) ==  1] 
    events = events[ak.num(events.SelTau) ==  1]   
    events = events[ak.num(events.SelMuon) ==  1]
    same_sign = (events.SelElectron[:,0].charge * events.SelMuon[:,0].charge) > 0
    events_SS = events[same_sign]
    #events_OS = events[~same_sign]
    #print(len(events_SS))
    tot_events_selected+= len(events_SS)
    # print("lepton (e mu tau) flav:")
    # print(events_SS.SelLepton.genPartFlav)
    # print("fake electron only")
    fake_e_sel = ak.any(events_SS.SelElectron.genPartFlav != 1, axis = -1) &  ak.any(events_SS.SelMuon.genPartFlav == 1, axis = -1) &   (ak.any(events_SS.SelTau.genPartFlav == 3, axis = -1) |  ak.any(events_SS.SelTau.genPartFlav == 4, axis = -1) |  ak.any(events_SS.SelTau.genPartFlav == 5, axis = -1) )
    fake_electrons = events_SS[fake_e_sel]

    tot_fake_e+=len(fake_electrons)
    # print("fake muon only")
    fake_mu_sel = ak.any(events_SS.SelElectron.genPartFlav == 1, axis = -1) &  ak.any(events_SS.SelMuon.genPartFlav != 1, axis = -1) &   (ak.any(events_SS.SelTau.genPartFlav == 3, axis = -1) |  ak.any(events_SS.SelTau.genPartFlav == 4, axis = -1) |  ak.any(events_SS.SelTau.genPartFlav == 5, axis = -1) )
    fake_muons = events_SS[fake_mu_sel]
    tot_fake_mu+=len(fake_muons)
    # print(fake_muons.SelLepton.genPartFlav)
    # print(len(fake_muons))

    # print("fake tau only")
    fake_tau_sel = ak.any(events_SS.SelElectron.genPartFlav == 1, axis = -1) &  ak.any(events_SS.SelMuon.genPartFlav == 1, axis = -1) &   (ak.any(events_SS.SelTau.genPartFlav != 3, axis = -1) &  ak.any(events_SS.SelTau.genPartFlav != 4, axis = -1) &  ak.any(events_SS.SelTau.genPartFlav != 5, axis = -1) )
    fake_taus = events_SS[fake_tau_sel]
    tot_fake_tau+=len(fake_taus)
    # print(fake_taus.SelLepton.genPartFlav)
    # print(len(fake_taus))


    # print("not one of the three above")
    unknown = events_SS[~fake_tau_sel &~fake_mu_sel & ~fake_e_sel] 
    tot_not_viable+= len(unknown)
    # print(unknown.SelLepton.genPartFlav)
    # print(len(unknown))
    # print(len(unknown)+len(fake_muons)+len(fake_taus)+len(fake_electrons))
   

# print(tot_events_selected)
# print(tot_fake_e)
# print(tot_fake_mu)
# print(tot_fake_tau)
# print(tot_not_viable)

print("\n\n\n\n\n")
print(u'\u2500' * 65)
print(f'| {"total (SS)":<10}', end='')
print(f'| {"fake e":<10}', end='')
print(f'| {"fake mu ":<11}', end='')
print(f'| {"fake tau  ":<12}', end='')
print(f'| {"others     |":<10}', end='\n')
print(u'\u2500' * 65)
print(f'| {tot_events_selected:<10}', end='')
print(f'| {tot_fake_e:<10}', end='')
print(f'| {tot_fake_mu:<11}', end='')
print(f'| {tot_fake_tau:<12}', end='')
print(f'| {tot_not_viable:<10} |', end='\n')
print(u'\u2500' * 65)
print("\n\n\n\n\n")


