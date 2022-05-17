
from faulthandler import is_enabled
from ssl import HAS_TLSv1_3
from statistics import multimode
import numpy as np
import awkward as ak
from coffea import processor, hist
from pyparsing import javaStyleComment
from helpers import delta_r, delta_phi, inv_mass_3p, cos_opening_angle, inv_mass, delta_eta

class HNLAnalysis(processor.ProcessorABC):
    def __init__(self):
        ds_axis = hist.Cat("ds", "Primary dataset")
        acc_dict = {var: hist.Hist("Counts", ds_axis, axis) for var, axis in self.get_var_axis_pairs()}
        print("using get_selection_test!!!")
        self.selections = self.get_selections_test()
        for stage in self.selections:
            acc_dict[f'n_ev_{stage}'] = processor.defaultdict_accumulator(int)
            acc_dict[f'sumw_{stage}'] = processor.defaultdict_accumulator(float)
        acc_dict['sel_array'] = processor.column_accumulator(np.ndarray((0, 12)))
        self._accumulator = processor.dict_accumulator(acc_dict)
    
    @property
    def accumulator(self):
        return self._accumulator
    
    @staticmethod
    def get_var_axis_pairs():
        
        mass_axis = hist.Bin("mass", r"$m_{\ell\ell}$ [GeV]", 300, 0., 1500.)
        pt_axis = hist.Bin("pt", r"$p_{T}$ [GeV]", 300, 0., 1500)
        mt_axis = hist.Bin("mt", r"Transverse mass [GeV]", 300, 0., 1500.)
        mc_axis = hist.Bin("mc", r"Combined mass [GeV]", 300, 0., 1500.)

        eta_axis = hist.Bin('eta', r'$\eta$', 30,  -3.1415927, 3.1415927)
        phi_axis = hist.Bin('phi', r'$\phi$', 30, -3.1415927, 3.1415927)
        dr_axis = hist.Bin("dr", r"$\Delta R$", 30, 0., 5.)
        charge_axis = hist.Bin('charge', r'Charge', 5, -2.5, 2.5)
 
        v_a_pairs = [
            ('pt_tau_1', pt_axis),
            ('eta_tau_1', eta_axis),
            ('phi_tau_1', phi_axis),
            ('charge_tau_1', charge_axis),
            ('mass_tau_1', mass_axis),

            ('pt_mu_1', pt_axis),
            ('eta_mu_1', eta_axis),
            ('phi_mu_1', phi_axis),
            ('charge_mu_1', charge_axis),
            ('mass_mu_1', mass_axis),

            ('pt_e_1', pt_axis),
            ('eta_e_1', eta_axis),
            ('phi_e_1', phi_axis),
            ('charge_e_1', charge_axis),
            ('mass_e_1', mass_axis),

            ('pt_jet_1', pt_axis),
            ('pt_jet_2', pt_axis),

        ('pt_sum_mutau', pt_axis),
        ('pt_sum_etau', pt_axis),
        ('pt_sum_emu', pt_axis),
        ('pt_sum_emutau', pt_axis),
        ('pt_sum_mutauMET', pt_axis),
        ('pt_sum_etauMET', pt_axis),
        ('pt_sum_emuMET', pt_axis),
        ('pt_sum_emutauMET', pt_axis),
        ('comb_mass_mutau', mc_axis),
        ('comb_mass_etau', mc_axis),
        ('comb_mass_emu', mc_axis),


        ('dr_mutau', dr_axis),
        ('dr_etau', dr_axis),
        ('dr_emu', dr_axis),

        ('dphi_mutau', phi_axis),
        ('dphi_etau', phi_axis),
        ('dphi_emu', phi_axis),

        ('deta_mutau', eta_axis),
        ('deta_etau', eta_axis),
        ('deta_emu', eta_axis),


        ('transverse_mass_mutau', mt_axis),
        ('transverse_mass_etau', mt_axis),
        ('transverse_mass_emu', mt_axis),
        ('transverse_mass_muMET', mt_axis),
        ('transverse_mass_eMET', mt_axis),
        ('transverse_mass_tauMET', mt_axis),
        ('MT_total', mt_axis),


        ('pt_MET', pt_axis),
        ('phi_MET', phi_axis),




        ]
        tags = ['_emutau']
        for tag in tags:
            v_a_pairs = [(name+tag, axis) for name, axis in v_a_pairs]
        tags_sign = ['_SS']
        for tag in tags_sign:
            v_a_pairs_1 = [(name+tag, axis) for name, axis in v_a_pairs]
        tags_sign = ['_OS']
        for tag in tags_sign:
            v_a_pairs_2 = [(name+tag, axis) for name, axis in v_a_pairs]

        return v_a_pairs_1 + v_a_pairs_2

    @staticmethod
    def get_selections():
        return [
            'all', 
            'incoming', 
            'lead_muon_id', 
            'isomu24',
            'bjet_veto'
        ]
    @staticmethod
    def get_selections_test():
        return [
            'all', 
            '3leptons', 
            'jet_veto', 
            'bjet_veto',
            'emutau'
        ]
    @staticmethod
    def get_selections_signal_emutau():
        return [
            'all', 
            'selection_electron',
            'selection_muon',
            'selection_tau',
            'selection_3leptons'
            'bjetveto',
            'jetveto'
        ]    
    
    # we will receive a NanoEvents instead of a coffea DataFrame
    def process(self, events):
        out = self.accumulator.identity()
        ds = events.metadata["dataset"]
        # To add a new column to events, use [] syntax on original (!) events - all others are views only
        # If we process a given dataset twice, need to work on views
        

        if 'Data' in ds:
            events['genWeight'] = events.run > 0
        
        events['weight'] = events.genWeight
        print("For dataset : "+ ds)
        out['sumw_all'][ds] += ak.sum(events.genWeight)
        events['SelElectron'] = events.Electron[(events.Electron.pt > 24.) & (events.Electron.mvaFall17V2Iso_WP80 > 0.5)  ]

        events['SelMuon'] = events.Muon[(events.Muon.pt > 24.) & (events.Muon.mediumPromptId) & (events.Muon.pfRelIso03_all < 0.2) & (np.abs(events.Muon.dxy) < 0.005)]
        muon1, electron1 = ak.unzip(ak.cartesian([events.SelMuon, events.SelElectron], nested=True))
        match1 = ak.any(muon1.jetIdx == electron1.jetIdx, axis=-1, mask_identity=False)
        notinjet1 = ak.any(muon1.jetIdx == -1) 
        events['SelMuon'] = events.SelMuon[(~(match1))]    
        events['SelTau'] = events.Tau[(events.Tau.pt > 20.) & (abs(events.Tau.eta) < 2.3)& (events.Tau.idDeepTau2017v2p1VSmu >=4) & (events.Tau.idDeepTau2017v2p1VSe >= 4) & (events.Tau.idDeepTau2017v2p1VSjet >=16)]


        tau2, electron2 = ak.unzip(ak.cartesian([events.SelTau, events.SelElectron], nested=True))
        match2 = ak.any(tau2.jetIdx == electron2.jetIdx, axis=-1, mask_identity=False)
        tau3, muon3 = ak.unzip(ak.cartesian([events.SelTau, events.SelMuon], nested=True))
        match3 = ak.any(tau3.jetIdx == muon3.jetIdx, axis=-1, mask_identity=False)
        notinjet2 = ak.any(tau2.jetIdx == int(-1))
        notinjet3 = ak.any(tau3.jetIdx == int(-1))
        
        events['SelTau'] = events.SelTau[((~(match2) & ~(match3)) )]  
        events = events[ak.num(events.SelElectron) + ak.num(events.SelMuon) + ak.num(events.SelTau) == 3]   
        out['sumw_3leptons'][ds] +=ak.sum(events.genWeight)


        #events = events[events.MET.pt < 75]
        events['SelLepton'] = ak.concatenate([events.SelElectron, events.SelMuon, events.SelTau],axis = -1)
        events = events[events.HLT.IsoMu24]

        
        
        
        events = self.jet_veto_selection(events)
        out['sumw_jet_veto'][ds] += ak.sum(events.genWeight)
        events = self.bjet_veto(events)
        out['sumw_bjet_veto'][ds] += ak.sum(events.genWeight)

        # Allow for a dataset to be processed multiple times, e.g. for background estimation
        self.process_ds(ds, out, events)
        # e.g. self.process_ds(ds, out, events, mode='InvertedMuonIso')

        return out
    
    def process_ds(self, ds, out, events, mode=''):

        # out['n_ev_all'][ds] += len(events)
        # out['n_ev_incoming'][ds] += len(events)

        # out['sumw_all'][ds] += ak.sum(events.genWeight)
        # #out['sumw_all'][ds] += len(events)

        # out['sumw_incoming'][ds] += ak.sum(events.genWeight)
        
        ### Example how to apply a signal-only selection

        # if 'HNL' in ds:
        #     # Require two muons from HNL decay
        #     hnl_muons = events.GenPart[(np.abs(events.GenPart.pdgId) == 13) & (events.GenPart.status==1) & (events.GenPart.distinctParentIdxG != -1)]
        #     hnl_muons = hnl_muons[hnl_muons.distinctParent.pdgId == 9900012]
        #     events = events[ak.num(hnl_muons) == 2]
        #     hnl_muons = hnl_muons[ak.num(hnl_muons) == 2]
        



        # out['n_ev_lead_muon_id'][ds] += len(events)
        # out['sumw_lead_muon_id'][ds] += ak.sum(events.genWeight)

        #events = events[events.HLT.IsoMu24]
        # out['n_ev_isomu24'][ds] += len(events)
        # out['sumw_isomu24'][ds] += ak.sum(events.genWeight)
        

        
        


        self.analyse_emutau(events, out, ds, mode)
        #self.analyse_etautau(events, out, ds, mode)
        #self.analyse_mutautau(events, out, ds, mode)

        # # # print information about run/lumi/event
        # if ds in ['W', 'WSS']:
        #     for ev in events:
        #         print(ds)
        #         print('run:', ev.run, ', lumi:', ev.luminosityBlock, ', ev:', ev.event)
    

    def save_tag(self, events, out, ds, mode, tag=''):
        self.save_global(events, out, ds, mode, tag)
        self.save_muon(events, out, ds, mode, tag)
        self.save_electron(events, out, ds, mode, tag)
        self.save_tau(events, out, ds, mode, tag)
        self.save_sum_pt(events, out, ds, mode, tag)
        self.save_combined_mass(events, out, ds, mode, tag)
        self.save_combined_var(events, out, ds, mode, tag)
        self.save_MET(events, out, ds, mode, tag)

        self.save_transverse_mass(events, out, ds, mode, tag)
        # self.save_inv_mass(events, out, ds, mode, tag)
        # self.save_mass_HNL(events, out, ds, mode, tag)
        # self.save_3_highest_pt_mass(events, out, ds, mode, tag)
        # self.save_ratio_final_state(events, out, ds, mode, tag)

        



    def save_global(self, events, out, ds, mode, tag=''):
        print("save_global " + ds)
        jet_1_events = events[ak.num(events.Jet)>=1]
        jet_2_events = events[ak.num(events.Jet)>=2]
        out[f'pt_jet_1{tag}'].fill(ds=ds, pt=ak.flatten(jet_1_events.Jet[:, 0].pt, axis=None), weight=jet_1_events.weight)
        out[f'pt_jet_2{tag}'].fill(ds=ds, pt=ak.flatten(jet_2_events.Jet[:, 1].pt, axis=None), weight=jet_2_events.weight)

    def save_mass_HNL(self, events, out, ds, mode, tag=''):
        print("TODO inv mass of HNL")

    

    def save_sum_pt(self, events, out, ds, mode, tag=''):
        print("save_sum_pt " + ds)
        lead_muon = events.SelMuon[:, 0]
        lead_electron= events.SelElectron[:, 0]
        lead_tau = events.SelTau[:, 0]
        MET_pt = events.MET.pt
        out[f'pt_sum_mutau{tag}'].fill(ds=ds, pt=ak.flatten((lead_muon.pt + lead_tau.pt), axis=None), weight=events.weight)
        out[f'pt_sum_etau{tag}'].fill(ds=ds, pt=ak.flatten((lead_electron.pt + lead_tau.pt), axis=None), weight=events.weight)
        out[f'pt_sum_emu{tag}'].fill(ds=ds, pt=ak.flatten((lead_muon.pt + lead_electron.pt), axis=None), weight=events.weight)
        out[f'pt_sum_emutau{tag}'].fill(ds=ds, pt=ak.flatten((lead_muon.pt + lead_electron.pt  + lead_tau.pt ), axis=None), weight=events.weight)
        out[f'pt_sum_mutauMET{tag}'].fill(ds=ds, pt=ak.flatten((lead_muon.pt + lead_tau.pt + MET_pt), axis=None), weight=events.weight)
        out[f'pt_sum_etauMET{tag}'].fill(ds=ds, pt=ak.flatten((lead_electron.pt + lead_tau.pt+ MET_pt), axis=None), weight=events.weight)
        out[f'pt_sum_emuMET{tag}'].fill(ds=ds, pt=ak.flatten((lead_muon.pt + lead_electron.pt+ MET_pt), axis=None), weight=events.weight)
        out[f'pt_sum_emutauMET{tag}'].fill(ds=ds, pt=ak.flatten((lead_muon.pt + lead_electron.pt  + lead_tau.pt+ MET_pt ), axis=None), weight=events.weight)
    
    
    def save_combined_mass(self, events, out, ds, mode, tag=''):
        print("save_combined_mass  " + ds)
        lead_muon = events.SelMuon[:, 0]
        lead_electron= events.SelElectron[:, 0]
        lead_tau = events.SelTau[:, 0]
        out[f'comb_mass_mutau{tag}'].fill(ds=ds, mc=ak.flatten((lead_muon + lead_tau).mass, axis=None), weight=events.weight)
        out[f'comb_mass_etau{tag}'].fill(ds=ds, mc=ak.flatten((lead_electron+ lead_tau).mass, axis=None), weight=events.weight)
        out[f'comb_mass_emu{tag}'].fill(ds=ds, mc=ak.flatten((lead_muon + lead_electron).mass, axis=None), weight=events.weight)

    def save_combined_var(self, events, out, ds, mode, tag=''):
        print("save_combined_var  " + ds)
        lead_muon = events.SelMuon[:, 0]
        lead_electron= events.SelElectron[:, 0]
        lead_tau = events.SelTau[:, 0]
        out[f'dr_mutau{tag}'].fill(ds=ds, dr=ak.flatten(delta_r(lead_tau,lead_muon), axis=None), weight=events.weight)
        out[f'dr_etau{tag}'].fill(ds=ds, dr=ak.flatten(delta_r(lead_electron,lead_tau), axis=None), weight=events.weight)
        out[f'dr_emu{tag}'].fill(ds=ds, dr=ak.flatten(delta_r(lead_electron,lead_muon), axis=None), weight=events.weight)

        out[f'dphi_mutau{tag}'].fill(ds=ds, phi=ak.flatten(delta_phi(lead_tau,lead_muon), axis=None), weight=events.weight)
        out[f'dphi_etau{tag}'].fill(ds=ds, phi=ak.flatten(delta_phi(lead_electron,lead_tau), axis=None), weight=events.weight)
        out[f'dphi_emu{tag}'].fill(ds=ds, phi=ak.flatten(delta_phi(lead_electron,lead_muon), axis=None), weight=events.weight)

        out[f'deta_mutau{tag}'].fill(ds=ds, eta=ak.flatten(delta_eta(lead_tau,lead_muon), axis=None), weight=events.weight)
        out[f'deta_etau{tag}'].fill(ds=ds, eta=ak.flatten(delta_eta(lead_electron,lead_tau), axis=None), weight=events.weight)
        out[f'deta_emu{tag}'].fill(ds=ds, eta=ak.flatten(delta_eta(lead_electron,lead_muon), axis=None), weight=events.weight)


    def save_transverse_mass(self, events, out, ds, mode, tag=''):
        print("save_transverse_mass  " + ds)
        lead_muon = events.SelMuon[:, 0]
        lead_electron= events.SelElectron[:, 0]
        lead_tau = events.SelTau[:, 0]
        muon_mass_sq = lead_muon.mass**2
        tau_mass_sq = lead_tau.mass**2
        electron_mass_sq = lead_electron.mass**2
        muon_pt = lead_muon.pt
        electron_pt = lead_electron.pt
        tau_pt = lead_tau.pt
        muon_phi = lead_muon.phi
        electron_phi = lead_electron.phi
        tau_phi = lead_tau.phi
        MET_phi = events.MET.phi
        MET_pt = events.MET.pt
        muon_energy = np.sqrt(muon_mass_sq + muon_pt**2)
        electron_energy = np.sqrt(electron_mass_sq + electron_pt**2)
        tau_energy = np.sqrt(tau_mass_sq + tau_pt**2)
        cos_emu = np.cos(abs(electron_phi - muon_phi))
        cos_etau = np.cos(abs(electron_phi - tau_phi))
        cos_mutau = np.cos(abs(tau_phi - muon_phi))
        cos_eMET = np.cos(abs(electron_phi - MET_phi))
        cos_muMET = np.cos(abs(MET_phi - muon_phi))
        cos_tauMET = np.cos(abs(tau_phi - MET_phi))


        transverse_mass_emu = np.sqrt(electron_mass_sq + muon_mass_sq + 2*(electron_energy*muon_energy - electron_pt*muon_pt*cos_emu))
        transverse_mass_etau = np.sqrt(electron_mass_sq + tau_mass_sq + 2*(electron_energy*tau_energy - electron_pt*tau_pt*cos_etau))
        transverse_mass_mutau = np.sqrt(tau_mass_sq + muon_mass_sq + 2*(muon_energy*tau_energy - tau_pt*muon_pt*cos_mutau))
        transverse_mass_eMET = np.sqrt(electron_mass_sq  + 2*(electron_energy*MET_pt - electron_pt*MET_pt*cos_eMET))
        transverse_mass_muMET =np.sqrt(muon_mass_sq  + 2*(muon_energy*MET_pt - muon_pt*MET_pt*cos_muMET))
        transverse_mass_tauMET =np.sqrt(tau_mass_sq  + 2*(tau_energy*MET_pt - tau_pt*MET_pt*cos_tauMET))


        MT_total = np.sqrt(transverse_mass_emu**2 +transverse_mass_etau**2 +transverse_mass_mutau**2 +transverse_mass_eMET**2 +transverse_mass_muMET**2 +transverse_mass_tauMET**2  )        
        out[f'transverse_mass_mutau{tag}'].fill(ds=ds, mt=ak.flatten(transverse_mass_mutau, axis=None), weight=events.weight)
        out[f'transverse_mass_etau{tag}'].fill(ds=ds, mt=ak.flatten(transverse_mass_etau, axis=None), weight=events.weight)
        out[f'transverse_mass_emu{tag}'].fill(ds=ds, mt=ak.flatten(transverse_mass_emu, axis=None), weight=events.weight)
        out[f'transverse_mass_muMET{tag}'].fill(ds=ds, mt=ak.flatten(transverse_mass_muMET, axis=None), weight=events.weight)
        out[f'transverse_mass_eMET{tag}'].fill(ds=ds, mt=ak.flatten(transverse_mass_eMET, axis=None), weight=events.weight)
        out[f'transverse_mass_tauMET{tag}'].fill(ds=ds, mt=ak.flatten(transverse_mass_tauMET, axis=None), weight=events.weight)
        out[f'MT_total{tag}'].fill(ds=ds, mt=MT_total, weight=events.weight)

    def save_muon(self, events, out, ds, mode, tag=''):
        print("save_muon " + ds)
        lead_muon = events.SelMuon[:, 0]
        out[f'pt_mu_1{tag}'].fill(ds=ds, pt=ak.flatten(lead_muon.pt, axis=None), weight=events.weight)
        out[f'eta_mu_1{tag}'].fill(ds=ds, eta=ak.flatten(lead_muon.eta, axis=None), weight=events.weight)
        out[f'phi_mu_1{tag}'].fill(ds=ds, phi=ak.flatten(lead_muon.phi, axis=None), weight=events.weight)
        out[f'charge_mu_1{tag}'].fill(ds=ds, charge=ak.flatten(lead_muon.charge, axis=None), weight=events.weight)
        out[f'mass_mu_1{tag}'].fill(ds=ds, mass=ak.flatten(lead_muon.mass, axis=None), weight=events.weight)


        # events_temp = events[ak.num(events.SelMuon)>=2]
        # sublead_muon = events_temp.SelMuon[:,1]
        # out[f'pt_mu_2{tag}'].fill(ds=ds, pt=ak.flatten(sublead_muon.pt, axis=None), weight=events_temp.weight)
        # out[f'eta_mu_2{tag}'].fill(ds=ds, eta=ak.flatten(sublead_muon.eta, axis=None), weight=events_temp.weight)
        # out[f'phi_mu_2{tag}'].fill(ds=ds, phi=ak.flatten(sublead_muon.phi, axis=None), weight=events_temp.weight)
        # out[f'charge_mu_2{tag}'].fill(ds=ds, charge=ak.flatten(sublead_muon.charge, axis=None), weight=events_temp.weight)
    
    
    def save_MET(self, events, out, ds, mode, tag=''):
        print("save_MET " + ds)
        out[f'pt_MET{tag}'].fill(ds=ds, pt=ak.flatten(events.MET.pt, axis=None), weight=events.weight)
        out[f'phi_MET{tag}'].fill(ds=ds, phi=ak.flatten(events.MET.phi, axis=None), weight=events.weight)


    
    def save_electron(self, events, out, ds, mode, tag=''):
        print("save_electron" + ds)
        lead_electron = events.SelElectron[:,0]
        out[f'pt_e_1{tag}'].fill(ds=ds, pt=ak.flatten(lead_electron.pt, axis=None), weight=events.weight)
        out[f'eta_e_1{tag}'].fill(ds=ds, eta=ak.flatten(lead_electron.eta, axis=None), weight=events.weight)
        out[f'phi_e_1{tag}'].fill(ds=ds, phi=ak.flatten(lead_electron.phi, axis=None), weight=events.weight)
        out[f'charge_e_1{tag}'].fill(ds=ds, charge=ak.flatten(lead_electron.charge, axis=None), weight=events.weight)
        out[f'mass_e_1{tag}'].fill(ds=ds, mass=ak.flatten(lead_electron.mass, axis=None), weight=events.weight)

        # events_temp = events[ak.num(events.SelElectron)>=2]
        # sublead_electron = events_temp.SelElectron[:,1]
        # out[f'pt_e_2{tag}'].fill(ds=ds, pt=ak.flatten(sublead_electron.pt, axis=None), weight=events_temp.weight)
        # out[f'eta_e_2{tag}'].fill(ds=ds, eta=ak.flatten(sublead_electron.eta, axis=None), weight=events_temp.weight)
        # out[f'phi_e_2{tag}'].fill(ds=ds, phi=ak.flatten(sublead_electron.phi, axis=None), weight=events_temp.weight)
        # out[f'charge_e_2{tag}'].fill(ds=ds, charge=ak.flatten(sublead_electron.charge, axis=None), weight=events_temp.weight)



    def save_tau(self, events, out, ds, mode, tag=''):
        print("save_tau" + ds)
        lead_tau = events.SelTau[:,0]


        out[f'pt_tau_1{tag}'].fill(ds=ds, pt=ak.flatten(lead_tau.pt, axis=None), weight=events.weight)
        out[f'eta_tau_1{tag}'].fill(ds=ds, eta=ak.flatten(lead_tau.eta, axis=None), weight=events.weight)
        out[f'phi_tau_1{tag}'].fill(ds=ds, phi=ak.flatten(lead_tau.phi, axis=None), weight=events.weight)
        out[f'charge_tau_1{tag}'].fill(ds=ds, charge=ak.flatten(lead_tau.charge, axis=None), weight=events.weight)
        out[f'mass_tau_1{tag}'].fill(ds=ds, mass=ak.flatten(lead_tau.mass, axis=None), weight=events.weight)

        # events_temp = events[ak.num(events.SelTau)>=2]
        # sublead_tau = events_temp.SelTau[:,1]

        # out[f'pt_tau_2{tag}'].fill(ds=ds, pt=ak.flatten(sublead_tau.pt, axis=None), weight=events_temp.weight)
        # out[f'eta_tau_2{tag}'].fill(ds=ds, eta=ak.flatten(sublead_tau.eta, axis=None), weight=events_temp.weight)
        # out[f'phi_tau_2{tag}'].fill(ds=ds, phi=ak.flatten(sublead_tau.phi, axis=None), weight=events_temp.weight)
        # out[f'charge_tau_2{tag}'].fill(ds=ds, charge=ak.flatten(sublead_tau.charge, axis=None), weight=events_temp.weight)

    




    
    def save_inv_mass(self, events, out, ds, mode, tag=''):
        print("save_inv_mass" + ds)
        events_muon_tau = events[ak.num(events.SelMuon) >= 1]
        events_muon_tau = events_muon_tau[ak.num(events_muon_tau.SelTau) >= 1]
        lead_muon = events_muon_tau.SelMuon[:, 0]
        lead_tau = events_muon_tau.SelTau[:,0]
        out[f'm_mu_tau{tag}'].fill(ds=ds, mass=(lead_muon+lead_tau).mass, weight=events_muon_tau.weight)
        out[f'dr_mu_tau{tag}'].fill(ds=ds, dr = delta_r(lead_tau,lead_muon), weight=events_muon_tau.weight)


        events_muon_electron = events[ak.num(events.SelMuon) >= 1]
        events_muon_electron = events_muon_electron[ak.num(events_muon_electron.SelElectron) >= 1]

        lead_electron = events_muon_electron.SelElectron[:, 0]
        lead_muon = events_muon_electron.SelMuon[:,0]
        out[f'm_e_mu{tag}'].fill(ds=ds, mass=(lead_electron+lead_muon).mass, weight=events_muon_electron.weight)
        out[f'dr_e_mu{tag}'].fill(ds=ds, dr = delta_r(lead_electron,lead_muon), weight=events_muon_electron.weight)


        events_electron_tau = events[ak.num(events.SelElectron) >= 1]
        events_electron_tau = events_electron_tau[ak.num(events_electron_tau.SelTau) >= 1]

        lead_tau = events_electron_tau.SelTau[:,0]
        lead_electron = events_electron_tau.SelElectron[:, 0]

        out[f'm_e_tau{tag}'].fill(ds=ds, mass=(lead_electron+lead_tau).mass, weight=events_electron_tau.weight)
        out[f'dr_e_tau{tag}'].fill(ds=ds, dr = delta_r(lead_electron,lead_tau), weight=events_electron_tau.weight)

    def save_3_highest_pt_mass(self, events, out, ds, mode, tag=''):
        print("save_3highest" + ds) 
        out[f'm_lll{tag}'].fill(ds=ds, mass=(events.SelLepton[:,0]+events.SelLepton[:,1] +events.SelLepton[:,2]).mass, weight=events.weight)


        


    #def fill_numbers(self, out, ds, events, label):
        # out[f'n_ev_{label}'][ds] += len(events)
        # out[f'sumw_{label}'][ds] += ak.sum(events.genWeight)

    def jet_veto_selection(self, events, dr_cut=0.5):
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
    # def jet_veto_selection(self, events, mode, dr_cut=0.5):
    #     jets = events.Jet[events.Jet.pt > 30.]
    #     lead_muon = events.MuonInvIso[:, 0] if 'inv_iso' in mode else events.SelMuon[:, 0]
    #     l1 = events.SelLepton[:,0]
    #     l2 = events.SelLepton[:,1]
    #     l3 = events.SelLepton[:,2]
    #     dr1 = jets.delta_r(l1)
    #     dr2 = jets.delta_r(l2)
    #     dr3 = jets.delta_r(l3)
    #     jets = jets[dr1 > 0.5]
    #     jets = jets[dr2 > 0.5]
    #     jets = jets[dr3 > 0.5]
    #     events  = events[ak.num(jets) == 0]
    #     ret
    def bjet_veto(self, events):
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
    
    def zero_out(self, events):
        events = events[events.run < 0]
        return events


    ### The following is example code

    # def fill_var_array(self, events, out, mode):
    #     # The following is potentially lots of duplication
    #     # Check if there's a good way how these can be added to "events"

    #     # Objects we need
        
    #     lead_muon = events.MuonInvIso[:, 0] if 'inv_iso' in mode else events.SelMuon[:, 0]


    #     # Calculate all interesting vars
    #     m = inv_mass_3p(lead_muon, dsa_1, dsa_2)
    #     dphi_1 = np.abs(delta_phi(dsa_1, lead_muon))
    #     dphi_2 = np.abs(delta_phi(dsa_2, lead_muon))
    #     dr = delta_r(dsa_1, dsa_2)
    #     weight = events.weight

    #     arrays = [ak.to_numpy(a, False) for a in (cos_oa, m, dxy_sig, sv_prob, dphi_1, dphi_2, dr, weight)]
    #     # Save vars to event
    #     arr = np.stack(arrays, axis=1)

    #     # Not sure if this is the best way, better add directly?
    #     out['sel_array'] += processor.column_accumulator(arr)

    


    def analyse_etautau(self, events, out, ds, mode):
        
 
        events = events[ak.num(events.SelElectron) ==  1] 
        events = events[ak.num(events.SelTau) ==  2]   
        events = events[ak.num(events.SelMuon) ==  0] 
        
        # Save histograms
        self.save_tag(events, out, ds, mode, '_etautau')


        ### Example - e.g. to write out array
        # self.fill_var_array(events, out, mode)

        # if 'blind' in mode:
        #     events = self.zero_out(events)

        # if ds in ['W', 'WSS']:
        #     events_check = events[dimuons[:, 0].mass < 4.]
        #     for ev in events_check:
        #         print(ds)
        #         print('run:', ev.run, ', lumi:', ev.luminosityBlock, ', ev:', ev.event)


        if len(events) == 0:
            return

        events = self.bjet_veto(events)
        self.fill_numbers(out, ds, events, 'bjet_veto')
        
        return events
    
    def analyse_mutautau(self, events, out, ds, mode):
        
    
            events = events[ak.num(events.SelElectron) ==  0] 
            events = events[ak.num(events.SelTau) ==  2]   
            events = events[ak.num(events.SelMuon) ==  1] 

            # Save histograms
            self.save_tag(events, out, ds, mode, '_mutautau')


            ### Example - e.g. to write out array
            # self.fill_var_array(events, out, mode)

            # if 'blind' in mode:
            #     events = self.zero_out(events)

            # if ds in ['W', 'WSS']:
            #     events_check = events[dimuons[:, 0].mass < 4.]
            #     for ev in events_check:
            #         print(ds)
            #         print('run:', ev.run, ', lumi:', ev.luminosityBlock, ', ev:', ev.event)


            if len(events) == 0:
                return

            events = self.bjet_veto(events)
            #self.fill_numbers(out, ds, events, 'bjet_veto')
            
            return events
        


    
    
    def analyse_emutau(self, events, out, ds, mode):
        events = events[ak.num(events.SelElectron) ==  1] 
        events = events[ak.num(events.SelTau) ==  1]   
        events = events[ak.num(events.SelMuon) ==  1]
        out['sumw_emutau'][ds] += ak.sum(events.genWeight)
        print("**********************")
        same_sign = (events.SelElectron[:,0].charge * events.SelMuon[:,0].charge) > 0
        events_SS = events[same_sign]
        events_OS = events[~same_sign]
        # Save histograms
        self.save_tag(events_SS, out, ds, mode, '_emutau_SS')
        self.save_tag(events_OS, out, ds, mode, '_emutau_OS')

        ### Example - e.g. to write out array
        # self.fill_var_array(events, out, mode)

        # if 'blind' in mode:
        #     events = self.zero_out(events)

        # if ds in ['W', 'WSS']:
        #     events_check = events[dimuons[:, 0].mass < 4.]
        #     for ev in events_check:
        #         print(ds)
        #         print('run:', ev.run, ', lumi:', ev.luminosityBlock, ', ev:', ev.event)


        if len(events) == 0:
            return

        #events = self.bjet_veto(events)
        #self.fill_numbers(out, ds, events, 'bjet_veto')
        
        return events
    

    def postprocess(self, accumulator):
        return accumulator

    def save_ratio_final_state(self, events, out, ds, mode, tag=''):


        taus = events.GenPart[(np.abs(events.GenPart.pdgId) == 15)]
        taus = taus[taus.hasFlags(['isLastCopy'])]
        counts = ak.num(taus)
        taus = ak.flatten(taus)
    
        taus_from_HNL = abs(taus.distinctParent.pdgId) == 9900012
        taus_from_W_from_HNL =  abs(taus.distinctParent.pdgId) ==24
        taus_from_W_from_q =  abs(taus.distinctParent.pdgId) <=6
        taus_from_W_from_g = abs(taus.distinctParent.pdgId) == 21
        taus_to_e = ak.any(abs(taus.distinctChildren.pdgId) == 11,axis = -1)
        taus_to_m = ak.any(abs(taus.distinctChildren.pdgId) == 13,axis = -1)
        taus_to_h = ak.any(abs(taus.distinctChildren.pdgId) > 40,axis = -1)

        taus_from_Z_or_hadron = ~(taus_from_HNL | taus_from_W_from_g | taus_from_W_from_q | taus_from_W_from_HNL)
        taus_from_Z_or_hadron_match =ak.any(ak.unflatten(taus_from_Z_or_hadron, counts), axis = -1)
        taus_from_Z_or_hadron = ak.count_nonzero(taus_from_Z_or_hadron_match)
    

        taus_from_HNL_to_e = taus_from_HNL & taus_to_e
        taus_from_HNL_to_e = ak.any(ak.unflatten(taus_from_HNL_to_e, counts), axis = -1)
        taus_from_HNL_to_mu = taus_from_HNL & taus_to_m
        taus_from_HNL_to_mu = ak.any(ak.unflatten(taus_from_HNL_to_mu, counts), axis = -1)
        taus_from_HNL_to_h = taus_from_HNL & taus_to_h
        taus_from_HNL_to_h = ak.any(ak.unflatten(taus_from_HNL_to_h, counts), axis = -1)
        taus_from_HNL = ak.any(ak.unflatten(taus_from_HNL, counts), axis = -1)
     
        taus_from_W_from_HNL_to_e = taus_from_W_from_HNL & taus_to_e
        taus_from_W_from_HNL_to_e = ak.any(ak.unflatten(taus_from_W_from_HNL_to_e, counts), axis = -1)
        taus_from_W_from_HNL_to_mu = taus_from_W_from_HNL & taus_to_m
        taus_from_W_from_HNL_to_mu = ak.any(ak.unflatten(taus_from_W_from_HNL_to_mu, counts), axis = -1)
        taus_from_W_from_HNL_to_h = taus_from_W_from_HNL & taus_to_h
        taus_from_W_from_HNL_to_h = ak.any(ak.unflatten(taus_from_W_from_HNL_to_h, counts), axis = -1)
        taus_from_W_from_HNL = ak.any(ak.unflatten(taus_from_W_from_HNL, counts), axis = -1)
     
        taus_from_W_from_q_to_e = taus_from_W_from_q & taus_to_e
        taus_from_W_from_q_to_e = ak.any(ak.unflatten(taus_from_W_from_q_to_e, counts), axis = -1)
        taus_from_W_from_q_to_mu = taus_from_W_from_q & taus_to_m
        taus_from_W_from_q_to_mu = ak.any(ak.unflatten(taus_from_W_from_q_to_mu, counts), axis = -1)
        taus_from_W_from_q_to_h = taus_from_W_from_q & taus_to_h
        taus_from_W_from_q_to_h = ak.any(ak.unflatten(taus_from_W_from_q_to_h, counts), axis = -1)
        taus_from_W_from_q = ak.any(ak.unflatten(taus_from_W_from_q, counts), axis = -1)

        taus_from_W_from_g_to_e = taus_from_W_from_g & taus_to_e
        taus_from_W_from_g_to_e = ak.any(ak.unflatten(taus_from_W_from_g_to_e, counts), axis = -1)
        taus_from_W_from_g_to_mu = taus_from_W_from_g & taus_to_m
        taus_from_W_from_g_to_mu = ak.any(ak.unflatten(taus_from_W_from_g_to_mu, counts), axis = -1)
        taus_from_W_from_g_to_h = taus_from_W_from_g & taus_to_h
        taus_from_W_from_g_to_h = ak.any(ak.unflatten(taus_from_W_from_g_to_h, counts), axis = -1)
        taus_from_W_from_g = ak.any(ak.unflatten(taus_from_W_from_g, counts), axis = -1)

        electrons = events.GenPart[(np.abs(events.GenPart.pdgId) == 11)]
        electrons = electrons[electrons.hasFlags(['isLastCopy'])]
        electrons_from_HNL = (ak.any(abs(electrons.distinctParent.pdgId) == 9900012, axis = -1))
        electrons_from_W_from_HNL = ( ak.any( abs(electrons.distinctParent.pdgId) ==24, axis = -1))
        electrons_from_W_from_q =  (( ak.any( abs(electrons.distinctParent.pdgId) <=6, axis = -1) |ak.any( abs(electrons.distinctParent.pdgId) ==21, axis = -1)  )& ak.any(abs(electrons.distinctParent.distinctChildren.pdgId) == 9900012)) 



        muons = events.GenPart[(np.abs(events.GenPart.pdgId) == 13)]
        muons = muons[muons.hasFlags(['isLastCopy'])]
        muons_from_HNL = (ak.any(abs(muons.distinctParent.pdgId) == 9900012, axis = -1))
        muons_from_W_from_HNL = ( ak.any( abs(muons.distinctParent.pdgId) ==24, axis = -1))
        muons_from_W_from_q =  (( ak.any( abs(muons.distinctParent.pdgId) <=6, axis = -1) |ak.any( abs(muons.distinctParent.pdgId) ==21, axis = -1)  )& ak.any(abs(muons.distinctParent.distinctChildren.pdgId) == 9900012)) 

        HNL = events.GenPart[(np.abs(events.GenPart.pdgId) == 9900012)]
        HNL = HNL[HNL.hasFlags(['isLastCopy'])]

        HNL_to_Z = ak.any(ak.any(abs(HNL.distinctChildren.pdgId) == 23, axis = -1), axis = -1)
        Znu = ak.count_nonzero(HNL_to_Z)
        others = ak.count_nonzero(HNL_to_Z | taus_from_Z_or_hadron_match)
        
        #ttt = ak.count_nonzero(taus_from_W_from_q & taus_from_W_from_HNL & taus_from_HNL)
        tetete = ak.count_nonzero(taus_from_W_from_q_to_e & taus_from_W_from_HNL_to_e & taus_from_HNL_to_e)
        tetetm = ak.count_nonzero(taus_from_W_from_q_to_e & taus_from_W_from_HNL_to_e & taus_from_HNL_to_mu)
        teteth = ak.count_nonzero(taus_from_W_from_q_to_e & taus_from_W_from_HNL_to_e & taus_from_HNL_to_h)

        tetmte = ak.count_nonzero(taus_from_W_from_q_to_e & taus_from_W_from_HNL_to_mu & taus_from_HNL_to_e)
        tetmtm = ak.count_nonzero(taus_from_W_from_q_to_e & taus_from_W_from_HNL_to_mu & taus_from_HNL_to_mu)
        tetmth = ak.count_nonzero(taus_from_W_from_q_to_e & taus_from_W_from_HNL_to_mu & taus_from_HNL_to_h)

        tethte = ak.count_nonzero(taus_from_W_from_q_to_e & taus_from_W_from_HNL_to_h & taus_from_HNL_to_e)
        tethtm = ak.count_nonzero(taus_from_W_from_q_to_e & taus_from_W_from_HNL_to_h & taus_from_HNL_to_mu)
        tethth = ak.count_nonzero(taus_from_W_from_q_to_e & taus_from_W_from_HNL_to_h & taus_from_HNL_to_h)
###############
        tmtete = ak.count_nonzero(taus_from_W_from_q_to_mu & taus_from_W_from_HNL_to_e & taus_from_HNL_to_e)
        tmtetm = ak.count_nonzero(taus_from_W_from_q_to_mu & taus_from_W_from_HNL_to_e & taus_from_HNL_to_mu)
        tmteth = ak.count_nonzero(taus_from_W_from_q_to_mu & taus_from_W_from_HNL_to_e & taus_from_HNL_to_h)

        tmtmte = ak.count_nonzero(taus_from_W_from_q_to_mu & taus_from_W_from_HNL_to_mu & taus_from_HNL_to_e)
        tmtmtm = ak.count_nonzero(taus_from_W_from_q_to_mu & taus_from_W_from_HNL_to_mu & taus_from_HNL_to_mu)
        tmtmth = ak.count_nonzero(taus_from_W_from_q_to_mu & taus_from_W_from_HNL_to_mu & taus_from_HNL_to_h)

        tmthte = ak.count_nonzero(taus_from_W_from_q_to_mu & taus_from_W_from_HNL_to_h & taus_from_HNL_to_e)
        tmthtm = ak.count_nonzero(taus_from_W_from_q_to_mu & taus_from_W_from_HNL_to_h & taus_from_HNL_to_mu)
        tmthth = ak.count_nonzero(taus_from_W_from_q_to_mu & taus_from_W_from_HNL_to_h & taus_from_HNL_to_h)

###############
        thtete = ak.count_nonzero(taus_from_W_from_q_to_h & taus_from_W_from_HNL_to_e & taus_from_HNL_to_e)
        thtetm = ak.count_nonzero(taus_from_W_from_q_to_h & taus_from_W_from_HNL_to_e & taus_from_HNL_to_mu)
        thteth = ak.count_nonzero(taus_from_W_from_q_to_h & taus_from_W_from_HNL_to_e & taus_from_HNL_to_h)

        thtmte = ak.count_nonzero(taus_from_W_from_q_to_h & taus_from_W_from_HNL_to_mu & taus_from_HNL_to_e)
        thtmtm = ak.count_nonzero(taus_from_W_from_q_to_h & taus_from_W_from_HNL_to_mu & taus_from_HNL_to_mu)
        thtmth = ak.count_nonzero(taus_from_W_from_q_to_h & taus_from_W_from_HNL_to_mu & taus_from_HNL_to_h)

        ththte = ak.count_nonzero(taus_from_W_from_q_to_h & taus_from_W_from_HNL_to_h & taus_from_HNL_to_e)
        ththtm = ak.count_nonzero(taus_from_W_from_q_to_h & taus_from_W_from_HNL_to_h & taus_from_HNL_to_mu)
        ththth = ak.count_nonzero(taus_from_W_from_q_to_h & taus_from_W_from_HNL_to_h & taus_from_HNL_to_h)
###########################

        #tte =  ak.count_nonzero(taus_from_W_from_q & taus_from_W_from_HNL & electrons_from_HNL)
        tetee =  ak.count_nonzero(taus_from_W_from_q_to_e & taus_from_W_from_HNL_to_e & electrons_from_HNL)
        tetme =  ak.count_nonzero(taus_from_W_from_q_to_e & taus_from_W_from_HNL_to_mu & electrons_from_HNL)
        tethe =  ak.count_nonzero(taus_from_W_from_q_to_e & taus_from_W_from_HNL_to_h & electrons_from_HNL)
    ####
        tmtee =  ak.count_nonzero(taus_from_W_from_q_to_mu & taus_from_W_from_HNL_to_e & electrons_from_HNL)
        tmtme =  ak.count_nonzero(taus_from_W_from_q_to_mu & taus_from_W_from_HNL_to_mu & electrons_from_HNL)
        tmthe =  ak.count_nonzero(taus_from_W_from_q_to_mu & taus_from_W_from_HNL_to_h & electrons_from_HNL)
        ####
        thtee =  ak.count_nonzero(taus_from_W_from_q_to_h & taus_from_W_from_HNL_to_e & electrons_from_HNL)
        thtme =  ak.count_nonzero(taus_from_W_from_q_to_h & taus_from_W_from_HNL_to_mu & electrons_from_HNL)
        ththe =  ak.count_nonzero(taus_from_W_from_q_to_h & taus_from_W_from_HNL_to_h & electrons_from_HNL)

###############################

        #ttm =  ak.count_nonzero(taus_from_W_from_q & taus_from_W_from_HNL & muons_from_HNL)
        tetem =  ak.count_nonzero(taus_from_W_from_q_to_e & taus_from_W_from_HNL_to_e & muons_from_HNL)
        tethm =  ak.count_nonzero(taus_from_W_from_q_to_e & taus_from_W_from_HNL_to_h & muons_from_HNL)
        tetmm =  ak.count_nonzero(taus_from_W_from_q_to_e & taus_from_W_from_HNL_to_mu & muons_from_HNL)
############
        tmtem =  ak.count_nonzero(taus_from_W_from_q_to_mu & taus_from_W_from_HNL_to_e & muons_from_HNL)
        tmthm =  ak.count_nonzero(taus_from_W_from_q_to_mu & taus_from_W_from_HNL_to_h & muons_from_HNL)
        tmtmm =  ak.count_nonzero(taus_from_W_from_q_to_mu & taus_from_W_from_HNL_to_mu & muons_from_HNL)
###########
        thtem =  ak.count_nonzero(taus_from_W_from_q_to_h & taus_from_W_from_HNL_to_e & muons_from_HNL)
        ththm =  ak.count_nonzero(taus_from_W_from_q_to_h & taus_from_W_from_HNL_to_h & muons_from_HNL)
        thtmm =  ak.count_nonzero(taus_from_W_from_q_to_h & taus_from_W_from_HNL_to_mu & muons_from_HNL)
##############3

        #tet =  ak.count_nonzero(taus_from_W_from_q & electrons_from_W_from_HNL & taus_from_HNL)
        teete =  ak.count_nonzero(taus_from_W_from_q_to_e & electrons_from_W_from_HNL & taus_from_HNL_to_e)
        teetm =  ak.count_nonzero(taus_from_W_from_q_to_e & electrons_from_W_from_HNL & taus_from_HNL_to_mu)
        teeth =  ak.count_nonzero(taus_from_W_from_q_to_e & electrons_from_W_from_HNL & taus_from_HNL_to_h)
        #################
        tmete =  ak.count_nonzero(taus_from_W_from_q_to_mu & electrons_from_W_from_HNL & taus_from_HNL_to_e)
        tmetm =  ak.count_nonzero(taus_from_W_from_q_to_mu & electrons_from_W_from_HNL & taus_from_HNL_to_mu)
        tmeth =  ak.count_nonzero(taus_from_W_from_q_to_mu & electrons_from_W_from_HNL & taus_from_HNL_to_h)
        #################
        thete =  ak.count_nonzero(taus_from_W_from_q_to_h & electrons_from_W_from_HNL & taus_from_HNL_to_e)
        thetm =  ak.count_nonzero(taus_from_W_from_q_to_h & electrons_from_W_from_HNL & taus_from_HNL_to_mu)
        theth =  ak.count_nonzero(taus_from_W_from_q_to_h & electrons_from_W_from_HNL & taus_from_HNL_to_h)
        #################
#################################

        #tee =  ak.count_nonzero(taus_from_W_from_q & electrons_from_W_from_HNL & electrons_from_HNL)
        teee =  ak.count_nonzero(taus_from_W_from_q_to_e & electrons_from_W_from_HNL & electrons_from_HNL)
        tmee =  ak.count_nonzero(taus_from_W_from_q_to_mu & electrons_from_W_from_HNL & electrons_from_HNL)
        thee =  ak.count_nonzero(taus_from_W_from_q_to_h & electrons_from_W_from_HNL & electrons_from_HNL)
#######################

        #tem =  ak.count_nonzero(taus_from_W_from_q & electrons_from_W_from_HNL & muons_from_HNL)
        teem =  ak.count_nonzero(taus_from_W_from_q_to_e & electrons_from_W_from_HNL & muons_from_HNL)
        tmem =  ak.count_nonzero(taus_from_W_from_q_to_mu & electrons_from_W_from_HNL & muons_from_HNL)
        them =  ak.count_nonzero(taus_from_W_from_q_to_h & electrons_from_W_from_HNL & muons_from_HNL)
 ##################
        #tmt =  ak.count_nonzero(taus_from_W_from_q & muons_from_W_from_HNL & taus_from_HNL)
        temte =  ak.count_nonzero(taus_from_W_from_q_to_e & muons_from_W_from_HNL & taus_from_HNL_to_e)
        temtm =  ak.count_nonzero(taus_from_W_from_q_to_e & muons_from_W_from_HNL & taus_from_HNL_to_mu)
        temth =  ak.count_nonzero(taus_from_W_from_q_to_e & muons_from_W_from_HNL & taus_from_HNL_to_h)
############################3
        tmmte =  ak.count_nonzero(taus_from_W_from_q_to_mu & muons_from_W_from_HNL & taus_from_HNL_to_e)
        tmmtm =  ak.count_nonzero(taus_from_W_from_q_to_mu & muons_from_W_from_HNL & taus_from_HNL_to_mu)
        tmmth =  ak.count_nonzero(taus_from_W_from_q_to_mu & muons_from_W_from_HNL & taus_from_HNL_to_h)
############################3

        thmte =  ak.count_nonzero(taus_from_W_from_q_to_h& muons_from_W_from_HNL & taus_from_HNL_to_e)
        thmtm =  ak.count_nonzero(taus_from_W_from_q_to_h & muons_from_W_from_HNL & taus_from_HNL_to_mu)
        thmth =  ak.count_nonzero(taus_from_W_from_q_to_h & muons_from_W_from_HNL & taus_from_HNL_to_h)
############################3

        #tme =  ak.count_nonzero(taus_from_W_from_q & muons_from_W_from_HNL & electrons_from_HNL)
        teme =  ak.count_nonzero(taus_from_W_from_q_to_e & muons_from_W_from_HNL & electrons_from_HNL)
        tmme =  ak.count_nonzero(taus_from_W_from_q_to_mu & muons_from_W_from_HNL & electrons_from_HNL)
        thme =  ak.count_nonzero(taus_from_W_from_q_to_h & muons_from_W_from_HNL & electrons_from_HNL)
####################################
        #tmm =  ak.count_nonzero(taus_from_W_from_q & muons_from_W_from_HNL & muons_from_HNL)
        temm =  ak.count_nonzero(taus_from_W_from_q_to_e & muons_from_W_from_HNL & muons_from_HNL)
        tmmm =  ak.count_nonzero(taus_from_W_from_q_to_mu & muons_from_W_from_HNL & muons_from_HNL)
        thmm =  ak.count_nonzero(taus_from_W_from_q_to_h & muons_from_W_from_HNL & muons_from_HNL)
################################
        #ett =  ak.count_nonzero(electrons_from_W_from_q & taus_from_W_from_HNL & taus_from_HNL)
        etete =  ak.count_nonzero(electrons_from_W_from_q & taus_from_W_from_HNL_to_e & taus_from_HNL_to_e)
        etetm =  ak.count_nonzero(electrons_from_W_from_q & taus_from_W_from_HNL_to_e & taus_from_HNL_to_mu)
        eteth =  ak.count_nonzero(electrons_from_W_from_q & taus_from_W_from_HNL_to_e & taus_from_HNL_to_h)
################
        etmte =  ak.count_nonzero(electrons_from_W_from_q & taus_from_W_from_HNL_to_mu & taus_from_HNL_to_e)
        etmtm =  ak.count_nonzero(electrons_from_W_from_q & taus_from_W_from_HNL_to_mu & taus_from_HNL_to_mu)
        etmth =  ak.count_nonzero(electrons_from_W_from_q & taus_from_W_from_HNL_to_mu & taus_from_HNL_to_h)
################
        ethte =  ak.count_nonzero(electrons_from_W_from_q & taus_from_W_from_HNL_to_h & taus_from_HNL_to_e)
        ethtm =  ak.count_nonzero(electrons_from_W_from_q & taus_from_W_from_HNL_to_h & taus_from_HNL_to_mu)
        ethth =  ak.count_nonzero(electrons_from_W_from_q & taus_from_W_from_HNL_to_h & taus_from_HNL_to_h)
################
        #ete =  ak.count_nonzero(electrons_from_W_from_q & taus_from_W_from_HNL & electrons_from_HNL)
        etee =  ak.count_nonzero(electrons_from_W_from_q & taus_from_W_from_HNL_to_e & electrons_from_HNL)
        etme =  ak.count_nonzero(electrons_from_W_from_q & taus_from_W_from_HNL_to_mu & electrons_from_HNL)
        ethe =  ak.count_nonzero(electrons_from_W_from_q & taus_from_W_from_HNL_to_h & electrons_from_HNL)

######################33
        #etm =  ak.count_nonzero(electrons_from_W_from_q & taus_from_W_from_HNL & muons_from_HNL)
        etem =  ak.count_nonzero(electrons_from_W_from_q & taus_from_W_from_HNL_to_e & muons_from_HNL)
        etmm =  ak.count_nonzero(electrons_from_W_from_q & taus_from_W_from_HNL_to_mu & muons_from_HNL)
        ethm =  ak.count_nonzero(electrons_from_W_from_q & taus_from_W_from_HNL_to_h & muons_from_HNL)
##############################

        #eet =  ak.count_nonzero(electrons_from_W_from_q & electrons_from_W_from_HNL & taus_from_HNL)
        eete =  ak.count_nonzero(electrons_from_W_from_q & electrons_from_W_from_HNL & taus_from_HNL_to_e)
        eetm =  ak.count_nonzero(electrons_from_W_from_q & electrons_from_W_from_HNL & taus_from_HNL_to_mu)
        eeth =  ak.count_nonzero(electrons_from_W_from_q & electrons_from_W_from_HNL & taus_from_HNL_to_h)

        ##########################
        eee =  ak.count_nonzero(electrons_from_W_from_q & electrons_from_W_from_HNL & electrons_from_HNL)
        eem =  ak.count_nonzero(electrons_from_W_from_q & electrons_from_W_from_HNL & muons_from_HNL)
######################33
        #emt =  ak.count_nonzero(electrons_from_W_from_q & muons_from_W_from_HNL & taus_from_HNL)
        emte =  ak.count_nonzero(electrons_from_W_from_q & muons_from_W_from_HNL & taus_from_HNL_to_e)
        emtm =  ak.count_nonzero(electrons_from_W_from_q & muons_from_W_from_HNL & taus_from_HNL_to_mu)
        emth =  ak.count_nonzero(electrons_from_W_from_q & muons_from_W_from_HNL & taus_from_HNL_to_h)

         ##################33
        eme =  ak.count_nonzero(electrons_from_W_from_q & muons_from_W_from_HNL & electrons_from_HNL)
        emm =  ak.count_nonzero(electrons_from_W_from_q & muons_from_W_from_HNL & muons_from_HNL)



#########################
        #mtt =  ak.count_nonzero(muons_from_W_from_q & taus_from_W_from_HNL & taus_from_HNL)
        mtete =  ak.count_nonzero(muons_from_W_from_q & taus_from_W_from_HNL_to_e & taus_from_HNL_to_e)
        mtetm =  ak.count_nonzero(muons_from_W_from_q & taus_from_W_from_HNL_to_e & taus_from_HNL_to_mu)
        mteth =  ak.count_nonzero(muons_from_W_from_q & taus_from_W_from_HNL_to_e & taus_from_HNL_to_h)
        #################3
        mtmte =  ak.count_nonzero(muons_from_W_from_q & taus_from_W_from_HNL_to_mu & taus_from_HNL_to_e)
        mtmtm =  ak.count_nonzero(muons_from_W_from_q & taus_from_W_from_HNL_to_mu & taus_from_HNL_to_mu)
        mtmth =  ak.count_nonzero(muons_from_W_from_q & taus_from_W_from_HNL_to_mu & taus_from_HNL_to_h)

        #################3
        mthte =  ak.count_nonzero(muons_from_W_from_q & taus_from_W_from_HNL_to_h & taus_from_HNL_to_e)
        mthtm =  ak.count_nonzero(muons_from_W_from_q & taus_from_W_from_HNL_to_h & taus_from_HNL_to_mu)
        mthth =  ak.count_nonzero(muons_from_W_from_q & taus_from_W_from_HNL_to_h & taus_from_HNL_to_h)

        #################3



        #mte =  ak.count_nonzero(muons_from_W_from_q & taus_from_W_from_HNL & electrons_from_HNL)
        mtee =  ak.count_nonzero(muons_from_W_from_q & taus_from_W_from_HNL_to_e & electrons_from_HNL)
        mtme =  ak.count_nonzero(muons_from_W_from_q & taus_from_W_from_HNL_to_mu & electrons_from_HNL)
        mthe =  ak.count_nonzero(muons_from_W_from_q & taus_from_W_from_HNL_to_h & electrons_from_HNL)

        ##########################
        #mtm =  ak.count_nonzero(muons_from_W_from_q & taus_from_W_from_HNL & muons_from_HNL)
        mtem =  ak.count_nonzero(muons_from_W_from_q & taus_from_W_from_HNL_to_e & muons_from_HNL)
        mtmm =  ak.count_nonzero(muons_from_W_from_q & taus_from_W_from_HNL_to_mu & muons_from_HNL)
        mthm =  ak.count_nonzero(muons_from_W_from_q & taus_from_W_from_HNL_to_h & muons_from_HNL)

        ####################3


        #met =  ak.count_nonzero(muons_from_W_from_q & electrons_from_W_from_HNL & taus_from_HNL)
        mete =  ak.count_nonzero(muons_from_W_from_q & electrons_from_W_from_HNL & taus_from_HNL_to_e)
        metm =  ak.count_nonzero(muons_from_W_from_q & electrons_from_W_from_HNL & taus_from_HNL_to_mu)
        meth =  ak.count_nonzero(muons_from_W_from_q & electrons_from_W_from_HNL & taus_from_HNL_to_h)

        ##################3
        mee =  ak.count_nonzero(muons_from_W_from_q & electrons_from_W_from_HNL & electrons_from_HNL)
        mem =  ak.count_nonzero(muons_from_W_from_q & electrons_from_W_from_HNL & muons_from_HNL)
####################3
        #mmt =  ak.count_nonzero(muons_from_W_from_q & muons_from_W_from_HNL & taus_from_HNL)
        mmte =  ak.count_nonzero(muons_from_W_from_q & muons_from_W_from_HNL & taus_from_HNL_to_e)
        mmtm =  ak.count_nonzero(muons_from_W_from_q & muons_from_W_from_HNL & taus_from_HNL_to_mu)
        mmth =  ak.count_nonzero(muons_from_W_from_q & muons_from_W_from_HNL & taus_from_HNL_to_h)
        ####################3
        mme =  ak.count_nonzero(muons_from_W_from_q & muons_from_W_from_HNL & electrons_from_HNL)
        mmm =  ak.count_nonzero(muons_from_W_from_q & muons_from_W_from_HNL & muons_from_HNL)


        total = len(events) - others
    

        h3 = ththth
        h2e1 = tethth + thteth + ththe + theth + ethth  +ththte 
        h2m1 = tmthth + ththm + thmth + mthth +ththtm +thtmth
        m2h1 = tmtmth + tmthtm + thtmtm + thmm + mthm + mmth + mthtm + mtmth + tmmth +thtmm + tmthm + thmtm
        e2h1= tethte + teteth + thtete + thee + ethe + eeth + ethte + eteth + thete + teeth + thtee +tethe
        h1m1e1 = tetmth + tethtm + tmteth + tmthte + thtetm + thtmte + tmthe + tethm + thtme + thtem + them + thme + ethm + mthe + meth + emth + mthte + mteth + ethtm + etmth + thmte + temth + thetm + tmeth

        e3 =eee + tetete + eete + etee + etete + teee + teete +tetee
        e2m1 = tetmte + tetetm + tmtete + tetme + tmtee + teetm + eem + eme + mee + teem + etem +tetem +mete + mtee + mtete + emte + eetm + etme + etmte + etetm + teme + temte + tmee + tmete
        m2e1 = tetmtm + tmtetm + tmtmte + tmtme + tetmm + tmem + emm + mem + mme + mmte + metm + mtem + mtme + mtmte + mtetm + emtm + etmm + etmtm + temm + tmme + tmmte + temtm + tmtem + tmetm
        m3 = tmtmtm + mmm + mmtm + mtmm + mtmtm + tmmm + tmmtm + tmtmm


        # t3 = ttt
        # t2e1 = tte + tet + ett
        # t2m1 = ttm + tmt + mtt
        # t1e1m1 = tme + tem + met + mte + etm + emt

        # e3 = eee
        # e2t1 = eet + ete + tee
        # e2m1 = eem = eme + mee

        # m3 = mmm
        # m2t1 = mmt + mtm + tmm
        # m2e1 = mme + mem + emm
        total_res = m3 + h3 + e3 + m2h1+ m2e1 + h2e1 + h2m1 + e2m1 + e2h1 + h1m1e1
        #total = total_res



        taus_to_e = ak.count_nonzero(ak.any(abs(taus.distinctChildren.pdgId) == 11, axis = -1))
        taus_to_m = ak.count_nonzero(ak.any(abs(taus.distinctChildren.pdgId) == 13, axis = -1))
        taus_to_h = ak.count_nonzero(ak.any(abs(taus.distinctChildren.pdgId) >40, axis = -1))
        taus_to_q = ak.count_nonzero(ak.any(abs(taus.distinctChildren.pdgId) <=6, axis = -1))

        total_decay_tau = taus_to_e  + taus_to_h + taus_to_m + taus_to_q
   


        nameoftxt = 'C:\\Users\\lucas\\Desktop\\PDM\\analysis\\ratio.txt' 
        with open(nameoftxt, 'a') as f:
            f.write("for " + str(ds) + ", ratio are : \n")
            f.write("h3 = " + str(h3)+"/"+str(total)+" = "+ str(h3/total)+" ----->" + str(100*h3/total) + "%\n")
            f.write("e3 = " + str(e3)+"/"+str(total)+" = "+ str(e3/total)+" ----->" + str(100*e3/total)+ "%\n")
            f.write("m3 = " + str(m3)+"/"+str(total)+" = "+ str(m3/total)+" ----->" + str(100*m3/total)+ "%\n")
            f.write("h2e1 = " + str(h2e1)+"/"+str(total)+" = "+ str(h2e1/total)+" ----->" + str(100*h2e1/total)+ "%\n")
            f.write("h2m1 = " + str(h2m1)+"/"+str(total)+" = "+ str(h2m1/total)+" ----->" + str(100*h2m1/total)+ "%\n")
            f.write("e2m1 = " + str(e2m1)+"/"+str(total)+" = "+ str(e2m1/total)+" ----->" + str(100*e2m1/total)+ "%\n")
            f.write("m2e1 = " + str(m2e1)+"/"+str(total)+" = "+ str(m2e1/total)+" ----->" + str(100*m2e1/total)+ "%\n")
            f.write("m2h1 = " + str(m2h1)+"/"+str(total)+" = "+ str(m2h1/total)+" ----->" + str(100*m2h1/total)+ "%\n")
            f.write("e2h1 = " + str(e2h1)+"/"+str(total)+" = "+ str(e2h1/total)+" ----->" + str(100*e2h1/total)+ "%\n")
            f.write("h1e1m1 = " + str(h1m1e1)+"/"+str(total)+" = "+ str(h1m1e1/total)+" ----->" + str(100*h1m1e1/total)+ "%\n")
            f.write("total = " + str(total_res)+"/"+str(total)+" = "+ str(total_res/total)+" ----->" + str(100*total_res/total)+ "%\n")
            f.write("\n")
            f.write("---------------------------------------------------------------------------------\n")
            f.write("\n")
            f.write("tau to muons = " + str(taus_to_m)+"/"+str(total_decay_tau) + " ----> " +str(100*taus_to_m/total_decay_tau)+"%\n")
            f.write("tau to electrons = " + str(taus_to_e)+"/"+str(total_decay_tau) + " ----> " +str(100*taus_to_e/total_decay_tau)+"%\n")
            f.write("tau to hadrons = " + str(taus_to_h + taus_to_q)+"/"+str(total_decay_tau) + " ----> " +str(100*(taus_to_h + taus_to_q)/total_decay_tau)+"%\n")
            f.write("\n")
            f.write("---------------------------------------------------------------------------------\n")
            f.write("\n")

       

