import pickle
import numpy as np
import uproot3
from coffea import hist
from samples import signal_samples
from helpers import files_from_dir, files_from_dirs

from ROOT import TFile



def HNL_plot_analysis(tag, rebining):
    # Read in reults from coffea run
    #tag = 'bck_WZ_ZZ_TT_tighter_sel_electron'
    #tag = '220505_bck_HNL_v1'

    with open(f'result_{tag}.pkl', 'rb') as f:
        result = pickle.load(f)
    with open (f'counter_{tag}.pkl', 'rb') as f:
        event_counter = pickle.load(f)

    print("Results are downloaded!")

    # Lumi and cross sections for plotting
    xsecs = { # pb
        #HNL
        #"HNL85": 1.,
        "HNL100":1.,
        "HNL125":1.,
        "HNL150":1.,
        "HNL200":1.,
        "HNL250":1.,
        "HNL300":1.,
        "HNL350":1.,
        "HNL400":1.,
        "HNL450":1.,
        "HNL500":1.,
        "HNL600":1.,
        "HNL700":1.,
        "HNL800":1.,
        "HNL900":1.,
        #"HNL1000":1.,

        ##Drell-Yann
        'DYJets_To_LL_M_50_madgraphMLM': 5398.0,
        #'DYJets_To_LL_M_50_amcatnloFXFX:': 6404.0,
        #'DY1Jets_To_LL_M_50' : 928.3,
        # 'DY2Jets_To_LL_M_50' : 293.6,
        # 'DY3Jets_To_LL_M_50' : 86.53,
        # 'DY4Jets_To_LL_M_50' : 41.28,
        #'DYJets_To_LL_M_10to50' : 15890.0,
        #Electroweak
        'EWK_WMinus2Jets_W_To_LNu_M_50':32.05,
        'EWK_WPlus2Jets_W_To_LNu_M_50':39.05,
        'EWK_Z2Jets_Z_To_LL_M_50':6.215,
        #TTbar
        'TT_To_2L2Nu':687.1,
        'TT_To_SemiLeptonic':687.1,
        'TT_To_Hadronic':687.1,
        # #SingleTop
        'ST_tW_antitop_5f_inclusiveDecays':34.97,
        'ST_tW_top_5f_inclusiveDecays':34.91,
        'ST_t_channel_antitop_4f_InclusiveDecays':69.09,
        'ST_t_channel_top_4f_InclusiveDecays':115.3,
        # W+jets
        #'WJets_To_LNu':53870.0,
        'WJ1ets_To_LNu':8927.0,
        'W2Jets_To_LNu':2809.0,
        'W3Jets_To_LNu':826.3,
        'W4Jets_To_LNu':544.3,
        #DiBoson
        #'ZZ':12.17,
        'ZZ_To_4L':1.325,
        'ZZ_To_2L2Nu':0.9738,
        #'WW':75.95,
        'WW_To_2L2Nu':11.09,
        #'WZ':27.59,
        'WZ_To_3LNu':5.213,
        #TT+bosons
        'TTWJets_To_LNu':0.2161,
        'TTZ_To_LLNuNu_M_10':0.2439,
        'TTWW':0.007003,
        'TTZZ':0.001386,
        'TTWZ':0.002453,


        #Tribosons
        "WWW_4F": 0.00001,
        "ZZZ": 0.00001,
        "WZZ": 0.00001,
        "WWZ_4F":0.00001,


        # #DATA
        # 'Data':1.,
        
        
        
        
        
        ###OLD ONES
        #'WZ':5.213,
        #'ZZ_To_4L':1.325,
        #'ZZ':12.14,
        #'ZZ_old':1.325, #ZZ to 4L old one
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
        # better to use this one : https://cms-gen-dev.cern.ch/xsdb/?columns=67108863&currentPage=0&pageSize=10&searchQuery=DAS%3DZZ_TuneCP5_13TeV-pythia8
        #'DY':6077.22,
        #'DY50':6077.22,
        #'DY10to50':6077.22,
        #'WZ_old':5.213,
        #'tt2l2nu':87.315,
        #'TTZ':87.315,
        #'W':20508.9*3.,
        #'TTWJets':20508.9*3.,

    }

    lumi = 60000. # pb-1
    plot_scale = {
        #HNL
        #"HNL85": 1.,
        "HNL100":10.,
        "HNL125":10.,
        "HNL150":10.,
        "HNL200":10.,
        "HNL250":10.,
        "HNL300":10.,
        "HNL350":10.,
        "HNL400":10.,
        "HNL450":10.,
        "HNL500":10.,
        "HNL600":10.,
        "HNL700":10.,
        "HNL800":10.,
        "HNL900":10.,
        #"HNL1000":1.,

        #Drell-Yann
        'DYJets_To_LL_M_50_madgraphMLM': 1.,
        'DYJets_To_LL_M_50_amcatnloFXFX:': 1.,
        'DY1Jets_To_LL_M_50' : 1.,
        'DY2Jets_To_LL_M_50' : 1.,
        'DY3Jets_To_LL_M_50' : 1.,
        'DY4Jets_To_LL_M_50' : 1.,
        'DYJets_To_LL_M_10to50' : 1.,
        #Electroweak
        'EWK_WMinus2Jets_W_To_LNu_M_50':1.,
        'EWK_WPlus2Jets_W_To_LNu_M_50':1.,
        'EWK_Z2Jets_Z_To_LL_M_50':1.,
        #TTbar
        'TT_To_2L2Nu':1.,
        'TT_To_SemiLeptonic':1.,
        'TT_To_Hadronic':1.,
        #SingleTop
        'ST_tW_antitop_5f_inclusiveDecays':1.,
        'ST_tW_top_5f_inclusiveDecays':1.,
        'ST_t_channel_antitop_4f_InclusiveDecays':1.,
        'ST_t_channel_top_4f_InclusiveDecays':1.,
        # W+jets
        'WJets_To_LNu':1.,
        'WJ1ets_To_LNu':1.,
        'W2Jets_To_LNu':1.,
        'W3Jets_To_LNu':1.,
        'W4Jets_To_LNu':1.,
        #DiBoson
        'ZZ':1.,
        'ZZ_To_4L':1.,
        'ZZ_To_2L2Nu':1.,
        'WW':1.,
        'WW_To_2L2Nu':1.,
        'WZ':1.,
        'WZ_To_3LNu':1.,
        #TT+bosons
        'TTWJets_To_LNu':1.,
        'TTZ_To_LLNuNu_M_10':1.,
        'TTWW':1.,
        'TTZZ':1.,
        'TTWZ':1.,
        # #DATA
        # 'Data':1.,
        #Tribosons
        "WWW_4F": 1.,
        "ZZZ": 1.,
        "WZZ": 1.,
        "WWZ_4F":1.,

        #     'WZ':1.,
        # 'WZ_old':1.,

        
        
    }
    event_counter['sumw']['Data'] = lumi

    with open('table_summary.txt', 'w') as f:
                f.write("var\tcross section\tN events expected\tSoW of runs tree\tSoW of full sample\tSoW after seletion\n")
    for s in xsecs.keys():
        print("##############################")
        print(str(s))
        print("#######################")
        exp = lumi*xsecs[s]
        all = result['sumw_all'][s]
        after_sel = result['sumw_emutau'][s]
        print(str(s)+"\t"+str(xsecs[s])+ "\t"+ str(exp) + "\t"+str(event_counter["sumw"][s])+"\t"+str(all)+"\t"+str(after_sel))

        with open('table_summary.txt', 'a') as f:
            
            f.write(str(s)+"\t"+str(xsecs[s])+ "\t"+ str(exp) + "\t"+str(event_counter["sumw"][s])+"\t"+str(all)+"\t"+str(after_sel)+"\n")
    

HNL_plot_analysis("220505_optimize_IDwp_run_13",10)    