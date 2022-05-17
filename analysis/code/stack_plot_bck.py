from datetime import datetime
import os
import re
import pickle

import numpy as np
import matplotlib.pyplot as plt
import uproot3
from coffea import hist
import time
from samples import signal_samples
from helpers import files_from_dir, files_from_dirs
#from exp_limits import exp_limits_mu

from HNLAnalysis import HNLAnalysis
from ROOT import TFile, TH1D, THStack,kGreen,kAzure, kCyan,kMagenta,kRed, kBlue, TCanvas, TLegend, gPad, gStyle, kRainBow, gSystem, kIsland, kGistEarth,TColor, kViridis,kOrange,kTeal, kPink,kSpring,kViolet
from cycler import cycler
import warnings

# Read in reults from coffea run
#tag = 'bck_WZ_ZZ_TT_tighter_sel_electron'
tag = '220428_all_bck'

with open(f'result_{tag}.pkl', 'rb') as f:
    result = pickle.load(f)
with open (f'counter_{tag}.pkl', 'rb') as f:
    event_counter = pickle.load(f)


print("Results are downloaded!")



# # Lumi and cross sections for plotting
# xsecs = { # pb
#     #HNL
#     # 'HNL100':1.,
#     # 'HNL500':1.,
#     # 'HNL1000':1.,
#     # 'HNL100_sing':1.,
#     # 'HNL500_sing':1.,
#     # 'HNL1000_sing':1.,
#     # 'HNL100_2highest':1.,
#     # 'HNL500_2highest':1.,
#     # 'HNL1000_2highest':1.,
#     # 'HNL100_+missingpT':1.,
#     # 'HNL500_+missingpT':1.,
#     # 'HNL1000_+missingpT':1.,
#     #Drell-Yann
#     'DYJets_To_LL_M_50_madgraphMLM': 5398.0,
#     'DYJets_To_LL_M_50_amcatnloFXFX: ': 6404.0,
#     'DYJets_To_LL_M_50_amcatnloFXFX:': 6404.0,
#     'DYJets_To_LL_M_50_amcatnloFXFX ': 6404.0,
#     'DYJets_To_LL_M_50_amcatnloFXFX': 6404.0,

#     'DY1Jets_To_LL_M_50' : 928.3,
#     'DY2Jets_To_LL_M_50' : 293.6,
#     'DY3Jets_To_LL_M_50' : 86.53,
#     'DY4Jets_To_LL_M_50' : 41.28,
#     'DYJets_To_LL_M_10to50 ' : 15890.0,
#     #Electroweak
#     'EWK_WMinus2Jets_W_To_LNu_M_50 ':32.05,
#     'EWK_WPlus2Jets_W_To_LNu_M_50  ':39.05,
#     'EWK_Z2Jets_Z_To_LL_M_50 ':6.215,
#     #TTbar
#     'TT_To_2L2Nu ':687.1,
#     'TT_To_SemiLeptonic ':687.1,
#     'TT_To_Hadronic ':687.1,
#     # #SingleTop
#     'ST_tW_antitop_5f_inclusiveDecays ':34.97,
#     'ST_tW_top_5f_inclusiveDecays ':34.91,
#     'ST_t_channel_antitop_4f_InclusiveDecays ':69.09,
#     'ST_t_channel_top_4f_InclusiveDecays ':115.3,
#     # W+jets
#     'WJets_To_LNu ':53870.0,
#     'WJ1ets_To_LNu ':8927.0,
#     'W2Jets_To_LNu ':2809.0,
#     'W3Jets_To_LNu ':826.3,
#     'W4Jets_To_LNu ':544.3,
#     #DiBoson
#     'ZZ ':12.17,
#     'ZZ_To_4L ':1.325,
#     'ZZ_To_2L2Nu ':0.9738,
#     'WW ':75.95,
#     'WW_To_2L2Nu ':11.09,
#     'WZ ':27.59,
#     'WZ_To_3LNu ':5.213,
#     #TT+bosons
#     'TTWJets_To_LNu ':0.2161,
#     'TTZ_To_LLNuNu_M_10 ':0.2439,
#     'TTWW  ':0.007003,
#     'TTZZ ':0.001386,
#     'TTWZ ':0.002453,
#     # #DATA
#     # 'Data':1.,
    
    
    
    
    
#     ###OLD ONES
#     #'WZ':5.213,
#     #'ZZ_To_4L':1.325,
#     #'ZZ':12.14,
#     #'ZZ_old':1.325, #ZZ to 4L old one
#     # https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
#     # better to use this one : https://cms-gen-dev.cern.ch/xsdb/?columns=67108863&currentPage=0&pageSize=10&searchQuery=DAS%3DZZ_TuneCP5_13TeV-pythia8
#     #'DY':6077.22,
#     #'DY50':6077.22,
#     #'DY10to50':6077.22,
#     #'WZ_old':5.213,
#     #'tt2l2nu':87.315,
#     #'TTZ':87.315,
#     #'W':20508.9*3.,
#     #'TTWJets':20508.9*3.,

# }

# lumi = 60000. # pb-1
# plot_scale = {
#         #HNL
#     # 'HNL100':1.,
#     # 'HNL500':1.,
#     # 'HNL1000':1.,
#     # 'HNL100_sing':1.,
#     # 'HNL500_sing':1.,
#     # 'HNL1000_sing':1.,
#     # 'HNL100_2highest':1.,
#     # 'HNL500_2highest':1.,
#     # 'HNL1000_2highest':1.,
#     # 'HNL100_+missingpT':1.,
#     # 'HNL500_+missingpT':1.,
#     # 'HNL1000_+missingpT':1.,
#     #Drell-Yann
#     'DYJets_To_LL_M_50_madgraphMLM': 1.,
#     'DYJets_To_LL_M_50_amcatnloFXFX: ': 1.0,
#     'DYJets_To_LL_M_50_amcatnloFXFX:': 1.0,
#     'DYJets_To_LL_M_50_amcatnloFXFX ': 1.0,
#     'DYJets_To_LL_M_50_amcatnloFXFX': 1.0,
#     'DY1Jets_To_LL_M_50' : 1.,
#     'DY2Jets_To_LL_M_50' : 1.,
#     'DY3Jets_To_LL_M_50' : 1.,
#     'DY4Jets_To_LL_M_50' : 1.,
#     'DYJets_To_LL_M_10to50 ' : 1.,
#     #Electroweak
#     'EWK_WMinus2Jets_W_To_LNu_M_50 ':1.,
#     'EWK_WPlus2Jets_W_To_LNu_M_50  ':1.,
#     'EWK_Z2Jets_Z_To_LL_M_50 ':1.,
#     #TTbar
#     'TT_To_2L2Nu ':1.,
#     'TT_To_SemiLeptonic ':1.,
#     'TT_To_Hadronic ':1.,
#     #SingleTop
#     'ST_tW_antitop_5f_inclusiveDecays ':1.,
#     'ST_tW_top_5f_inclusiveDecays ':1.,
#     'ST_t_channel_antitop_4f_InclusiveDecays ':1.,
#     'ST_t_channel_top_4f_InclusiveDecays ':1.,
#     # W+jets
#     'WJets_To_LNu ':1.,
#     'WJ1ets_To_LNu ':1.,
#     'W2Jets_To_LNu ':1.,
#     'W3Jets_To_LNu ':1.,
#     'W4Jets_To_LNu ':1.,
#     #DiBoson
#     'ZZ ':1.,
#     'ZZ_To_4L ':1.,
#     'ZZ_To_2L2Nu ':1.,
#     'WW ':1.,
#     'WW_To_2L2Nu ':1.,
#     'WZ ':1.,
#     'WZ_To_3LNu ':1.,
#     #TT+bosons
#     'TTWJets_To_LNu ':1.,
#     'TTZ_To_LLNuNu_M_10 ':1.,
#     'TTWW  ':1.,
#     'TTZZ ':1.,
#     'TTWZ ':1.,
#     # #DATA
#     # 'Data':1.,


#     #     'WZ':1.,
#     # 'WZ_old':1.,

    
    
# }




# Lumi and cross sections for plotting
xsecs = { # pb
    #HNL
    # 'HNL100':1.,
    # 'HNL500':1.,
    # 'HNL1000':1.,
    # 'HNL100_sing':1.,
    # 'HNL500_sing':1.,
    # 'HNL1000_sing':1.,
    # 'HNL100_2highest':1.,
    # 'HNL500_2highest':1.,
    # 'HNL1000_2highest':1.,
    # 'HNL100_+missingpT':1.,
    # 'HNL500_+missingpT':1.,
    # 'HNL1000_+missingpT':1.,
    #Drell-Yann
    'DYJets_To_LL_M_50_madgraphMLM': 5398.0,
    'DYJets_To_LL_M_50_amcatnloFXFX:': 6404.0,
    'DY1Jets_To_LL_M_50' : 928.3,
    'DY2Jets_To_LL_M_50' : 293.6,
    'DY3Jets_To_LL_M_50' : 86.53,
    'DY4Jets_To_LL_M_50' : 41.28,
    'DYJets_To_LL_M_10to50' : 15890.0,
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
    'WJets_To_LNu':53870.0,
    'WJ1ets_To_LNu':8927.0,
    'W2Jets_To_LNu':2809.0,
    'W3Jets_To_LNu':826.3,
    'W4Jets_To_LNu':544.3,
    #DiBoson
    'ZZ':12.17,
    'ZZ_To_4L':1.325,
    'ZZ_To_2L2Nu':0.9738,
    'WW':75.95,
    'WW_To_2L2Nu':11.09,
    'WZ':27.59,
    'WZ_To_3LNu':5.213,
    #TT+bosons
    'TTWJets_To_LNu':0.2161,
    'TTZ_To_LLNuNu_M_10':0.2439,
    'TTWW':0.007003,
    'TTZZ':0.001386,
    'TTWZ':0.002453,
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
    # 'HNL100':1.,
    # 'HNL500':1.,
    # 'HNL1000':1.,
    # 'HNL100_sing':1.,
    # 'HNL500_sing':1.,
    # 'HNL1000_sing':1.,
    # 'HNL100_2highest':1.,
    # 'HNL500_2highest':1.,
    # 'HNL1000_2highest':1.,
    # 'HNL100_+missingpT':1.,
    # 'HNL500_+missingpT':1.,
    # 'HNL1000_+missingpT':1.,
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


    #     'WZ':1.,
    # 'WZ_old':1.,

    
    
}
event_counter['sumw']['Data'] = lumi




not_hnl_data = re.compile('(?!HNL)(?!Data)^((?!SS).)*$')
hnl = re.compile('(HNL)')
# data = re.compile('(Data)((?!InvIso).)*$')
data = re.compile('^Data$|^DataSS$|^DataInvTimeSameBX$|^DataSSInvTimeSameBX$')

cols_bkg = ['#fde6f7', '#fafbc3', '#c4fbd6', '#b3e5f4']
cols_sig = ['black', '#386cb0','#f0027f','#bf5b17','#666666', '#7fc97f','#beaed4','#fdc086','#ffff99']
# cols_sig = ['red', 'green', 'blue','purple']
cols_stack = [kBlue-9,kMagenta-10,kRed-6, kOrange-2,kPink+6,kViolet-2, kSpring+2, kAzure-9, kRed-6]
plt.rcParams['axes.prop_cycle'] = cycler(color=cols_sig)
  
#gStyle.SetPalette(len(cols_stack),cols_stack)
gStyle.SetPalette(kViridis)
TColor.InvertPalette()
def number_of_counts(h):
    return h.Integral()

var_axis_pair = HNLAnalysis.get_var_axis_pairs()
for v in var_axis_pair:
    try:

        var = v[0]
        scales = {s:plot_scale[s]*lumi*xsecs[s]/event_counter['sumw'][s] if event_counter['sumw'][s] else 1. for s in xsecs.keys()}
        #print(scales)

        back_scales = {s:1./v for s, v in scales.items()}
        result[var].scale(scales, axis='ds')
        if not len(result[var].values()):
            continue
        print("variable " + str(var)+" is processed")

        plt.gca().set_prop_cycle(color=cols_sig)
        warnings.filterwarnings("ignore")
        
        #loop over all histogram of bck
        name_root_output = 'test_stack_'+ var+'.root'
        fout = uproot3.recreate(name_root_output )
        
        print("!!!WARNING!!! amcat is not plotted")
        for d in result[var].axis("ds").identifiers():
            if("amcat" not in str(d) and "Wjets_To_LNu" not in str(d) and "DY1Jets_To_LL_M_50" not in str(d) and "DY2Jets_To_LL_M_50" not in str(d) and "DY3Jets_To_LL_M_50" not in str(d) and "DY4Jets_To_LL_M_50" not in str(d) and "ZZ" != str(d)):
                new_h = result[var][d].integrate("ds")
                print(d)
                fout[str(d)] = hist.export1d(new_h)
           
        print("\tHistograms are exported to ROOT format")
        print("\tHistograms is written in root file named : " + name_root_output)


        fout.close()
        
        fout = TFile(name_root_output, "UPDATE")
        list_histos = []

        for d in result[var].axis("ds").identifiers():
            if("amcat" not in str(d) and "Wjets_To_LNu" not in str(d) and "DY1Jets_To_LL_M_50" not in str(d) and "DY2Jets_To_LL_M_50" not in str(d) and "DY3Jets_To_LL_M_50" not in str(d) and "DY4Jets_To_LL_M_50" not in str(d) and "ZZ" != str(d)):
                list_histos.append(fout.Get(str(d)))
        sum_count_histos = 0
        for h in list_histos:
            sum_count_histos += number_of_counts(h)
        if(sum_count_histos == 0):
            print("\t!!!!!!!")
            print("\tThis variable was not computed in the analysis, thus discarded")
            print("\tThe root file is deleted")
            fout.Close()
            gSystem.Exec("rm "+name_root_output)
            continue
        
        hs = THStack("stacked_"+var,"stacked_"+var)
        hs.SetTitle("stacked_"+var)
        c = TCanvas("stacked_"+var,"stacked_"+var)
        
        #sort the histograms form the biggest to the smallest one (number of counts)
        list_histos.sort(reverse=True, key=number_of_counts)
        print("\tHistograms are sorted by their number of total counts")
        iter_color = 0
        number_of_important_color = len(cols_stack)
        for h in list_histos:
            h.SetTitle(h.GetName())
            h.GetXaxis().SetTitle(var)
            h.GetYaxis().SetTitle("counts")
            if(iter_color<number_of_important_color):
                h.SetFillColor(cols_stack[iter_color%len(cols_stack)])
                h.SetLineColor(1)
            else:  
                h.SetFillColor(11+iter_color-number_of_important_color)
                h.SetLineColor(11+iter_color-number_of_important_color)  
            iter_color+=1

            hs.Add(h)
        print("\tHistograms are drawn on stacked histogram")

        #hs.GetXaxis().SetTitle(var)
        #hs.GetYaxis().SetName("counts")
        hs.Draw("HIST ") # PFC PLC PMC
        #gPad.SetGrid(1,0)
        c.Modified()

        c.BuildLegend(0.6,0.2,0.9,0.9,"")
        c.Draw()
        c.Write()
        print("\tCanvas is written in root file named : " + name_root_output)

        input("Press Enter to continue...")

        hs.Write()
        name_plot = tag +"_"+ var
        print('   ')
        result[var].scale(back_scales, axis='ds') #Alternatively make a deepcopy for each round of plotting?
        warnings.filterwarnings("default")
    except AssertionError:
        pass


