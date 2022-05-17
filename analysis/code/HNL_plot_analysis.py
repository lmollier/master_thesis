from datetime import datetime
import os
from pickletools import StackObject
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
from HNLAnalysis_etautau import HNLAnalysis_etautau
from HNLAnalysis_mutautau import HNLAnalysis_mutautau

from ROOT import TFile, TH1D, THStack,kGreen,kAzure, kCyan,kMagenta,kRed, kBlack,kBlue, TCanvas, TLegend, gPad, gStyle, kRainBow, gSystem, kIsland, kGistEarth,TColor, kViridis,kOrange,kTeal, kPink,kSpring,kViolet
from cycler import cycler
import warnings



def HNL_plot_analysis(tag, rebining):
    # Read in reults from coffea run
    #tag = 'bck_WZ_ZZ_TT_tighter_sel_electron'
    #tag = '220505_bck_HNL_v1'

    with open(f'result_{tag}.pkl', 'rb') as f:
        result = pickle.load(f)
    with open (f'counter_{tag}.pkl', 'rb') as f:
        event_counter = pickle.load(f)

    print("Results are downloaded!")
    name_root_output_final = 'histograms_output_'+tag+"_rebin_"+str(rebining)+'.root'
    fout_output = TFile(name_root_output_final, "RECREATE")
    fout_output.Close()
    print("root file for output is create under name : " + name_root_output_final)
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


    not_hnl_data = re.compile('(?!HNL)(?!Data)^((?!SS).)*$')
    hnl = re.compile('(HNL)')
    # data = re.compile('(Data)((?!InvIso).)*$')
    data = re.compile('^Data$|^DataSS$|^DataInvTimeSameBX$|^DataSSInvTimeSameBX$')

    cols_bkg = ['#fde6f7', '#fafbc3', '#c4fbd6', '#b3e5f4']
    #cols_sig = ['black', '#386cb0','#f0027f','#bf5b17','#666666', '#7fc97f','#beaed4','#fdc086','#ffff99']
    cols_sig = [kBlack,kAzure-3,kPink+7]
    cols_stack = [kBlue-9,kMagenta-10,kRed-6, kOrange-2,kPink+6,kViolet-2, kSpring+2, kAzure-9, kRed-6]
    #plt.rcParams['axes.prop_cycle'] = cycler(color=cols_sig)
    
    #gStyle.SetPalette(len(cols_stack),cols_stack)
    gStyle.SetPalette(kViridis)
    TColor.InvertPalette()
    def number_of_counts(h):
        return h.Integral()

    if("etautau" in tag):
        var_axis_pair = HNLAnalysis_etautau.get_var_axis_pairs()
    elif("mutautau" in tag):
        var_axis_pair = HNLAnalysis_mutautau.get_var_axis_pairs()
    else:
        var_axis_pair = HNLAnalysis.get_var_axis_pairs()


    nb_var = len(var_axis_pair)
    curr_var_i = 0

    
    for v in var_axis_pair:
        try:
            print("======== In process : "+str(int(100*curr_var_i/nb_var))+" %  ========")
            curr_var_i+=1
            var = v[0]
            
            scales = {s:plot_scale[s]*lumi*xsecs[s]/event_counter['sumw'][s] if event_counter['sumw'][s] else 1. for s in xsecs.keys()}
            #print(scales)

            back_scales = {s:1./v for s, v in scales.items()}
            result[var].scale(scales, axis='ds')
            if not len(result[var].values()):
                continue
            print("variable " + str(var)+" is processed:")

            #plt.gca().set_prop_cycle(color=cols_sig)
            warnings.filterwarnings("ignore")
            
            #loop over all histogram of bck
            name_root_output = 'bck_and_HNL_'+tag+"_"+ var+"_rebin_"+str(rebining)+'.root'
            fout = uproot3.recreate(name_root_output )
            print("\tStack of background is made and HNL root TH1D are made...")
            for d in result[var].axis("ds").identifiers():
                if("HNL" not in str(d)):
                    if("ZZ" != str(d) and "WZ" != str(d) and "WW" != str(d)):
                        h = result[var][d].integrate("ds")
                        fout[str(d)] = hist.export1d(h)
                else:
                    h = result[var][d].integrate("ds")
                    fout[str(d)] = hist.export1d(h)
            print("\tHistograms are exported to ROOT format")
            print("\tHistograms are written in root file named : " + name_root_output)
            fout.close()
            
            fout = TFile(name_root_output, "READ")
            list_histos = []
            list_histos_HNL = []

            for d in result[var].axis("ds").identifiers():
                if("HNL" not in str(d)):
                    if("ZZ" != str(d) and "WZ" != str(d) and "WW" != str(d)):
                        temp = fout.Get(str(d))
                        if(not ("charge" in var or "phi" in var or "dr" in var or "eta_" in var)):
                            temp = temp.Rebin(rebining)
                        elif( "phi" in var or "dr" in var or "eta_" in var):
                            temp = temp.Rebin(2)
                            print("rebin by 2 var:" + str(var))
                        list_histos.append(temp)
                else:
                    temp = fout.Get(str(d))
                    if(not ("charge" in var or "phi" in var or "dr" in var or "eta_" in var)):
                        temp = temp.Rebin(rebining)
                    elif( "phi" in var or "dr" in var or "eta_" in var):
                            temp = temp.Rebin(2)
                            print("rebin by 2 var:" + str(var))
                    list_histos_HNL.append(temp)
            
            
            sum_count_histos = 0
            for h in list_histos:
                sum_count_histos += number_of_counts(h)
            for h in list_histos_HNL:
                sum_count_histos += number_of_counts(h)    
            if(sum_count_histos == 0):
                print("\t!!!!!!!")
                print("\tThis variable was not computed in the analysis, thus discarded")
                print("\tThe root file is deleted")
                fout.Close()
                gSystem.Exec("rm "+name_root_output)
                continue


            hs = THStack("stacked_bck_"+var,"stacked_bck_"+var)
            hs.SetTitle("stacked_bck_"+var)



            #sort the histograms form the biggest to the smallest one (number of counts)
            list_histos.sort(reverse=True, key=number_of_counts)
            print("\tHistograms are sorted by their number of total counts")
            iter_color = 0
            number_of_important_color = len(cols_stack)
            
            fout_output = TFile(name_root_output_final, "UPDATE")


            fout_output.mkdir(var+"/raw histos/")
            fout_output.cd(var+"/raw histos/")
            first_rare = True
            first_red = True
            for h in list_histos:
                name = h.GetName()
                h.SetTitle(name+ "_"+var)
                h.SetName(name+ "_"+var)
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
                h.Write()
                cat = False
                cat_i = 0
                if("WZ_To_3LNu" == name):
                    h_WZ = h.Clone("h_WZ")
                    cat = True
                    cat_i += 1
                if("ZZ_To_4L" == name):
                    h_ZZ = h.Clone("h_ZZ")
                    cat = True
                    cat_i += 1
                if("TTWJets_To_LNu" == name or "TTZ_To_LLNuNu_M_10" == name or "TTWW" == name \
                    or "TTZZ" == name or "TTWZ" == name or "WWW_4F" == name or "ZZZ" == name or "WZZ" == name or "WWZ_4F" == name):
                    cat = True
                    cat_i += 1
                    if(first_rare):
                        first_rare = False
                        h_rare = h.Clone("h_rare")
                    else:
                        h_rare.Add(h)

                if("DYJets_To_LL_M_50_madgraphMLM" == name or "EWK_WMinus2Jets_W_To_LNu_M_50" == name  or "EWK_WPlus2Jets_W_To_LNu_M_50" == name \
                    or "EWK_Z2Jets_Z_To_LL_M_50" == name or "TT_To_2L2Nu" == name or "ST_tW_top_5f_inclusiveDecays" == name \
                    or "TT_To_SemiLeptonic" == name  or "TT_To_Hadronic" == name  or "ST_tW_antitop_5f_inclusiveDecays" == name \
                    or "ST_t_channel_antitop_4f_InclusiveDecays" == name or "ST_t_channel_top_4f_InclusiveDecays" == name or "WJ1ets_To_LNu" == name \
                    or "W2Jets_To_LNu" == name or "W3Jets_To_LNu" == name or "W4Jets_To_LNu" == name or "ZZ_To_2L2Nu" == name   \
                    or "WW_To_2L2Nu" == name or "W2Jets_To_LNu" == name or "W2Jets_To_LNu" == name or "W2Jets_To_LNu" == name or "W2Jets_To_LNu" == name):
                    cat = True
                    cat_i += 1
                    if(first_red):
                        first_red = False
                        h_reducible = h.Clone("h_reducible")
                    else:
                        h_reducible.Add(h)
                if(not cat):
                    print("this dataset was not categorized : "+name)
                if(cat_i > 1):
                    print("this dataset was categorized twice or more : "+name)
            fout_output.mkdir(var+"/output/")
            fout_output.cd(var+"/output/")
            print("\tHistogram of all bck together is stored ")
            h_all_bck = hs.GetStack().Last()
            h_all_bck.SetTitle("h_all_bck_"+var)
            h_all_bck.SetName("h_all_bck_"+var)
            h_all_bck.SetFillColor(kBlue)
            h_all_bck.SetLineColor(kBlue)

            h_all_bck.Write()

            print("\tHistograms are drawn on stacked histogram either global or individual")
            c = TCanvas("stacked_bck_"+var,"stacked_bck_"+var)
            hs.Draw("HIST ") 
            c.Modified()
            c.BuildLegend(0.6,0.2,0.9,0.9,"")
            c.Draw()
            c.Write()
            print("\tCanvas of stacked individual is written in root file named : " + name_root_output_final)
            c_all_bck = TCanvas("all_bck_"+var,"all_bck_"+var)
            h_all_bck.Draw("HIST ") 
            c_all_bck.Modified()
            c_all_bck.BuildLegend(0.7,0.7,0.9,0.9,"")
            c_all_bck.Draw()
            c_all_bck.Write()
            print("\tCanvas of all global bck is written in root file named : " + name_root_output_final)



    
            fout_output.cd(var+"/raw histos/")
            for h in list_histos_HNL:
                h.SetTitle(h.GetName()+ "_"+var)
                h.SetName(h.GetName()+ "_"+var)
                h.Write()
            print("\tHistograms of HNL is written in root file named : " + name_root_output_final)


            fout_output.cd(var+"/output/")
            c_HNL = TCanvas("comp_bck_HNL_"+var,"comp_bck_HNL_"+var)
            c_HNL.SetTitle("comp_bck_HNL_"+var)
            c_HNL.SetName("comp_bck_HNL_"+var)

            h_all_bck.SetFillColor(kBlue)
            h_all_bck.SetFillStyle(3003)
            #h_all_bck.SetLineColor(kBlue)
            h_all_bck.Draw("HIST") 
            iter_color = 0
            for h in list_histos_HNL:
                h.SetLineColor(cols_sig[iter_color%len(cols_sig)])
                h.SetLineWidth(2)
                temp =  str(h.GetName()).split("_")[0]
                h.SetTitle(temp)
                h.SetFillStyle(0)
                iter_color+=1
                h.Draw("HIST SAME")
            c_HNL.Modified()
            c_HNL.BuildLegend(0.7,0.7,0.9,0.9,"")
            c_HNL.Draw()
            c_HNL.Write()

            fout_output.mkdir("all_output/")
            fout_output.cd("all_output/")
            c_HNL.SetName(var)
            c_HNL.Write()
            print("\tHistograms of HNL vs. all bck is written in root file named : " + name_root_output_final)

        
            fout_output.cd(var+"/output/")
            c_comp_bcks_HNL = TCanvas("comp_bcks_HNL_"+var,"comp_bcks_HNL_"+var)
            c_comp_bcks_HNL.SetTitle("comp_bcks_HNL_"+var)
            c_comp_bcks_HNL.SetName("comp_bcks_HNL_"+var)

            h_ZZ.SetFillColor(kBlack)
            h_ZZ.SetFillStyle(3003)
            h_ZZ.SetName("ZZ_To_4L")
            h_ZZ.SetTitle("ZZ_To_4L")
            h_ZZ.Write()

            h_WZ.SetFillColor(kRed)
            h_WZ.SetLineColor(kRed)

            h_WZ.SetFillStyle(3003)
            h_WZ.SetName("WZ_To_3L")
            h_WZ.SetTitle("WZ_To_3L")
            h_WZ.Write()

            h_rare.SetFillColor(kGreen)
            h_rare.SetLineColor(kGreen)

            h_rare.SetFillStyle(3003)
            h_rare.SetName("Rare")
            h_rare.SetTitle("Rare")
            h_rare.Write()
            h_reducible.SetFillColor(kBlue)
            h_reducible.SetLineColor(kBlue)

            h_reducible.SetFillStyle(3003)  
            h_reducible.SetName("Reducible")
            h_reducible.SetTitle("Reducible")
            h_reducible.Write()
            stack_diff_bcks = THStack("stack_diff_bcks_"+var,"stack_diff_bcks_"+var)

            stack_diff_bcks.Add(h_reducible)
            stack_diff_bcks.Add(h_WZ)

            stack_diff_bcks.Add(h_ZZ)
            stack_diff_bcks.Add(h_rare)

            stack_diff_bcks.Draw("HIST")
            
            iter_color = 0
            for h in list_histos_HNL:
                h.SetLineColor(cols_sig[iter_color%len(cols_sig)])
                h.SetLineWidth(2)
                h.SetFillStyle(0)
                iter_color+=1
                h.Draw("HIST SAME")
            c_comp_bcks_HNL.Modified()
            c_comp_bcks_HNL.BuildLegend(0.7,0.7,0.9,0.9,"")
            c_comp_bcks_HNL.Draw()
            c_comp_bcks_HNL.Write()

            print("\tHistograms of HNL vs. stacked type of bck is written in root file named : " + name_root_output_final)


            fout_output.Close()
            gSystem.Exec("rm "+name_root_output)
            print("\tRoot file is deleted ,named : " + name_root_output)


            result[var].scale(back_scales, axis='ds') #Alternatively make a deepcopy for each round of plotting?
            warnings.filterwarnings("default")
        except AssertionError:
            pass


    fout_output = TFile(name_root_output_final, "UPDATE")

