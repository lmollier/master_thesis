from cmath import inf
from xxlimited import new
from ROOT import TFile, gSystem
import math
import array as arr
print("started")


path_to_results = "/afs/cern.ch/user/l/lmollier/git_PDM/analysis/results/"

#emutau_ss
input_file_emutau = path_to_results+ "histograms_output_220505_optimize_IDwp_run_12_rebin_5.root"
f_in = TFile(input_file_emutau, "READ")

output_filename_temp =  path_to_results+ "rootfile_limit_setting_temp_v2.root"
f_temp = TFile(output_filename_temp, "RECREATE")
f_temp.mkdir("emutau_SS/")
f_temp.cd("emutau_SS/")
f_temp.Close()

output_filename =  path_to_results+ "rootfile_limit_setting_v2.root"
f_out = TFile(output_filename, "RECREATE")
f_out.mkdir("emutau_SS/")
f_out.cd("emutau_SS/")
f_out.Close()



print(">>> emutau_SS")
variables_emutau = ["MT_total_emutau_SS","pt_sum_emutau_emutau_SS","comb_mass_mutau_emutau_SS", "mass_tau_1_emutau_SS"]
HNL_mass_list = [100,125,150,200,250,300,350,400,450,500,600,700,800,900]
f_temp = TFile(output_filename_temp, "UPDATE")
for var in variables_emutau:
    print(var)
    temp = var
    temp = temp.replace("_emutau_SS","")
   
    f_temp.mkdir("emutau_SS/"+temp+"/")
    f_temp.cd("emutau_SS/"+temp+"/")
    all_process = []
    all_process_for_binning = []
    name_histo_bck_rare = var+"/output/Rare"
    h_bck_rare= f_in.Get(name_histo_bck_rare)
    new_name = "Rare"
    h_bck_rare.SetName(new_name)
    h_bck_rare.SetTitle(new_name)
    curr_maxbin = h_bck_rare.FindLastBinAbove()
    cut = curr_maxbin
    all_process.append(new_name)
    h_bck_rare.SetBinContent(h_bck_rare.GetNbinsX(),h_bck_rare.GetBinContent(h_bck_rare.GetNbinsX())+h_bck_rare.GetBinContent(1+h_bck_rare.GetNbinsX())) #adding the overflow bin to the last bin
    h_bck_rare.SetBinContent(h_bck_rare.GetNbinsX()+1,0) #remove the overflow bin (=0)
    h_bck_rare.Write()
    #print("rare is written")
    name_histo_bck_reducible = var+"/output/Reducible"
    h_bck_reducible= f_in.Get(name_histo_bck_reducible)
    new_name = "Reducible"
    h_bck_reducible.SetName(new_name)
    h_bck_reducible.SetTitle(new_name)
    curr_maxbin = h_bck_reducible.FindLastBinAbove()
    cut = curr_maxbin
    all_process.append(new_name)
    h_bck_reducible.SetBinContent(h_bck_reducible.GetNbinsX(),h_bck_reducible.GetBinContent(h_bck_reducible.GetNbinsX())+h_bck_reducible.GetBinContent(1+h_bck_reducible.GetNbinsX())) #adding the overflow bin to the last bin
    h_bck_reducible.SetBinContent(h_bck_reducible.GetNbinsX()+1,0) #remove the overflow bin (=0)
    h_bck_reducible.Write()
    #print("red. is written")

    name_histo_data_obs = var+"/output/h_all_bck_"+var
    h_bck_data_obs= f_in.Get(name_histo_data_obs)
    new_name = "data_obs"
    h_bck_data_obs.SetName(new_name)
    h_bck_data_obs.SetTitle(new_name)
    curr_maxbin = h_bck_data_obs.FindLastBinAbove()
    cut = curr_maxbin
    all_process.append(new_name)
    h_bck_data_obs.SetBinContent(h_bck_data_obs.GetNbinsX(),h_bck_data_obs.GetBinContent(h_bck_data_obs.GetNbinsX())+h_bck_data_obs.GetBinContent(1+h_bck_data_obs.GetNbinsX())) #adding the overflow bin to the last bin
    h_bck_data_obs.SetBinContent(h_bck_data_obs.GetNbinsX()+1,0) #remove the overflow bin (=0)
    h_bck_data_obs.Write()
    #print("data obs is written") 
    all_process_for_binning.append(new_name)

    name_histo_bck_WZ = var+"/output/WZ_To_3L"
    h_bck_WZ= f_in.Get(name_histo_bck_WZ)
    new_name = "WZ"
    h_bck_WZ.SetName(new_name)
    h_bck_WZ.SetTitle(new_name)
    curr_maxbin = h_bck_WZ.FindLastBinAbove()
    cut = curr_maxbin
    all_process.append(new_name)
    h_bck_WZ.SetBinContent(h_bck_WZ.GetNbinsX(),h_bck_WZ.GetBinContent(h_bck_WZ.GetNbinsX())+h_bck_WZ.GetBinContent(1+h_bck_WZ.GetNbinsX())) #adding the overflow bin to the last bin
    h_bck_WZ.SetBinContent(h_bck_WZ.GetNbinsX()+1,0) #remove the overflow bin (=0)    
    h_bck_WZ.Write()
    #print("WZ is written")

    name_histo_bck_ZZ = var+"/output/ZZ_To_4L"
    h_bck_ZZ= f_in.Get(name_histo_bck_ZZ)
    new_name = "ZZ"
    h_bck_ZZ.SetName(new_name)
    h_bck_ZZ.SetTitle(new_name)
    curr_maxbin = h_bck_ZZ.FindLastBinAbove()
    cut = curr_maxbin
    all_process.append(new_name)
    h_bck_ZZ.SetBinContent(h_bck_ZZ.GetNbinsX(),h_bck_ZZ.GetBinContent(h_bck_ZZ.GetNbinsX())+h_bck_ZZ.GetBinContent(1+h_bck_ZZ.GetNbinsX())) #adding the overflow bin to the last bin
    h_bck_ZZ.SetBinContent(h_bck_ZZ.GetNbinsX()+1,0) #remove the overflow bin (=0)    
    h_bck_ZZ.Write()
    #print("ZZ is written")


    for HNL_mass in HNL_mass_list:
        HNL_mass = "HNL"+str(HNL_mass)
        name_histo_sig = var+"/raw histos/"+ HNL_mass + "_"+ var
        h_sig = f_in.Get(name_histo_sig)
        temp = var
        temp = temp.replace("_emutau_SS","")
        new_name = HNL_mass
        h_sig.SetName(new_name)
        h_sig.SetTitle(new_name)
        curr_maxbin = h_sig.FindLastBinAbove()
        if(curr_maxbin<cut):
            cut = curr_maxbin
        h_sig.SetBinContent(h_sig.GetNbinsX(),h_sig.GetBinContent(h_sig.GetNbinsX())+h_sig.GetBinContent(1+h_sig.GetNbinsX())) #adding the overflow bin to the last bin
        h_sig.SetBinContent(h_sig.GetNbinsX()+1,0) #remove the overflow bin (=0)  
        h_sig.Write()
        #print("one sig is written")

        all_process.append(new_name)
        #all_process_for_binning.append(new_name)
    
    
    

    # for process in all_process:
    #     h = f_temp.Get("emutau_SS/"+temp+"/"+process)
    #     nx = cut
    #     xmin = h.GetBinLowEdge(1)
    #     xmax = h.GetBinLowEdge(cut)+h.GetBinWidth(cut)
    #     h.SetBins(nx, xmin, xmax)
    #     h.Write()
    # print("histograms are cut")
    print("The histograms are not cut")
    i = h_sig.GetNbinsX()
    # f_out = TFile(output_filename, "UPDATE")
    # f_out.cd("emutau_SS/")
    while(i>1):
        redo_binning = False
        for process in all_process_for_binning:
            h = f_temp.Get("emutau_SS/"+temp+"/"+process)
            
            if(h.GetBinContent(i) < 5):
                redo_binning = True
        if(redo_binning):
            nbins_old = h.GetNbinsX()
            old_binning = []
            new_binning = arr.array('d')
            for bin in range(1, nbins_old+2):
                old_binning.append(h.GetBinLowEdge(bin))
                if(bin != i):
                    new_binning .append(h.GetBinLowEdge(bin))
            
            for process in all_process:
                h = f_temp.Get("emutau_SS/"+temp+"/"+process)
                hnew = h.Rebin(len(new_binning)-1, "hnew", new_binning)
                h.Delete()
                hnew.SetTitle(process)
                hnew.SetName(process)
                hnew.Write()

    
        i = i - 1
    
    f_out = TFile(output_filename, "UPDATE")
    f_out.mkdir("emutau_SS/"+temp+"/")
    f_out.cd("emutau_SS/"+temp+"/")
    for process in all_process:
            h = f_temp.Get("emutau_SS/"+temp+"/"+process)
            h.Write()
    f_out.Close()  
f_temp.Close()  


















print(">>> emutau_OS")

#emutau_OS
f_in = TFile(input_file_emutau, "READ")

f_temp = TFile(output_filename_temp, "UPDATE")
f_temp.mkdir("emutau_OS/")
f_temp.cd("emutau_OS/")
f_temp.Close()

f_out = TFile(output_filename, "UPDATE")
f_out.mkdir("emutau_OS/")
f_out.cd("emutau_OS/")
f_out.Close()



variables_emutau = ["transverse_mass_etau_emutau_OS","dr_emu_emutau_OS","MT_total_emutau_OS","pt_tau_1_emutau_OS"]

HNL_mass_list = [100,125,150,200,250,300,350,400,450,500,600,700,800,900]
f_temp = TFile(output_filename_temp, "UPDATE")
for var in variables_emutau:
    print(var)
    temp = var
    temp = temp.replace("_emutau_OS","")
   
    f_temp.mkdir("emutau_OS/"+temp+"/")
    f_temp.cd("emutau_OS/"+temp+"/")
    all_process = []
    all_process_for_binning = []
    name_histo_bck_rare = var+"/output/Rare"
    h_bck_rare= f_in.Get(name_histo_bck_rare)
    new_name = "Rare"
    h_bck_rare.SetName(new_name)
    h_bck_rare.SetTitle(new_name)
    curr_maxbin = h_bck_rare.FindLastBinAbove()
    cut = curr_maxbin
    all_process.append(new_name)
    h_bck_rare.Write()
    #print("rare is written")
    name_histo_bck_reducible = var+"/output/Reducible"
    h_bck_reducible= f_in.Get(name_histo_bck_reducible)
    new_name = "Reducible"
    h_bck_reducible.SetName(new_name)
    h_bck_reducible.SetTitle(new_name)
    curr_maxbin = h_bck_reducible.FindLastBinAbove()
    cut = curr_maxbin
    all_process.append(new_name)
    h_bck_reducible.Write()
    #print("red. is written")

    name_histo_data_obs = var+"/output/h_all_bck_"+var
    h_bck_data_obs= f_in.Get(name_histo_data_obs)
    new_name = "data_obs"
    h_bck_data_obs.SetName(new_name)
    h_bck_data_obs.SetTitle(new_name)
    curr_maxbin = h_bck_data_obs.FindLastBinAbove()
    cut = curr_maxbin
    all_process.append(new_name)
    h_bck_data_obs.Write()
    #print("data obs is written") 
    all_process_for_binning.append(new_name)

    name_histo_bck_WZ = var+"/output/WZ_To_3L"
    h_bck_WZ= f_in.Get(name_histo_bck_WZ)
    new_name = "WZ"
    h_bck_WZ.SetName(new_name)
    h_bck_WZ.SetTitle(new_name)
    curr_maxbin = h_bck_WZ.FindLastBinAbove()
    cut = curr_maxbin
    all_process.append(new_name)
    h_bck_WZ.Write()
    #print("WZ is written")

    name_histo_bck_ZZ = var+"/output/ZZ_To_4L"
    h_bck_ZZ= f_in.Get(name_histo_bck_ZZ)
    new_name = "ZZ"
    h_bck_ZZ.SetName(new_name)
    h_bck_ZZ.SetTitle(new_name)
    curr_maxbin = h_bck_ZZ.FindLastBinAbove()
    cut = curr_maxbin
    all_process.append(new_name)
    h_bck_ZZ.Write()
    #print("ZZ is written")


    for HNL_mass in HNL_mass_list:
        HNL_mass = "HNL"+str(HNL_mass)
        name_histo_sig = var+"/raw histos/"+ HNL_mass + "_"+ var
        h_sig = f_in.Get(name_histo_sig)
        temp = var
        temp = temp.replace("_emutau_OS","")
        new_name = HNL_mass
        h_sig.SetName(new_name)
        h_sig.SetTitle(new_name)
        curr_maxbin = h_sig.FindLastBinAbove()
        if(curr_maxbin<cut):
            cut = curr_maxbin
        h_sig.Write()
        #print("one sig is written")

        all_process.append(new_name)
        #all_process_for_binning.append(new_name)
    
    
    

    # for process in all_process:
    #     h = f_temp.Get("emutau_SS/"+temp+"/"+process)
    #     nx = cut
    #     xmin = h.GetBinLowEdge(1)
    #     xmax = h.GetBinLowEdge(cut)+h.GetBinWidth(cut)
    #     h.SetBins(nx, xmin, xmax)
    #     h.Write()
    # print("histograms are cut")
    print("The histograms are not cut")
    i = h_sig.GetNbinsX()
    # f_out = TFile(output_filename, "UPDATE")
    # f_out.cd("emutau_SS/")
    while(i>1):
        redo_binning = False
        for process in all_process_for_binning:
            h = f_temp.Get("emutau_OS/"+temp+"/"+process)
            
            if(h.GetBinContent(i) < 5):
                redo_binning = True
        if(redo_binning):
            nbins_old = h.GetNbinsX()
            old_binning = []
            new_binning = arr.array('d')
            for bin in range(1, nbins_old+2):
                old_binning.append(h.GetBinLowEdge(bin))
                if(bin != i):
                    new_binning .append(h.GetBinLowEdge(bin))
            
            for process in all_process:
                h = f_temp.Get("emutau_OS/"+temp+"/"+process)
                hnew = h.Rebin(len(new_binning)-1, "hnew", new_binning)
                h.Delete()
                hnew.SetTitle(process)
                hnew.SetName(process)
                hnew.Write()

    
        i = i - 1
    
    f_out = TFile(output_filename, "UPDATE")
    f_out.mkdir("emutau_OS/"+temp+"/")
    f_out.cd("emutau_OS/"+temp+"/")
    for process in all_process:
            h = f_temp.Get("emutau_OS/"+temp+"/"+process)
            h.Write()
    f_out.Close()  
f_temp.Close()  






print(">>> etautau")

#etautau
input_file_etautau =  path_to_results+  "histograms_output_220509_optimize_IDwp_etautau_run_9_rebin_5.root"
f_in = TFile(input_file_etautau, "READ")

f_temp = TFile(output_filename_temp, "UPDATE")
f_temp.mkdir("etautau/")
f_temp.cd("etautau/")
f_temp.Close()

f_out = TFile(output_filename, "UPDATE")
f_out.mkdir("etautau/")
f_out.cd("etautau/")
f_out.Close()



variables_etautau = ["pt_sum_etauMET_etautau","transverse_mass_tau2tau_etautau","dr_etau2_etautau","MT_total_etautau"]
HNL_mass_list = [100,125,150,200,250,300,350,400,450,500,600,700,800,900]
f_temp = TFile(output_filename_temp, "UPDATE")
for var in variables_etautau:
    print(var)
    temp = var
    temp = temp.replace("_etautau","")
   
    f_temp.mkdir("etautau/"+temp+"/")
    f_temp.cd("etautau/"+temp+"/")
    all_process = []
    all_process_for_binning = []
    name_histo_bck_rare = var+"/output/Rare"
    h_bck_rare= f_in.Get(name_histo_bck_rare)
    new_name = "Rare"
    h_bck_rare.SetName(new_name)
    h_bck_rare.SetTitle(new_name)
    curr_maxbin = h_bck_rare.FindLastBinAbove()
    cut = curr_maxbin
    all_process.append(new_name)
    h_bck_rare.Write()
    #print("rare is written")
    name_histo_bck_reducible = var+"/output/Reducible"
    h_bck_reducible= f_in.Get(name_histo_bck_reducible)
    new_name = "Reducible"
    h_bck_reducible.SetName(new_name)
    h_bck_reducible.SetTitle(new_name)
    curr_maxbin = h_bck_reducible.FindLastBinAbove()
    cut = curr_maxbin
    all_process.append(new_name)
    h_bck_reducible.Write()
    #print("red. is written")

    name_histo_data_obs = var+"/output/h_all_bck_"+var
    h_bck_data_obs= f_in.Get(name_histo_data_obs)
    new_name = "data_obs"
    h_bck_data_obs.SetName(new_name)
    h_bck_data_obs.SetTitle(new_name)
    curr_maxbin = h_bck_data_obs.FindLastBinAbove()
    cut = curr_maxbin
    all_process.append(new_name)
    h_bck_data_obs.Write()
    #print("data obs is written") 
    all_process_for_binning.append(new_name)

    name_histo_bck_WZ = var+"/output/WZ_To_3L"
    h_bck_WZ= f_in.Get(name_histo_bck_WZ)
    new_name = "WZ"
    h_bck_WZ.SetName(new_name)
    h_bck_WZ.SetTitle(new_name)
    curr_maxbin = h_bck_WZ.FindLastBinAbove()
    cut = curr_maxbin
    all_process.append(new_name)
    h_bck_WZ.Write()
    #print("WZ is written")

    name_histo_bck_ZZ = var+"/output/ZZ_To_4L"
    h_bck_ZZ= f_in.Get(name_histo_bck_ZZ)
    new_name = "ZZ"
    h_bck_ZZ.SetName(new_name)
    h_bck_ZZ.SetTitle(new_name)
    curr_maxbin = h_bck_ZZ.FindLastBinAbove()
    cut = curr_maxbin
    all_process.append(new_name)
    h_bck_ZZ.Write()
    #print("ZZ is written")


    for HNL_mass in HNL_mass_list:
        HNL_mass = "HNL"+str(HNL_mass)
        name_histo_sig = var+"/raw histos/"+ HNL_mass + "_"+ var
        h_sig = f_in.Get(name_histo_sig)
        temp = var
        temp = temp.replace("_etautau","")
        new_name = HNL_mass
        h_sig.SetName(new_name)
        h_sig.SetTitle(new_name)
        curr_maxbin = h_sig.FindLastBinAbove()
        if(curr_maxbin<cut):
            cut = curr_maxbin
        h_sig.Write()
        #print("one sig is written")

        all_process.append(new_name)
        #all_process_for_binning.append(new_name)
    
    
    

    # for process in all_process:
    #     h = f_temp.Get("emutau_SS/"+temp+"/"+process)
    #     nx = cut
    #     xmin = h.GetBinLowEdge(1)
    #     xmax = h.GetBinLowEdge(cut)+h.GetBinWidth(cut)
    #     h.SetBins(nx, xmin, xmax)
    #     h.Write()
    # print("histograms are cut")
    print("The histograms are not cut")
    i = h_sig.GetNbinsX()
    # f_out = TFile(output_filename, "UPDATE")
    # f_out.cd("emutau_SS/")
    while(i>1):
        redo_binning = False
        for process in all_process_for_binning:
            h = f_temp.Get("etautau/"+temp+"/"+process)
            
            if(h.GetBinContent(i) < 5):
                redo_binning = True
        if(redo_binning):
            nbins_old = h.GetNbinsX()
            old_binning = []
            new_binning = arr.array('d')
            for bin in range(1, nbins_old+2):
                old_binning.append(h.GetBinLowEdge(bin))
                if(bin != i):
                    new_binning .append(h.GetBinLowEdge(bin))
            
            for process in all_process:
                h = f_temp.Get("etautau/"+temp+"/"+process)
                hnew = h.Rebin(len(new_binning)-1, "hnew", new_binning)
                h.Delete()
                hnew.SetTitle(process)
                hnew.SetName(process)
                hnew.Write()

    
        i = i - 1
    
    f_out = TFile(output_filename, "UPDATE")
    f_out.mkdir("etautau/"+temp+"/")
    f_out.cd("etautau/"+temp+"/")
    for process in all_process:
            h = f_temp.Get("etautau/"+temp+"/"+process)
            h.Write()
    f_out.Close()  
f_temp.Close()  



print(">>> mutautau")
#mutautau
input_file_mutautau =  path_to_results+  "histograms_output_220510_optimize_IDwp_mutautau_run_3_rebin_5.root"
f_in = TFile(input_file_mutautau, "READ")

f_temp = TFile(output_filename_temp, "UPDATE")
f_temp.mkdir("mutautau/")
f_temp.cd("mutautau/")
f_temp.Close()

f_out = TFile(output_filename, "UPDATE")
f_out.mkdir("mutautau/")
f_out.cd("mutautau/")
f_out.Close()



variables_mutautau = ["MT_total_mutautau","transverse_mass_tau2tau_mutautau","dr_mutau2_mutautau","pt_sum_mutau2tau_mutautau"]
HNL_mass_list = [100,125,150,200,250,300,350,400,450,500,600,700,800,900]
f_temp = TFile(output_filename_temp, "UPDATE")
for var in variables_mutautau:
    print(var)
    temp = var
    temp = temp.replace("_mutautau","")
   
    f_temp.mkdir("mutautau/"+temp+"/")
    f_temp.cd("mutautau/"+temp+"/")
    all_process = []
    all_process_for_binning = []
    name_histo_bck_rare = var+"/output/Rare"
    h_bck_rare= f_in.Get(name_histo_bck_rare)
    new_name = "Rare"
    h_bck_rare.SetName(new_name)
    h_bck_rare.SetTitle(new_name)
    curr_maxbin = h_bck_rare.FindLastBinAbove()
    cut = curr_maxbin
    all_process.append(new_name)
    h_bck_rare.Write()
    #print("rare is written")
    name_histo_bck_reducible = var+"/output/Reducible"
    h_bck_reducible= f_in.Get(name_histo_bck_reducible)
    new_name = "Reducible"
    h_bck_reducible.SetName(new_name)
    h_bck_reducible.SetTitle(new_name)
    curr_maxbin = h_bck_reducible.FindLastBinAbove()
    cut = curr_maxbin
    all_process.append(new_name)
    h_bck_reducible.Write()
    #print("red. is written")

    name_histo_data_obs = var+"/output/h_all_bck_"+var
    h_bck_data_obs= f_in.Get(name_histo_data_obs)
    new_name = "data_obs"
    h_bck_data_obs.SetName(new_name)
    h_bck_data_obs.SetTitle(new_name)
    curr_maxbin = h_bck_data_obs.FindLastBinAbove()
    cut = curr_maxbin
    all_process.append(new_name)
    h_bck_data_obs.Write()
    #print("data obs is written") 
    all_process_for_binning.append(new_name)

    name_histo_bck_WZ = var+"/output/WZ_To_3L"
    h_bck_WZ= f_in.Get(name_histo_bck_WZ)
    new_name = "WZ"
    h_bck_WZ.SetName(new_name)
    h_bck_WZ.SetTitle(new_name)
    curr_maxbin = h_bck_WZ.FindLastBinAbove()
    cut = curr_maxbin
    all_process.append(new_name)
    h_bck_WZ.Write()
    #print("WZ is written")

    name_histo_bck_ZZ = var+"/output/ZZ_To_4L"
    h_bck_ZZ= f_in.Get(name_histo_bck_ZZ)
    new_name = "ZZ"
    h_bck_ZZ.SetName(new_name)
    h_bck_ZZ.SetTitle(new_name)
    curr_maxbin = h_bck_ZZ.FindLastBinAbove()
    cut = curr_maxbin
    all_process.append(new_name)
    h_bck_ZZ.Write()
    #print("ZZ is written")


    for HNL_mass in HNL_mass_list:
        HNL_mass = "HNL"+str(HNL_mass)
        name_histo_sig = var+"/raw histos/"+ HNL_mass + "_"+ var
        h_sig = f_in.Get(name_histo_sig)
        temp = var
        temp = temp.replace("_mutautau","")
        new_name = HNL_mass
        h_sig.SetName(new_name)
        h_sig.SetTitle(new_name)
        curr_maxbin = h_sig.FindLastBinAbove()
        if(curr_maxbin<cut):
            cut = curr_maxbin
        h_sig.Write()
        #print("one sig is written")

        all_process.append(new_name)
        #all_process_for_binning.append(new_name)
    
    
    

    # for process in all_process:
    #     h = f_temp.Get("emutau_SS/"+temp+"/"+process)
    #     nx = cut
    #     xmin = h.GetBinLowEdge(1)
    #     xmax = h.GetBinLowEdge(cut)+h.GetBinWidth(cut)
    #     h.SetBins(nx, xmin, xmax)
    #     h.Write()
    # print("histograms are cut")
    print("The histograms are not cut")
    i = h_sig.GetNbinsX()
    # f_out = TFile(output_filename, "UPDATE")
    # f_out.cd("emutau_SS/")
    while(i>1):
        redo_binning = False
        for process in all_process_for_binning:
            h = f_temp.Get("mutautau/"+temp+"/"+process)
            
            if(h.GetBinContent(i) < 5):
                redo_binning = True
        if(redo_binning):
            nbins_old = h.GetNbinsX()
            old_binning = []
            new_binning = arr.array('d')
            for bin in range(1, nbins_old+2):
                old_binning.append(h.GetBinLowEdge(bin))
                if(bin != i):
                    new_binning .append(h.GetBinLowEdge(bin))
            
            for process in all_process:
                h = f_temp.Get("mutautau/"+temp+"/"+process)
                hnew = h.Rebin(len(new_binning)-1, "hnew", new_binning)
                h.Delete()
                hnew.SetTitle(process)
                hnew.SetName(process)
                hnew.Write()

    
        i = i - 1
    
    f_out = TFile(output_filename, "UPDATE")
    f_out.mkdir("mutautau/"+temp+"/")
    f_out.cd("mutautau/"+temp+"/")
    for process in all_process:
            h = f_temp.Get("mutautau/"+temp+"/"+process)
            h.Write()
    f_out.Close()  
f_temp.Close()  

# command_to_copy_rootfile = "cp /afs/cern.ch/user/l/lmollier/git_PDM/analysis/results/rootfile_limit_setting_v2.root /afs/cern.ch/user/l/lmollier/git_PDM/analysis/limit_setting/CMSSW_10_2_13/src/HNL/"
# gSystem.Exec("rm  /afs/cern.ch/user/l/lmollier/git_PDM/analysis/limit_setting/CMSSW_10_2_13/src/HNL/rootfile_limit_setting.root")
# gSystem.Exec(command_to_copy_rootfile)
print("File created. \n>>>DONE")

# #emutau_os
# input_file_emutau_OS =  "histograms_output_220505_optimize_IDwp_run_13_rebin_10.root"
# f_in = TFile(input_file_emutau_OS, "READ")
# f_out = TFile(output_filename, "UPDATE")

# f_out.mkdir("emutau_OS/")
# f_out.cd("emutau_OS/")
# variables_emutau = ["transverse_mass_etau_emutau_OS","comb_mass_emu_emutau_OS","MT_total_emutau_OS"]
# HNL_mass_list = [100,125,150,200,250,300,350,400,450,500,600,700,800,900]
# for var in variables_emutau:
#     name_histo_bck = var+"/output/h_all_bck_"+var
#     h_bck= f_in.Get(name_histo_bck)
#     temp =  var
#     temp = temp.replace("_emutau_OS","")
#     new_name = "bck_" +temp
#     h_bck.SetName(new_name)
#     h_bck.SetTitle(new_name)
#     curr_maxbin = h_bck.FindLastBinAbove()
#     cut = curr_maxbin
#     h_bck.Write()
#     for HNL_mass in HNL_mass_list:
#         HNL_mass = "HNL"+str(HNL_mass)
#         name_histo_sig = var+"/raw histos/"+ HNL_mass + "_"+ var
#         h_sig = f_in.Get(name_histo_sig)

#         temp = var
#         temp = temp.replace("_emutau_OS","")
#         new_name = HNL_mass + "_"+ temp
#         h_sig.SetName(new_name)
#         h_sig.SetTitle(new_name)
#         h_sig = f_in.Get(name_histo_sig)
#         curr_maxbin = h_sig.FindLastBinAbove()
#         if(curr_maxbin<cut):
#             cut = curr_maxbin
#         h_sig.Write()
# f_out.Close()      

# # #etautau
# input_file_etautau =  "histograms_output_220509_optimize_IDwp_etautau_run_9_rebin_10.root"
# f_in = TFile(input_file_etautau, "READ")
# f_out = TFile(output_filename, "UPDATE")
# f_out.mkdir("etautau/")
# f_out.cd("etautau/")
# variables_etautau = ["transverse_mass_tauMET_etautau","transverse_mass_tau2tau_etautau","dr_etau2_etautau","MT_total_etautau"]
# HNL_mass_list = [100,125,150,200,250,300,350,400,450,500,600,700,800,900]
# for var in variables_etautau:
#     name_histo_bck = var+"/output/h_all_bck_"+var
#     h_bck= f_in.Get(name_histo_bck)
#     temp = var
#     temp = temp.replace("_etautau","")
#     new_name = "bck_" + temp
#     h_bck.SetName(new_name)
#     h_bck.SetTitle(new_name)
#     curr_maxbin = h_bck.FindLastBinAbove()
#     cut = curr_maxbin
#     h_bck.Write()
#     for HNL_mass in HNL_mass_list:
#         HNL_mass = "HNL"+str(HNL_mass)
#         name_histo_sig = var+"/raw histos/"+ HNL_mass + "_"+ var
#         h_sig = f_in.Get(name_histo_sig)
#         temp =  var
#         temp = temp.replace("_etautau","")
#         new_name = HNL_mass + "_"+ temp
#         h_sig.SetName(new_name)
#         h_sig.SetTitle(new_name)
#         curr_maxbin = h_sig.FindLastBinAbove()
#         if(curr_maxbin<cut):
#             cut = curr_maxbin
#         h_sig.Write()
# f_out.Close()  

# # #mutautau
# input_file_mutautau =  "histograms_output_220510_optimize_IDwp_mutautau_run_1_rebin_5.root"
# f_in = TFile(input_file_mutautau, "READ")
# f_out = TFile(output_filename, "UPDATE")
# f_out.mkdir("mutautau/")
# f_out.cd("mutautau/")
# variables_mutautau = ["transverse_mass_tauMET_mutautau","comb_mass_tau2tau_mutautau","comb_mass_mutau2_mutautau","dr_mutau2_mutautau"]
# HNL_mass_list = [100,125,150,200,250,300,350,400,450,500,600,700,800,900]
# for var in variables_mutautau:
#     name_histo_bck = var+"/output/h_all_bck_"+var
    
#     h_bck= f_in.Get(name_histo_bck)
#     temp =  var
#     temp = temp.replace("_mutautau","")
#     new_name = "bck_" +temp
#     h_bck.SetName(new_name)
#     h_bck.SetTitle(new_name)
#     curr_maxbin = h_bck.FindLastBinAbove()
#     cut = curr_maxbin
#     h_bck.Write()
#     for HNL_mass in HNL_mass_list:
#         HNL_mass = "HNL"+str(HNL_mass)
#         name_histo_sig = var+"/raw histos/"+ HNL_mass + "_"+ var
        
#         h_sig = f_in.Get(name_histo_sig)
#         temp =  var
#         temp = temp.replace("_mutautau","")
#         new_name = HNL_mass + "_"+temp
#         h_sig.SetName(new_name)
#         h_sig.SetTitle(new_name)
#         curr_maxbin = h_sig.FindLastBinAbove()
#         if(curr_maxbin<cut):
#             cut = curr_maxbin
#         h_sig.Write()
# f_out.Close()  