
from ROOT import TFile, TCanvas, kRed, kBlue, kBlack, gStyle,TPaveText

print("started")

input_file =  "rootfile_limit_setting.root"
f_in = TFile(input_file, "READ")

states = ["emutau_SS", "emutau_OS", "etautau", "mutautau"]
output_file = "plots_relevant_var_all_final_states.root"
f_out = TFile(output_file,"RECREATE")

for state in states:
    f_out.mkdir(state+"/")
    variables = []
    if (state == "emutau_SS"):
        variables = ["MT_total","transverse_mass_etau","transverse_mass_mutau", "pt_sum_emutau"]
    elif (state == "emutau_OS"):
        variables =  ["pt_sum_mutau","comb_mass_emu","MT_total","pt_tau_1"]
    elif (state == "etautau"):
        variables = ["pt_sum_tau2tau","transverse_mass_tau2tau","dphi_tau2tau","MT_total"]
    elif (state == "mutautau"):
        variables =["MT_total","transverse_mass_tau2tau","dr_mutau2","pt_sum_mutau2tau"]
    for var in variables:
        f_out.cd(state+"/")
        print(var)
        print(state)
        hbck = f_in.Get(state + "/"+var+"/data_obs")
        h100 = f_in.Get(state + "/"+var+"/HNL100")
        h400 = f_in.Get(state + "/"+var+"/HNL400")
        h900 = f_in.Get(state + "/"+var+"/HNL900")
        name = var
        c = TCanvas(name, name)
        gStyle.SetOptTitle(0)
        gStyle.SetOptStat(0)
        hbck.Scale(100/hbck.GetMaximum())
        h100.Scale(100/h100.GetMaximum())
        h900.Scale(100/h900.GetMaximum())
        h400.Scale(100/h400.GetMaximum())
        h100.SetLineColor(kRed)
        h400.SetLineColor(kBlue)
        h900.SetLineColor(kBlack)
        h100.SetLineWidth(2)
        h400.SetLineWidth(2)
        h900.SetLineWidth(2)
        hbck.SetLineColor(kBlue)
        hbck.SetFillColor(kBlue)
        hbck.SetFillStyle(3003)
        
        hbck.Draw("HIST ")
        h100.Draw("HIST SAME ")
        h400.Draw("HIST SAME ")
        h900.Draw("HIST SAME ")
        
        c.BuildLegend(0.65,0.65,0.8,0.8,"")
        c.SetName(var)
        c.SetTitle(var)
        t = TPaveText(0.4, 0.9, 0.6, 1.0, "brNDC")
        t.AddText(var)
        t.Draw("SAME")
        c.Write()
