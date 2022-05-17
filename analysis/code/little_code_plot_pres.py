
from ROOT import TFile, TCanvas, kRed, kBlue, kBlack

print("started")


input_file =  "histograms_output_220505_optimize_IDwp_run_13_rebin_10.root"
f_in = TFile(input_file, "READ")

h100 = f_in.Get("pt_tau_1_emutau_SS/raw histos/HNL100_pt_tau_1_emutau_SS")
h400 = f_in.Get("pt_tau_1_emutau_SS/raw histos/HNL400_pt_tau_1_emutau_SS")
h900 = f_in.Get("pt_tau_1_emutau_SS/raw histos/HNL900_pt_tau_1_emutau_SS")

hbck = f_in.Get("pt_tau_1_emutau_SS/output/h_all_bck_pt_tau_1_emutau_SS")


hbck.Scale(8)
c = TCanvas("c","c")
h100.SetLineColor(kRed)
h400.SetLineColor(kBlue)
h900.SetLineColor(kBlack)
h100.SetLineWidth(2)
h400.SetLineWidth(2)
h900.SetLineWidth(2)
hbck.SetLineColor(kBlue)
hbck.SetFillColor(kBlue)
hbck.SetFillStyle(3003)
h100.SetTitle("HNL 100 GeV" )
h400.SetTitle("HNL 400 GeV")

h900.SetTitle("HNL 900 GeV")
hbck.SetTitle("Backgrounds")
hbck.GetXaxis().SetTitle("#{p_T}^#tau")
hbck.Draw(" HIST")

h100.Draw("HIST SAME")
h400.Draw("HIST SAME")
h900.Draw("HIST SAME")
c.BuildLegend(0.7,0.7,0.9,0.9,"")
input("press enter")