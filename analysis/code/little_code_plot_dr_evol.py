
from ROOT import TFile, TCanvas, kRed, kBlue, kBlack

print("started")


input_file =  "histograms_output_220510_optimize_IDwp_mutautau_run_1_rebin_10.root"
f_in = TFile(input_file, "READ")

h100 = f_in.Get("dr_mutau2_mutautau/raw histos/HNL100_dr_mutau2_mutautau")
h125 = f_in.Get("dr_mutau2_mutautau/raw histos/HNL125_dr_mutau2_mutautau")
h150 = f_in.Get("dr_mutau2_mutautau/raw histos/HNL150_dr_mutau2_mutautau")

hbck = f_in.Get("dr_mutau2_mutautau/output/h_all_bck_dr_mutau2_mutautau")


hbck.Scale(1)
c = TCanvas("c","c")
h100.SetLineColor(kRed)
h125.SetLineColor(kBlue)
h150.SetLineColor(kBlack)
h100.SetLineWidth(2)
h125.SetLineWidth(2)
h150.SetLineWidth(2)
hbck.SetLineColor(kBlue)
hbck.SetFillColor(kBlue)
hbck.SetFillStyle(3003)
hbck.Draw(" HIST")

h100.Draw("HIST SAME")
h125.Draw("HIST SAME")
h150.Draw("HIST SAME")
c.BuildLegend(0.7,0.7,0.9,0.9,"")
input("press enter")