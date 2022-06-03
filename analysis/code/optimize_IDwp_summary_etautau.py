
from ROOT import  TH2D
import numpy as np
from HNLAnalysis import HNLAnalysis

from array import array


tag = "best_ams_v2_220509_optimize_IDwp_etautau"
number_runs = [1,2,3,4,5,6,7,8,9,10,11,12,13]
parameters_runs = {"1":"8,2,2,90",
                    "2":"8,4,4,90",
                    "3":"8,2,4,90",
                    "4":"8,4,2,90",
                    "5":"8,8,2,90",
                    "6":"16,2,2,90",
                    "7":"16,4,2,90",
                    "8":"32,4,2,90",
                    "9":"32,4,2,80",
                    "10":"64,4,2,90",
                    "11":"64,8,2,90",
                    "12":"64,8,4,90",
                    "13":"64,8,2,80",
                    }

nb_runs = len(number_runs)
runs_file  = []
path = "/afs/cern.ch/user/l/lmollier/git_PDM/analysis/results/"
for i in number_runs:
    runs_file.append(path+tag+"_run_"+str(i)+"rebin_5.txt")
HNL_mass = [100,125,150,200,250,300,350,400,450,500,600,700,800,900]
nb_HNL_mass = len(HNL_mass)

h_ams = TH2D("h_ams","h_ams",nb_runs,1,nb_runs+1,nb_HNL_mass,1,nb_HNL_mass+1)

#print(runs_file)
first = True
output_name = tag+"_summary.txt"
for run_name in runs_file:
    curr_run = run_name.split("_")[-2]
    print(curr_run)
    curr_run = curr_run.split("r")[0]
    print(curr_run)
    with open(run_name, 'r') as f:
        lines = f.readlines()
        if(first):
            first = False
            with open(output_name, 'w') as f_out:
                f_out.write("\t"+lines[1])
                f_out.write("run "+str(curr_run)+"\t"+lines[2]+ "\n")

        else:
            with open(output_name, 'a') as f_out:
                f_out.write("run "+str(curr_run)+"\t"+lines[2]+ "\n")
        ams_score = lines[2].split("\t")[:-1]
        for curr_mass in range(0, len(ams_score)):
            h_ams.Fill(int(curr_run),curr_mass+1, float(ams_score[curr_mass]))
for i in range(1, len(HNL_mass)+1):
    h_ams.GetYaxis().SetBinLabel(i,"HNL "+str(HNL_mass[i-1]))
for i in range(1, len(number_runs)+1):
    r = str(number_runs[i-1])
    h_ams.GetXaxis().SetBinLabel(i,"run "+ r + "("+parameters_runs[r]+")")
#h_ams.GetXaxis().SetTitle("run number (tVSj,tVSe,tVSm,eWP)")
h_ams.SetTitle("HNL mass vs. runs (tVSj, tVSe, tVSm, eWP)")
h_ams.SetStats(0)
h_ams.Draw("TEXT COL")
input("press enter")
            

