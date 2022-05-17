from re import S
from ROOT import TFile, TGraph, kBlue, gStyle, kBlack, TCanvas, kPink, gPad, TGaxis , kRed, TH1D,kGreen
import numpy as np
from HNLAnalysis import HNLAnalysis
from HNLAnalysis_etautau import HNLAnalysis_etautau
from HNLAnalysis_mutautau import HNLAnalysis_mutautau

from array import array


###
#return the list of the cuts with the ams score
def compute_cuts(h_bck, h_sig):
    #bins 0 = underflow 
    len_bins = h_bck.GetNbinsX()
    list_ams = [] #list of the value of cut and corresponding ams score 
    total_sig_event = h_sig.Integral(0,len_bins+1)
    for i in range(1, len_bins+1):
        min_bin = i
        max_bin = len_bins+1
        #print("*****")
        #print("min bin = " + str(min_bin)+ " : "+ str(h_bck.GetBinLowEdge(i)))
        b = h_bck.Integral(min_bin, max_bin) #Integral from bin number min_bin to bin number max_bin (both included)
        s = h_sig.Integral(min_bin, max_bin)
        curr_ams = compute_AMS(b,s)
        if(curr_ams != 0 and s > 0*total_sig_event):  #requiring at least 10% of the events after cut
            list_ams.append(((h_bck.GetBinLowEdge(min_bin), h_bck.GetBinLowEdge(max_bin)),curr_ams))
        #print("in range : [ " + str(h_bck.GetBinLowEdge(min_bin)) + " , "+str(h_bck.GetBinLowEdge(max_bin) + h_bck.GetBinWidth(max_bin)) + " ]: b = "+ str(b) + ", s = " + str(s)+", ams = "+ str(curr_ams) )
    return list_ams    

###
#return the list of the cuts with the ams score for angular variables
def compute_cuts_angular(h_bck, h_sig):
    #bins 0 = underflow 
    len_bins = h_bck.GetNbinsX()
    list_ams = [] #list of the value of cut and corresponding ams score 
    total_sig_event = h_sig.Integral(0,len_bins+1)
    for i in range(1, len_bins+1):
        min_bin = i
        for j in range(i, len_bins+1):
            max_bin = j
            b = h_bck.Integral(min_bin, max_bin) #Integral from bin number min_bin to bin number max_bin (both included)
            s = h_sig.Integral(min_bin, max_bin)
            curr_ams = compute_AMS(b,s)
            if(curr_ams != 0 and s > 0.1*total_sig_event):
                list_ams.append(((h_bck.GetBinLowEdge(min_bin), h_bck.GetBinLowEdge(max_bin)),curr_ams))
    return list_ams


def compute_AMS(b, s):
    #print("b = " + str(b))
    #print("s = " + str(s))
    if(b == 0 or b< 0 or s <=0): #
        ams = 0
    else:
        temp = ((s+b)*np.log(1+(s/b))) - s
        temp1 = 2*temp
        ams = np.sqrt(temp1)
    #print("ams = " + str(ams))
    return ams


def ams_score(e):
    return e[1]

def best_var_ams_score(e):
    return e[1][1]



def compute_ams_for_HNLmass(HNL_mass, signs, variables, name_file, output_file):
    print("AMS analysis for HNL mass : " + HNL_mass + " (" + signs+")")
    overall_best_cut = []
    for var in variables:
        
        if(signs in var or "tautau" in name_file):
            # print("\tfor variable "+var)
            name_histo_bck = var+"/output/h_all_bck_"+var
            name_histo_sig = var+"/raw histos/"+ HNL_mass + "_"+ var
            # print(name_histo_bck)
            # print(name_histo_sig)
            f = TFile(name_file, "READ")
            h_bck= f.Get(name_histo_bck)
            h_sig = f.Get(name_histo_sig)

            # print("\t\thistograms are imported:")
            # print("\t\t"+str(h_bck))
            # print("\t\t"+str(h_sig))
            if("eta_" in name_histo_bck or "phi_" in name_histo_bck or "dr_" in name_histo_bck):
                list_cut = compute_cuts_angular(h_bck, h_sig)
            else:
                list_cut = compute_cuts(h_bck, h_sig)
            list_cut_natural_order = list_cut.copy()
            list_cut.sort(reverse=True, key=ams_score)
            best_cut = list_cut[0]
            overall_best_cut.append((var,best_cut, HNL_mass, signs))
            # print("\t\tFor variable : "+ var+ " the best ams score is : " + str(best_cut[1]) + " , with cut at :" + str(best_cut[0]))
            f_out = TFile(output_file, "UPDATE")
            c = TCanvas(HNL_mass + "_ams_"+var,HNL_mass + "_ams_"+var)
            n = len(list_cut_natural_order)
            x = array( 'd' )
            y = array( 'd' )
            for i in range(0, len(list_cut_natural_order)):
                x.append(list_cut_natural_order[i][0][0])
                y.append(list_cut_natural_order[i][1])
            g = TGraph(n, x,y)
            graph_max = max(y)
            f_out.mkdir(HNL_mass+"_"+signs+"/")
            f_out.cd(HNL_mass+"_"+signs+"/")
            g.SetName(HNL_mass + "_ams_"+var)
            g.SetTitle(HNL_mass + "_ams_"+var)
            gStyle.SetLineWidth(2)
            gStyle.SetLineColor(kBlack)


            h_bck.SetFillColor(kBlue)
            h_bck.SetLineColor(kBlue)
            h_sig.SetFillColor(kPink)
            h_sig.SetLineColor(kPink)

            h_bck.SetFillStyle(3003)
            h_sig.SetFillStyle(3003)
            rightmax = 1.1*h_bck.GetMaximum()
            scale = graph_max/rightmax

            rightmax_sig = 2*h_sig.GetMaximum()
            scale_sig = graph_max/rightmax_sig


            h_bck.Scale(scale)
            h_sig.Scale(scale_sig)
            
            max_x = max(h_bck.GetBinLowEdge(h_bck.GetNbinsX()+1), h_sig.GetBinLowEdge(h_sig.GetNbinsX()+1))
            min_x = min(h_bck.GetBinLowEdge(1), h_sig.GetBinLowEdge(1))
            min_y = 0
            max_y = graph_max
            if("eta_" in name_histo_bck or "phi_" in name_histo_bck or "dr_" in name_histo_bck):
                g.Draw("A*")
            else:
                g.Draw("AC*")
            g.GetXaxis().SetLimits(int(min_x),int(max_x))
            g.GetYaxis().SetRangeUser(0,max_y*1.1)
            g.GetYaxis().SetName("AMS score")
            g.GetYaxis().SetTitle("AMS score")
            h_bck.Draw("HIST SAME")
            h_sig.Draw("HIST SAME")
            #draw an axis on the right side
            axis = TGaxis(max_x,0, max_x, max_y,0,rightmax,10, "+L")
            axis.SetLabelOffset(0.035)
            axis.SetLabelColor(kBlue)
            axis.SetLineColor(kBlue)
            axis.SetName("counts bck")
            axis_sig = TGaxis(max_x,0, max_x, max_y*1.05,0,rightmax_sig,10, "+L")
            axis_sig.SetLineColor(kPink)
            axis_sig.SetLabelColor(kPink)
            axis_sig.SetName("counts sig")


            #g.Draw("C* SAME")
            axis.Draw()
            axis_sig.Draw()
            c.BuildLegend(0.65,0.65,0.85,0.85,"")
            c.Write()
            f_out.Close()
    overall_best_cut.sort(reverse=True, key=best_var_ams_score)
    overall_best_cut = overall_best_cut[0]

    print("\tThe best overall cut for HNL mass "+ overall_best_cut[2] +"("+signs+") is made with var : " + overall_best_cut[0]+ " , with ams score = " + str(overall_best_cut[1][1])+". The cut is made as events between "+ str(overall_best_cut[1][0][0]) + " and " + str(overall_best_cut[1][0][1]))

    return overall_best_cut




# variables_1 = ["pt", "eta", "phi", "charge"]
# variables_2 = ["_e", "_tau", "_mu"]
# variables_3 = ["_1"]
# variables_4 = ["_emutau"]
# variables_5 = ["_SS", "_OS"]
# variables =[]
# for i in range(0, len(variables_1)):
#     for j in range(0, len(variables_2)):
#         for k in range(0, len(variables_3)):
#             for l in range(0, len(variables_4)):
#                 for m in range(0, len(variables_5)):
#                     variables.append(variables_1[i]+variables_2[j]+variables_3[k]+variables_4[l]+variables_5[m])

def AMS_calculation(tag, rebining):
    if("etautau" in tag):
        var_axis_pair = HNLAnalysis_etautau.get_var_axis_pairs()
    elif("mutautau" in tag):
        var_axis_pair = HNLAnalysis_mutautau.get_var_axis_pairs()
    else:    
        var_axis_pair = HNLAnalysis.get_var_axis_pairs()
    variables = []
    for v in var_axis_pair:
        if("mass_mu_1" not in v[0] and "mass_e_1" not in v[0] and "mass_tau_1" not in v[0] ):
                variables.append(v[0])

    #tag = '220504_bck_HNL_v1'
    name_file = 'histograms_output_'+tag+"_rebin_"+str(rebining)+'.root'
    print("file analysed : " + name_file)
    output_file = 'ams_output_'+tag+"_rebin_"+str(rebining)+'.root'
    f_out = TFile(output_file, "RECREATE")
    f_out.Close()

    list_HNL_mass = [
    "HNL100",
    "HNL125",
    "HNL150",
    "HNL200",
    "HNL250",
    "HNL300",
    "HNL350",
    "HNL400",
    "HNL450",
    "HNL500",
    "HNL600",
    "HNL700",
    "HNL800",
    "HNL900",




    #"HNL1000"
    ]

    list_mass_vs_best_ams = []
    for mass in list_HNL_mass:
        if("tautau" in name_file):
            list_mass_vs_best_ams.append(compute_ams_for_HNLmass(mass, "", variables, name_file,output_file))

        else:    
            list_mass_vs_best_ams.append(compute_ams_for_HNLmass(mass, "SS", variables, name_file,output_file))
            list_mass_vs_best_ams.append(compute_ams_for_HNLmass(mass, "OS", variables,name_file,output_file))

    print("output file is : "+output_file)

    print("\n\n\n\n\n")
    print("################")
    print("################")
    print("FINAL OUTPUT: ")
    print("requiring at least 10 percent of signal events after cut")
    for i in range(0, len(list_mass_vs_best_ams)):
        print("The best overall cut for HNL mass "+ list_mass_vs_best_ams[i][2] +" (" +list_mass_vs_best_ams[i][3]+") is made with var : " + list_mass_vs_best_ams[i][0]+ " , with ams score = " + str(list_mass_vs_best_ams[i][1][1])+". The cut is made as events between "+ str(list_mass_vs_best_ams[i][1][0][0]) + " and " + str(list_mass_vs_best_ams[i][1][0][1]))

    print("################")
    print("################")


    f_out = TFile(output_file, "UPDATE")
    f_out.mkdir("best_variables/")
    f_out.cd("best_variables/")
    for i in range(0, len(list_mass_vs_best_ams)):
        mass = list_mass_vs_best_ams[i][2]
        var =  list_mass_vs_best_ams[i][0]
        signs = list_mass_vs_best_ams[i][3]
        best_hist = f_out.Get(mass +"_"+signs+ "/"+mass+"_ams_"+var)
        name = mass+"_"+signs+"_"+var
        best_hist.SetName(name)
        best_hist.SetTitle(name)
        best_hist.Write()

    f_out.mkdir("best_cut/")
    f_out.cd("best_cut/")
    for i in range(0, len(list_mass_vs_best_ams)):
        mass = list_mass_vs_best_ams[i][2]
        var =  list_mass_vs_best_ams[i][0]
        signs = list_mass_vs_best_ams[i][3]
        min_cut = list_mass_vs_best_ams[i][1][0][0]
        max_cut = list_mass_vs_best_ams[i][1][0][1]
        best_ams_score = list_mass_vs_best_ams[i][1][1]
        best_hist = f_out.Get(mass +"_"+signs+ "/"+mass+"_ams_"+var)
        name = mass+"_"+signs+"_"+var

        h_cut = TH1D("h_cut", "h_cut", 1, min_cut, max_cut)
        for i in range(0, int(2*best_ams_score)):
            h_cut.Fill(min_cut)
        h_cut.SetFillColor(kGreen-9)
        h_cut.SetFillStyle(3351)
        gStyle.SetHatchesSpacing(8)
        gStyle.SetMarkerStyle(1)
        h_cut.Draw("SAME ][")
        best_hist.SetName(name)
        best_hist.SetTitle(name)
        best_hist.Update()
        best_hist.Write()







    print("\n\n\n\n\n")

    print(f'| {"HNL mass":<20}', end='')
    print(f'| {"signs":<6}', end='')
    print(f'| {"variable   ":<32}', end='')
    print(f'| {"cut             |":<25}', end='')
    print(f'| {"AMS score|":<15}', end='\n')

    print(u'\u2500' * 74)
    for x in list_mass_vs_best_ams:
        mass =  x[2]
        var = x[0]
        signs = x[3]
        ams = str(x[1][1])
        cut = str(x[1][0])
        print(f'| {mass:<20}', end='')
        print(f'| {signs:<6}', end='')
        print(f'| {var:<32} ', end='')
        print(f'| {cut:<15} |', end='')
        print(f'| {ams:<15} |', end='\n')

        print(u'\u2500' * 74)
        

    name_txt_file_output  = 'ams_output_'+tag+"_rebin_"+str(rebining)+'.txt'
    with open(name_txt_file_output, 'w') as f:
        f.write("Here are the best cut to make for each HNL mass. The selection are meant to be between than what is shown in cut.\nRequiring at least 10 percent of signal event after cut. \nThe histograms can be found in the root file called : " + output_file + " \n")
        for i in range(0, len(list_mass_vs_best_ams)):
            f.write("The best overall cut for HNL mass "+ list_mass_vs_best_ams[i][2] +" ("+list_mass_vs_best_ams[i][3]+") is made with var : " + list_mass_vs_best_ams[i][0]+ " , with ams score = " + str(list_mass_vs_best_ams[i][1][1])+". The cut is made as events > "+ str(list_mass_vs_best_ams[i][1][0]) + "\n")
        

        f.write("\n\n\n\n")
        f.write("HNL mass \t")
        f.write("variable \t\t\t")
        f.write("cut\t")
        f.write("ams score\n")

        for i in range(0, len(list_mass_vs_best_ams)):
            f.write("\n")
            f.write(list_mass_vs_best_ams[i][2]+" ("+list_mass_vs_best_ams[i][3]+")" + "\t " + list_mass_vs_best_ams[i][0]+ "\t\t\t"+ str(list_mass_vs_best_ams[i][1][0])+"\t"+  str(list_mass_vs_best_ams[i][1][1]))

        
        
    name_txt_file_output_for_run  = 'best_ams'+tag+"_rebin_"+str(rebining)+'.txt'
    with open(name_txt_file_output_for_run, 'w') as f:
        f.write("HNL mass \t")
        f.write("ams score\n")
        for i in range(0, len(list_mass_vs_best_ams)):
            mass = list_mass_vs_best_ams[i][2] +"("+list_mass_vs_best_ams[i][3]+")\t"
            f.write(mass)
        f.write("\n")
        for i in range(0, len(list_mass_vs_best_ams)):
            ams =str(list_mass_vs_best_ams[i][1][1])+"\t"
            f.write(ams)       
