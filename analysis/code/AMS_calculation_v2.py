from re import S
from ROOT import TFile, TGraph, kBlue, gStyle, kBlack, TCanvas, kPink, gPad, TGaxis , kRed, TH1D,kGreen
import numpy as np
from HNLAnalysis import HNLAnalysis
from HNLAnalysis_etautau import HNLAnalysis_etautau
from HNLAnalysis_mutautau import HNLAnalysis_mutautau
from array import array

path_to_results = "/afs/cern.ch/user/l/lmollier/git_PDM/analysis/results/"

###
#return the list of the cuts with the ams score
def compute_total_ams(h_bck, h_sig):
    n_bins = h_sig.GetNbinsX()
    total_ams = 0 
    i = n_bins+1
    while(i>=0):
        b=0
        s=0
        while(b<5 and i>=0):
            b += h_bck.Integral(i,i)
            s += h_sig.Integral(i,i)
            i-=1
        if(b>=5):
            total_ams += compute_AMS(b,s)**2
        i-=1
    return np.sqrt(total_ams)


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



def best_var_ams_score(e):
    return e[1]



def compute_ams_for_HNLmass(HNL_mass, signs, variables, name_file, output_file, rebining):
    print("AMS analysis for HNL mass : " + HNL_mass + " (" + signs+")")
    var_ams = []
    HNL_scale = {
        "HNL100":1. ,
        "HNL125": 1. ,
        "HNL150":1. ,
        "HNL200":1. ,
        "HNL250":1. ,
        "HNL300":1. ,
        "HNL350": 1. ,
        "HNL400": 1. ,
        "HNL450": 1. ,
        "HNL500": 1. ,
        "HNL600": 1. ,
        "HNL700": 1. ,
        "HNL800": 1. ,
        "HNL900": 1. ,
        }
    if("etautau" in name_file):
         HNL_scale = {
        "HNL100": 1. ,
        "HNL125": 1. ,
        "HNL150": 1. ,
        "HNL200": 1. ,
        "HNL250": 1. ,
        "HNL300": 1. ,
        "HNL350": 1. ,
        "HNL400": 1. ,
        "HNL450": 1. ,
        "HNL500": 1. ,
        "HNL600": 1. ,
        "HNL700": 1. ,
        "HNL800": 1. ,
        "HNL900": 1. ,
        }
    elif("mutautau" in name_file):
        HNL_scale = {
        "HNL100": 1. ,
        "HNL125": 1. ,
        "HNL150": 1. ,
        "HNL200": 1. ,
        "HNL250": 1. ,
        "HNL300": 1. ,
        "HNL350": 1. ,
        "HNL400": 1. ,
        "HNL450": 1. ,
        "HNL500": 1. ,
        "HNL600": 1. ,
        "HNL700": 1. ,
        "HNL800": 1. ,
        "HNL900": 1. ,
        }
    
    # HNL_scale = {
    #     "HNL100": 0.14,
    #     "HNL125": 0.075,
    #     "HNL150": 0.05,
    #     "HNL200": 0.03,
    #     "HNL250": 0.02,
    #     "HNL300": 0.01,
    #     "HNL350": 0.0075,
    #     "HNL400": 0.007,
    #     "HNL450": 0.005,
    #     "HNL500": 0.0035,
    #     "HNL600": 0.0025,
    #     "HNL700": 0.002,
    #     "HNL800": 0.0015,
    #     "HNL900": 0.001,
    #     }
    # if("etautau" in name_file):
    #      HNL_scale = {
    #     "HNL100": 0.075,
    #     "HNL125": 0.025,
    #     "HNL150": 0.02,
    #     "HNL200": 0.0075,
    #     "HNL250": 0.004,
    #     "HNL300": 0.0035,
    #     "HNL350": 0.002,
    #     "HNL400": 0.002,
    #     "HNL450": 0.0015,
    #     "HNL500": 0.0013,
    #     "HNL600": 0.00075,
    #     "HNL700": 0.0006,
    #     "HNL800": 0.0005,
    #     "HNL900": 0.0005,
    #     }
    # elif("mutautau" in name_file):
    #     HNL_scale = {
    #     "HNL100": 0.06,
    #     "HNL125": 0.04,
    #     "HNL150": 0.02,
    #     "HNL200": 0.01,
    #     "HNL250": 0.005,
    #     "HNL300": 0.005,
    #     "HNL350": 0.002,
    #     "HNL400": 0.002,
    #     "HNL450": 0.0015,
    #     "HNL500": 0.0015,
    #     "HNL600": 0.0011,
    #     "HNL700": 0.0008,
    #     "HNL800": 0.0006,
    #     "HNL900": 0.0005,
    #     }
    
    for var in variables:
        
        if(signs in var or "tautau" in name_file):
            #print("\tfor variable "+var)
            
            name_histo_bck = var+"/output/h_all_bck_"+var
            name_histo_sig = var+"/raw histos/"+ HNL_mass + "_"+ var
            # print(name_histo_bck)
            #print(name_histo_sig)
            f = TFile(name_file, "READ")
            h_bck= f.Get(name_histo_bck)
            h_sig = f.Get(name_histo_sig)
            # if("phi" in var):
            #     print("#################")
            #     print("#################")
            #     print("#################")
            #     print(name_histo_sig)
            #     print(h_sig)
            h_sig.Scale(HNL_scale[HNL_mass])
            # if(not ("charge" in var or "phi" in var or "dr" in var or "eta_" in var)):   
            #     h_sig.Rebin(int(rebining))
            #     h_bck.Rebin(int(rebining))
            # # print("\t\thistograms are imported:")
            # print("\t\t"+str(h_bck))
            # print("\t\t"+str(h_sig))
            ams_score = compute_total_ams(h_bck, h_sig)
            var_ams.append((var, ams_score , HNL_mass, signs))
            f_out = TFile(output_file, "UPDATE")
            c = TCanvas(HNL_mass + "_ams_"+var,HNL_mass + "_ams_"+var)
            f_out.mkdir(HNL_mass+"_"+signs+"/")
            f_out.cd(HNL_mass+"_"+signs+"/")
            gStyle.SetLineWidth(2)
            gStyle.SetLineColor(kBlack)
            h_bck.SetFillColor(kBlue)
            h_bck.SetLineColor(kBlue)
            h_sig.SetFillColor(kPink)
            h_sig.SetLineColor(kPink)
            h_bck.SetFillStyle(3003)
            h_sig.SetFillStyle(3003)
            h_bck.Draw("HIST")
            max_sig = 2*h_sig.GetMaximum()
            max_bck = h_bck.GetMaximum()
            scale_sig = max_bck/max_sig
            h_sig.Scale(scale_sig)
            max_x = h_bck.GetBinLowEdge(h_bck.GetNbinsX()+1)
            axis = TGaxis(max_x,0, max_x, max_bck,0,max_sig,10, "+L")
            axis.SetLabelOffset(0.035)
            axis.SetLabelColor(kRed)
            axis.SetLineColor(kRed)
            axis.SetName("counts sig")
            axis.Draw("SAME")
            h_sig.Draw("HIST SAME")

            
            c.BuildLegend(0.65,0.65,0.85,0.85,"")
            c.Write()
            f_out.Close()
    var_ams.sort(reverse=True, key=best_var_ams_score)
    var1_ams = var_ams[0]
    var2_ams = var_ams[1]
    var3_ams = var_ams[2]
    var4_ams = var_ams[3]
    print("\tThe best overall cut for HNL mass "+ var1_ams[2] +"("+signs+") is made with var : " + var1_ams[0]+ " , with ams score = " + str(var1_ams[1]))
    print("\tThe second best overall cut for HNL mass "+ var2_ams[2] +"("+signs+") is made with var : " + var2_ams[0]+ " , with ams score = " + str(var2_ams[1]))
    print("\tThe third best overall cut for HNL mass "+ var3_ams[2] +"("+signs+") is made with var : " + var3_ams[0]+ " , with ams score = " + str(var3_ams[1]))
    print("\tThe fourth best overall cut for HNL mass "+ var4_ams[2] +"("+signs+") is made with var : " + var4_ams[0]+ " , with ams score = " + str(var4_ams[1]))

    return (var1_ams, var2_ams, var3_ams,var4_ams)



def AMS_calculation_v2(tag, rebining):

    if("etautau" in tag):
        var_axis_pair = HNLAnalysis_etautau.get_var_axis_pairs()
    elif("mutautau" in tag):
        var_axis_pair = HNLAnalysis_mutautau.get_var_axis_pairs()
    else:    
        var_axis_pair = HNLAnalysis.get_var_axis_pairs()
    variables = []
    for v in var_axis_pair:
        #if("mass_mu_1" not in v[0] and "mass_e_1" not in v[0] and "mass_tau_1" not in v[0] ):
                variables.append(v[0])

    #tag = '220504_bck_HNL_v1'
    name_file = path_to_results+'histograms_output_'+tag+'.root'
    name_file = path_to_results+'histograms_output_'+tag+"_rebin_"+str(rebining)+'.root'
    print("file analysed : " + name_file)
    output_file = path_to_results+'ams_v2_output_'+tag+"_rebin_"+str(rebining)+'.root'
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
    ]


    list_mass_vs_best_ams = []
    for mass in list_HNL_mass:
        if("tautau" in name_file):
            list_mass_vs_best_ams.append(compute_ams_for_HNLmass(mass, "", variables, name_file,output_file, rebining))

        else:    
            list_mass_vs_best_ams.append(compute_ams_for_HNLmass(mass, "SS", variables, name_file,output_file, rebining))
            list_mass_vs_best_ams.append(compute_ams_for_HNLmass(mass, "OS", variables,name_file,output_file, rebining))

    print("output file is : "+output_file)
    list_mass_vs_best_ams_second  = []
    list_mass_vs_best_ams_third  = []
    list_mass_vs_best_ams_first  = []
    list_mass_vs_best_ams_fourth  =[]
    print(len(list_mass_vs_best_ams))
    for i in range(0, len(list_mass_vs_best_ams)):
        list_mass_vs_best_ams_second.append(list_mass_vs_best_ams[i][1])
        list_mass_vs_best_ams_third.append(list_mass_vs_best_ams[i][2])
        list_mass_vs_best_ams_first.append(list_mass_vs_best_ams[i][0])
        list_mass_vs_best_ams_fourth.append(list_mass_vs_best_ams[i][3])

    print("\n\n\n\n\n")
    print("################")
    print("################")
    print(list_mass_vs_best_ams_first)
    print("FINAL OUTPUT: ")
    for i in range(0, len(list_mass_vs_best_ams)):
        print("########")
        print("The best overall cut for HNL mass "+ list_mass_vs_best_ams_first[i][2] +" (" +list_mass_vs_best_ams_first[i][3]+") is made with var : " + list_mass_vs_best_ams_first[i][0]+ " , with ams score = " + str(list_mass_vs_best_ams_first[i][1]))
        print("The 2nd overall cut for HNL mass "+ list_mass_vs_best_ams_second[i][2] +" (" +list_mass_vs_best_ams_second[i][3]+") is made with var : " + list_mass_vs_best_ams_second[i][0]+ " , with ams score = " + str(list_mass_vs_best_ams_second[i][1]))
        print("The 3rd overall cut for HNL mass "+ list_mass_vs_best_ams_third[i][2] +" (" +list_mass_vs_best_ams_third[i][3]+") is made with var : " + list_mass_vs_best_ams_third[i][0]+ " , with ams score = " + str(list_mass_vs_best_ams_third[i][1]))
        print("The 3rd overall cut for HNL mass "+ list_mass_vs_best_ams_third[i][2] +" (" +list_mass_vs_best_ams_third[i][3]+") is made with var : " + list_mass_vs_best_ams_third[i][0]+ " , with ams score = " + str(list_mass_vs_best_ams_third[i][1]))
        print("The 4th overall cut for HNL mass "+ list_mass_vs_best_ams_fourth[i][2] +" (" +list_mass_vs_best_ams_fourth[i][3]+") is made with var : " + list_mass_vs_best_ams_fourth[i][0]+ " , with ams score = " + str(list_mass_vs_best_ams_fourth[i][1]))

    print("################")
    print("################")


    f_out = TFile(output_file, "UPDATE")
    f_out.mkdir("best_variables/")
    f_out.cd("best_variables/")
    for i in range(0, len(list_mass_vs_best_ams_first)):
        mass = list_mass_vs_best_ams_first[i][2]
        var =  list_mass_vs_best_ams_first[i][0]
        signs = list_mass_vs_best_ams_first[i][3]
        best_hist = f_out.Get(mass +"_"+signs+ "/"+mass+"_ams_"+var)
        name = mass+"_"+signs+"_"+var
        best_hist.SetName(name)
        best_hist.SetTitle(name)
        best_hist.Write()

    f_out.mkdir("2nd_best_variables/")
    f_out.cd("2nd_best_variables/")
    for i in range(0, len(list_mass_vs_best_ams_second)):
        mass = list_mass_vs_best_ams_second[i][2]
        var =  list_mass_vs_best_ams_second[i][0]
        signs = list_mass_vs_best_ams_second[i][3]
        best_hist = f_out.Get(mass +"_"+signs+ "/"+mass+"_ams_"+var)
        name = mass+"_"+signs+"_"+var
        best_hist.SetName(name)
        best_hist.SetTitle(name)
        best_hist.Write()


    f_out.mkdir("3rd_best_variables/")
    f_out.cd("3rd_best_variables/")
    for i in range(0, len(list_mass_vs_best_ams_third)):
        mass = list_mass_vs_best_ams_third[i][2]
        var =  list_mass_vs_best_ams_third[i][0]
        signs = list_mass_vs_best_ams_third[i][3]
        best_hist = f_out.Get(mass +"_"+signs+ "/"+mass+"_ams_"+var)
        name = mass+"_"+signs+"_"+var
        best_hist.SetName(name)
        best_hist.SetTitle(name)
        best_hist.Write()


    f_out.mkdir("4th_best_variables/")
    f_out.cd("4th_best_variables/")
    for i in range(0, len(list_mass_vs_best_ams_fourth)):
        mass = list_mass_vs_best_ams_fourth[i][2]
        var =  list_mass_vs_best_ams_fourth[i][0]
        signs = list_mass_vs_best_ams_fourth[i][3]
        best_hist = f_out.Get(mass +"_"+signs+ "/"+mass+"_ams_"+var)
        name = mass+"_"+signs+"_"+var
        best_hist.SetName(name)
        best_hist.SetTitle(name)
        best_hist.Write()


    print("\n\n\n\n\n")

    print(f'| {"HNL mass":<20}', end='')
    print(f'| {"signs":<6}', end='')
    print(f'| {"variable 1":<32}', end='')
    print(f'| {"AMS score":<15}', end='')
    print(f'| {"variable 2":<32}', end='')
    print(f'| {"AMS score":<15}', end='')
    print(f'| {"variable 3":<32}', end='')
    print(f'| {"AMS score":<15}', end='')
    print(f'| {"variable 4":<32}', end='')
    print(f'| {"AMS score|":<15}', end='\n')

    print(u'\u2500' * 74)
    for i in range(0, len(list_mass_vs_best_ams_first)):
        x = list_mass_vs_best_ams_first[i]
        mass =  x[2]
        var = x[0]
        signs = x[3]
        ams = str(x[1])

        x2 = list_mass_vs_best_ams_second[i]
        var2 = x2[0]
        ams2 = str(x2[1])
        x3 = list_mass_vs_best_ams_third[i]
        var3 = x3[0]
        ams3 = str(x3[1])
        x4= list_mass_vs_best_ams_fourth[i]
        var4 = x4[0]
        ams4 = str(x4[1])

        print(f'| {mass:<20}', end='')
        print(f'| {signs:<6}', end='')
        print(f'| {var:<32} ', end='')
        print(f'| {ams:<40} |', end='')

        print(f'| {var2:<32} ', end='')
        print(f'| {ams2:<40} |', end='')

        print(f'| {var3:<32} ', end='')
        print(f'| {ams3:<40} |', end='')

        print(f'| {var4:<32} ', end='')
        print(f'| {ams4:<40} |', end='\n')
        print(u'\u2500' * 74)
        

    name_txt_file_output  = path_to_results+'ams_v2_output_'+tag+"rebin_"+str(rebining)+'.txt'
    with open(name_txt_file_output, 'w') as f:
        f.write("Here are the best ams to make for each HNL mass. \nThe histograms can be found in the root file called : " + output_file + " \n")
        for i in range(0, len(list_mass_vs_best_ams_first)):
            f.write("The best overall cut for HNL mass "+ list_mass_vs_best_ams_first[i][2] +" ("+list_mass_vs_best_ams_first[i][3]+") is made with var : " + list_mass_vs_best_ams_first[i][0]+ " , with ams score = " + str(list_mass_vs_best_ams_first[i][1])+"\n")
            f.write("The second best overall cut for HNL mass "+ list_mass_vs_best_ams_second[i][2] +" ("+list_mass_vs_best_ams_second[i][3]+") is made with var : " + list_mass_vs_best_ams_second[i][0]+ " , with ams score = " + str(list_mass_vs_best_ams_second[i][1])+"\n")
            f.write("The third best overall cut for HNL mass "+ list_mass_vs_best_ams_third[i][2] +" ("+list_mass_vs_best_ams_third[i][3]+") is made with var : " + list_mass_vs_best_ams_third[i][0]+ " , with ams score = " + str(list_mass_vs_best_ams_third[i][1])+"\n")
            f.write("The fourth best overall cut for HNL mass "+ list_mass_vs_best_ams_fourth[i][2] +" ("+list_mass_vs_best_ams_fourth[i][3]+") is made with var : " + list_mass_vs_best_ams_fourth[i][0]+ " , with ams score = " + str(list_mass_vs_best_ams_fourth[i][1])+"\n")
            f.write("###")

        f.write("\n\n\n\n")
        f.write("HNL mass \t")
        f.write("variable \t\t\t")
        f.write("ams score")
        f.write("variable 2\t\t\t")
        f.write("ams score 2")
        f.write("variable 3 \t\t\t")
        f.write("ams score 3")
        f.write("variable 4 \t\t\t")
        f.write("ams score 4\n")
        for i in range(0, len(list_mass_vs_best_ams_first)):
            f.write("\n")
            f.write(list_mass_vs_best_ams_first[i][2]+" ("+list_mass_vs_best_ams_first[i][3]+")" + "\t " + list_mass_vs_best_ams_first[i][0]+ "\t\t\t"+"\t"+  str(list_mass_vs_best_ams_first[i][1])+ "\t " + list_mass_vs_best_ams_second[i][0]+ "\t\t\t"+"\t"+  str(list_mass_vs_best_ams_second[i][1])+ "\t " + list_mass_vs_best_ams_third[i][0]+ "\t\t\t"+"\t"+  str(list_mass_vs_best_ams_third[i][1])+ "\t " + list_mass_vs_best_ams_fourth[i][0]+ "\t\t\t"+"\t"+  str(list_mass_vs_best_ams_fourth[i][1]))

        
        
    name_txt_file_output_for_run  = path_to_results+'best_ams_v2_'+tag+"rebin_"+str(rebining)+'.txt'
    with open(name_txt_file_output_for_run, 'w') as f:
        f.write("HNL mass \t")
        f.write("ams score\n")
        for i in range(0, len(list_mass_vs_best_ams_first)):
            mass = list_mass_vs_best_ams_first[i][2] +"("+list_mass_vs_best_ams_first[i][3]+")\t"
            f.write(mass)
        f.write("\n")
        for i in range(0, len(list_mass_vs_best_ams_first)):
            ams =str(list_mass_vs_best_ams_first[i][1])+"\t"
            f.write(ams)       

    name_txt_file_output_second_for_run  = path_to_results+'second_best_ams_v2_'+tag+"rebin_"+str(rebining)+'.txt'
    with open(name_txt_file_output_second_for_run, 'w') as f:
        f.write("HNL mass \t")
        f.write("ams score\n")
        for i in range(0, len(list_mass_vs_best_ams_second)):
            mass = list_mass_vs_best_ams_second[i][2] +"("+list_mass_vs_best_ams_second[i][3]+")\t"
            f.write(mass)
        f.write("\n")
        for i in range(0, len(list_mass_vs_best_ams_second)):
            ams =str(list_mass_vs_best_ams_second[i][1])+"\t"
            f.write(ams)

    name_txt_file_output_third_for_run  = path_to_results+'third_best_ams_v2_'+tag+"rebin_"+str(rebining)+'.txt'
    with open(name_txt_file_output_third_for_run, 'w') as f:
        f.write("HNL mass \t")
        f.write("ams score\n")
        for i in range(0, len(list_mass_vs_best_ams_third)):
            mass = list_mass_vs_best_ams_third[i][2] +"("+list_mass_vs_best_ams_third[i][3]+")\t"
            f.write(mass)
        f.write("\n")
        for i in range(0, len(list_mass_vs_best_ams_third)):
            ams =str(list_mass_vs_best_ams_third[i][1])+"\t"
            f.write(ams)


    name_txt_file_output_fourth_for_run  = path_to_results+'fourth_best_ams_v2_'+tag+"rebin_"+str(rebining)+'.txt'
    with open(name_txt_file_output_fourth_for_run, 'w') as f:
        f.write("HNL mass \t")
        f.write("ams score\n")
        for i in range(0, len(list_mass_vs_best_ams_fourth)):
            mass = list_mass_vs_best_ams_fourth[i][2] +"("+list_mass_vs_best_ams_fourth[i][3]+")\t"
            f.write(mass)
        f.write("\n")
        for i in range(0, len(list_mass_vs_best_ams_fourth)):
            ams =str(list_mass_vs_best_ams_fourth[i][1])+"\t"
            f.write(ams)            