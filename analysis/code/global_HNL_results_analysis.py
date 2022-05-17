#from AMS_calculation import AMS_calculation
from AMS_calculation_v2 import AMS_calculation_v2
from AMS_calculation_v3 import AMS_calculation_v3
from HNL_plot_analysis import HNL_plot_analysis
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('tag', help='tag ')
parser.add_argument('which', help='which analysis to make, either 1 for only HNL_plot_analysis, 2 for only AMS_calculation and 12 for both ')
parser.add_argument('rebining', help='rebin the histograms', default = 1)

args = parser.parse_args()

tag = args.tag
rebining = args.rebining
which = args.which
print("the rebining is " + str(rebining))
if(which == "12"):
    print("HNL_plot_analysis & AMS_calculation will be performed")
if(which == "1"):
    print("Only HNL_plot_analysis will be performed")
if(which == "2"):
    print("Only AMS_calculation will be performed")

print("tag = " + tag)
if("1" in which):
    print("HNL_plot_analysis is beginning")
    HNL_plot_analysis(tag, int(rebining))
if("2" in which):
    print("AMS_calculation_v2 is beginning")
    AMS_calculation_v2(tag, int(rebining))