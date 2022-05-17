#!/usr/bin/env python
import os

import pickle
import argparse

from coffea import processor
from coffea.nanoevents import NanoAODSchema
NanoAODSchema.warn_missing_crossrefs = True

from HNLAnalysis_etautau import HNLAnalysis_etautau
from CountEvents import CountEvents

from samples import signal_samples
from helpers import files_from_dir, files_from_dirs, one_file_from_dir, HNL_from_dir

parser = argparse.ArgumentParser()
parser.add_argument('tag', help='tag that will be added to produced pkl files')
parser.add_argument('--test', action='store_true', default=False, help='test with a subset of the files')
#parser.add_argument('-d', '--local_dir', default='/media/sf_PDM/data/HeavyNeutrino_trilepton/', help='folder which contains subfolder with custom NanoAOD files')
parser.add_argument('-d', '--local_dir', default='/eos/user/p/pdebryas/HNL/skimmed_samples/nanoAOD/2018/HNL_tau/')
args = parser.parse_args()

local_dir = args.local_dir
bck_dir = '/eos/user/l/lmollier/PDM/data/bck_v3/MC_background/'
print("backgrounds taken from : " + bck_dir)
 #Drell-Yann
local_dir_bck_DYJets_To_LL_M_50_madgraphMLM = bck_dir+'DYJets_To_LL_M-50-madgraphMLM'
local_dir_bck_DYJets_To_LL_M_50_amcatnloFXFX = bck_dir+'DYJets_To_LL_M-50-amcatnloFXFX'
local_dir_bck_DY1Jets_To_LL_M_50 = bck_dir+'DY1Jets_To_LL_M-50'
local_dir_bck_DY2Jets_To_LL_M_50 = bck_dir+'DY2Jets_To_LL_M-50'
local_dir_bck_DY3Jets_To_LL_M_50 = bck_dir+'DY3Jets_To_LL_M-50'
local_dir_bck_DY4Jets_To_LL_M_50 = bck_dir+'DY4Jets_To_LL_M-50'
local_dir_bck_DYJets_To_LL_M_10to50 = bck_dir+'DYJets_To_LL_M-10to50'
#ELECRTOWEAK
local_dir_bck_EWK_WMinus2Jets_W_To_LNu_M_50 = bck_dir+'EWK_WMinus2Jets_W_To_LNu_M-50'
local_dir_bck_EWK_WPlus2Jets_W_To_LNu_M_50 = bck_dir+'EWK_WPlus2Jets_W_To_LNu_M-50'
local_dir_bck_EWK_Z2Jets_Z_To_LL_M_50 = bck_dir+'EWK_Z2Jets_Z_To_LL_M-50'
#TTbar
local_dir_bck_TT_To_2L2Nu = bck_dir+'TT_To_2L2Nu'
local_dir_bck_TT_To_SemiLeptonic = bck_dir+'TT_To_SemiLeptonic'
local_dir_bck_TT_To_Hadronic = bck_dir+'TT_To_Hadronic'

#singletop
local_dir_bck_ST_tW_antitop_5f_inclusiveDecays = bck_dir+'ST_tW_antitop_5f_inclusiveDecays'
local_dir_bck_ST_tW_top_5f_inclusiveDecays = bck_dir+'ST_tW_top_5f_inclusiveDecays'
local_dir_bck_ST_t_channel_antitop_4f_InclusiveDecays = bck_dir+'ST_t-channel_antitop_4f_InclusiveDecays'
local_dir_bck_ST_t_channel_top_4f_InclusiveDecays = bck_dir+'ST_t-channel_top_4f_InclusiveDecays'
#W+jets
local_dir_bck_WJets_To_LNu = bck_dir+'WJets_To_LNu'
local_dir_bck_WJ1ets_To_LNu = bck_dir+'WJ1ets_To_LNu'
local_dir_bck_W2Jets_To_LNu = bck_dir+'W2Jets_To_LNu'
local_dir_bck_W3Jets_To_LNu = bck_dir+'W3Jets_To_LNu'
local_dir_bck_W4Jets_To_LNu = bck_dir+'W4Jets_To_LNu'
#DiBoson
local_dir_bck_ZZ = bck_dir+'ZZ'
local_dir_bck_ZZ_To_4L = bck_dir+'ZZ_To_4L'
local_dir_bck_ZZ_To_2L2Nu = bck_dir+'ZZ_To_2L2Nu'
local_dir_bck_WW = bck_dir+'WW'
local_dir_bck_WW_To_2L2Nu = bck_dir+'WW_To_2L2Nu'
local_dir_bck_WZ = bck_dir+'WZ'
local_dir_bck_WZ_To_3LNu = bck_dir+'WZ_To_3LNu'
#TT + bosons
local_dir_bck_TTWJets_To_LNu = bck_dir+'TTWJets_To_LNu'
local_dir_bck_TTZ_To_LLNuNu_M_10 = bck_dir+'TTZ_To_LLNuNu_M-10'
local_dir_bck_TTWW = bck_dir+'TTWW'
local_dir_bck_TTZZ = bck_dir+'TTZZ'
local_dir_bck_TTWZ = bck_dir+'TTWZ'
#Tribosons
local_dir_bck_WWW_4F= bck_dir+'WWW_4F'
local_dir_bck_ZZZ= bck_dir+'ZZZ'
local_dir_bck_WZZ= bck_dir+'WZZ'
local_dir_bck_WWZ_4F= bck_dir+'WWZ_4F'


local_dir_bck_ZZ_old = bck_dir+'ZZ_old'


###
#OLD ONES
###
# local_dir_bck_zz = '/media/sf_PDM/data/backgrounds/zz'
# local_dir_bck_wz = '/media/sf_PDM/data/backgrounds/wz'
#local_dir_bck_tt2l2nu = '/media/sf_PDM/data/backgrounds/tt2l2nu'
# local_dir_bck_tt2l2nu = bck_dir+'tt2l2nu'
# local_dir_bck_ZZ_To_4L = bck_dir+'ZZ_To_4L'
# local_dir_bck_ZZ = bck_dir+'ZZ'
# local_dir_bck_WZ = bck_dir+'WZ_To_3LNu'
# local_dir_bck_WW = bck_dir+'WW_To_2L2Nu'
# local_dir_bck_TTZ = bck_dir+'TTZ_To_LLNuNu'
# local_dir_bck_TTWJets = bck_dir+'TTWJets_To_LNu'
# local_dir_bck_DYJets_To_LL_M_50 = bck_dir+'DYJets_To_LL_M-50'
# local_dir_bck_DYJets_To_LL_M_10to50 = bck_dir+'DYJets_To_LL_M-10to50'
# local_dir_bck_WZ_old = bck_dir+'WZ_old'
# local_dir_bck_ZZ_old = bck_dir+'ZZ_old'



test = args.test
tag = args.tag


samples = {
#HNL

# "HNL85": HNL_from_dir(local_dir,85),
"HNL100": HNL_from_dir(local_dir, 100),
"HNL125": HNL_from_dir(local_dir, 125),
"HNL150": HNL_from_dir(local_dir, 150),
"HNL200": HNL_from_dir(local_dir, 200),
"HNL250": HNL_from_dir(local_dir, 250),
"HNL300": HNL_from_dir(local_dir, 300),
"HNL350": HNL_from_dir(local_dir, 350),
"HNL400": HNL_from_dir(local_dir, 400),
"HNL450": HNL_from_dir(local_dir, 450),
"HNL500": HNL_from_dir(local_dir,500),
"HNL600": HNL_from_dir(local_dir,600),
"HNL700": HNL_from_dir(local_dir,700),
"HNL800": HNL_from_dir(local_dir,800),
"HNL900": HNL_from_dir(local_dir,900),
"HNL1000": HNL_from_dir(local_dir,1000),

# #"HNL20": HNL_from_dir(local_dir, 20),
# # "HNL30": HNL_from_dir(local_dir, 30),
# # "HNL40": HNL_from_dir(local_dir, 40),
# # "HNL50": HNL_from_dir(local_dir,50),
# # "HNL60": HNL_from_dir(local_dir,60),
# # "HNL70": HNL_from_dir(local_dir,70),
# # "HNL75": HNL_from_dir(local_dir,75),




#Drell-Yann
"DYJets_To_LL_M_50_madgraphMLM": files_from_dir(local_dir_bck_DYJets_To_LL_M_50_madgraphMLM),
# "DYJets_To_LL_M_50_amcatnloFXFX": files_from_dir(local_dir_bck_DYJets_To_LL_M_50_amcatnloFXFX),
# "DY1Jets_To_LL_M_50": files_from_dir(local_dir_bck_DY1Jets_To_LL_M_50),
# "DY2Jets_To_LL_M_50": files_from_dir(local_dir_bck_DY2Jets_To_LL_M_50 ),
# "DY3Jets_To_LL_M_50": files_from_dir(local_dir_bck_DY3Jets_To_LL_M_50 ),
# "DY4Jets_To_LL_M_50": files_from_dir(local_dir_bck_DY4Jets_To_LL_M_50 ),
#"DYJets_To_LL_M_10to50":  files_from_dir(local_dir_bck_DYJets_To_LL_M_10to50 ),
#ELECRTOWEAK
"EWK_WMinus2Jets_W_To_LNu_M_50": files_from_dir(local_dir_bck_EWK_WMinus2Jets_W_To_LNu_M_50),
"EWK_WPlus2Jets_W_To_LNu_M_50": files_from_dir(local_dir_bck_EWK_WPlus2Jets_W_To_LNu_M_50 ),
"EWK_Z2Jets_Z_To_LL_M_50": files_from_dir(local_dir_bck_EWK_Z2Jets_Z_To_LL_M_50),
#TTbar
"TT_To_2L2Nu": files_from_dir(local_dir_bck_TT_To_2L2Nu),
"TT_To_SemiLeptonic": files_from_dir(local_dir_bck_TT_To_SemiLeptonic),
"TT_To_Hadronic": files_from_dir(local_dir_bck_TT_To_Hadronic),
#singletop
"ST_tW_antitop_5f_inclusiveDecays": files_from_dir(local_dir_bck_ST_tW_antitop_5f_inclusiveDecays),
"ST_tW_top_5f_inclusiveDecays": files_from_dir(local_dir_bck_ST_tW_top_5f_inclusiveDecays),
"ST_t_channel_antitop_4f_InclusiveDecays": files_from_dir(local_dir_bck_ST_t_channel_antitop_4f_InclusiveDecays),
"ST_t_channel_top_4f_InclusiveDecays": files_from_dir(local_dir_bck_ST_t_channel_top_4f_InclusiveDecays),
#W+jets
#"WJets_To_LNu": files_from_dir(local_dir_bck_WJets_To_LNu),
"WJ1ets_To_LNu": files_from_dir(local_dir_bck_WJ1ets_To_LNu),
"W2Jets_To_LNu": files_from_dir(local_dir_bck_W2Jets_To_LNu),
"W3Jets_To_LNu": files_from_dir(local_dir_bck_W3Jets_To_LNu),
"W4Jets_To_LNu": files_from_dir(local_dir_bck_W4Jets_To_LNu),
#DiBoson
#"ZZ": files_from_dir(local_dir_bck_ZZ),
"ZZ_To_4L": files_from_dir(local_dir_bck_ZZ_To_4L),
"ZZ_To_2L2Nu": files_from_dir(local_dir_bck_ZZ_To_2L2Nu),
#"WW": files_from_dir(local_dir_bck_WW),
"WW_To_2L2Nu": files_from_dir(local_dir_bck_WW_To_2L2Nu),
#"WZ": files_from_dir(local_dir_bck_WZ),
"WZ_To_3LNu": files_from_dir(local_dir_bck_WZ_To_3LNu),
#TT + bosons
"TTWJets_To_LNu": files_from_dir(local_dir_bck_TTWJets_To_LNu),
"TTZ_To_LLNuNu_M_10": files_from_dir(local_dir_bck_TTZ_To_LLNuNu_M_10),
"TTWW": files_from_dir(local_dir_bck_TTWW ),
 "TTZZ": files_from_dir(local_dir_bck_TTZZ),
"TTWZ": files_from_dir(local_dir_bck_TTWZ),
#tribosons
"WWW_4F": files_from_dir(local_dir_bck_WWW_4F),
"ZZZ": files_from_dir(local_dir_bck_ZZZ),
"WZZ": files_from_dir(local_dir_bck_WZZ),
"WWZ_4F": files_from_dir(local_dir_bck_WWZ_4F),


#"ZZ_old": files_from_dir(local_dir_bck_ZZ_old),

}

if test:
    for k, v in samples.items():
        samples[k] = v[:1]

result = processor.run_uproot_job(
    samples,
    "Events",
    HNLAnalysis_etautau(),
    #processor.futures_executor,
    processor.iterative_executor, # may be better for debugging
    {"schema": NanoAODSchema, 'workers': 6},
)
event_counter = processor.run_uproot_job(
    samples,
    'Runs',
    CountEvents(),
    #processor.futures_executor,
    processor.iterative_executor, # may be better for debugging

    {"schema": NanoAODSchema, 'workers': 6},
)



with open(f'counter_{tag}.pkl', 'wb') as f:
    pickle.dump(event_counter, f)

with open(f'result_{tag}.pkl', 'wb') as f:
    pickle.dump(result, f)


# old_name = 'C:\\Users\\lucas\\Desktop\\PDM\\analysis\\ratio.txt'
# new_name = "C:\\Users\\lucas\\Desktop\\PDM\\analysis\\ratio_"+str(tag)+".txt"
# os.rename(old_name,new_name)
