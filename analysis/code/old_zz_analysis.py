import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import numpy as np
from matplotlib import pyplot as plt 

fname = '/eos/user/l/lmollier/PDM/data/bck_v2/MC_background/ZZ_old/nano.root'
events = NanoEventsFactory.from_root(fname, schemaclass=NanoAODSchema).events()
print("done")
z =  events.GenPart[(np.abs(events.GenPart.pdgId) == 23)]
z = z[z.hasFlags(['isLastCopy'])]
#always 2Z in the events
child_z = z.distinctChildren
print(child_z.pdgId)
child_z = child_z[child_z.hasFlags(['isLastCopy'])]
print(child_z.pdgId)
is_e = abs(child_z.pdgId) == 11
is_m = abs(child_z.pdgId) == 13
is_t = abs(child_z.pdgId) == 15

child_z = child_z[~is_e & ~is_m & ~is_t]
print(child_z.pdgId)
all_child = ak.flatten(ak.flatten(child_z))
print(all_child)



