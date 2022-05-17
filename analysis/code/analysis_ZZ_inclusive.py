import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from coffea import hist
import numpy as np
from matplotlib import pyplot as plt 
import mplhep as hep
import uproot3
from ROOT import TFile
plt.style.use([hep.style.CMS])

fname = '/eos/user/l/lmollier/PDM/data/bck_v2/MC_background/ZZ/ZZ.root'
events = NanoEventsFactory.from_root(fname, schemaclass=NanoAODSchema).events()
print("done")
z =  events.GenPart[(np.abs(events.GenPart.pdgId) == 23)]
z = z[z.hasFlags(['isLastCopy'])]
z_child = z.distinctChildren[z.distinctChildren.hasFlags(['isLastCopy'])] 
print(z_child.pdgId)
z_child = ak.flatten(z_child)
print(z_child.pdgId)
print(z_child[:,0].pdgId)
print(len(z_child))
z_child = z_child[ak.num(z_child) == 2]
print(len(z_child))
z_mass_from_child = (z_child[:,0] + z_child[:,1]).mass
print(z_mass_from_child)

mass = z_mass_from_child
#mass = z.mass
#mass = ak.flatten(mass)
mass_axis = hist.Bin("mass", r"$m_{z}$ [GeV]", 300, 0., 150)
h = hist.Hist("Z mass", mass_axis)
h.fill(mass=mass)
fout = uproot3.recreate('inclZZ_massZ.root' )
fout['Z_mass'] = hist.export1d(h)
fout.close()

fout = TFile("inclZZ_massZ.root", "READ")
hnew = fout.Get("Z_mass")
hnew.Draw("HIST")
input("press enter")
fout.Close()
# mass = ak.to_numpy(mass)
# bins = []
# for i in range(0,400):
#     bins.append(i/2)
# plt.hist(mass, bins = bins)
# plt.title("Z mass of inclusive ZZ")
# plt.xlabel(r'$m_{Z}$')
# plt.ylabel("counts")
# plt.show()
