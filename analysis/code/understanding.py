import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

fname = "data_test.root"



events = NanoEventsFactory.from_root(fname, schemaclass=NanoAODSchema).events()
events = events[2:5]
muons = events.Muon
electrons = events.Electron
print(electrons)
print(muons)
cart = ak.cartesian([muons,electrons])
print(cart)
cart_n = ak.cartesian([muons,electrons], nested = True)
print(cart_n)