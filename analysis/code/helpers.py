import os
import re
import numpy as np

def delta_r2(v1, v2):
    '''Calculates deltaR squared between two particles v1, v2 whose
    eta and phi methods return arrays
    '''
    dphi = (v1.phi - v2.phi + np.pi) % (2 * np.pi) - np.pi
    deta = v1.eta - v2.eta
    dr2 = dphi**2 + deta**2
    return dr2

def delta_r(v1, v2):
    '''Calculates deltaR between two particles v1, v2 whose
    eta and phi methods return arrays.
    
    Note: Prefer delta_r2 for cuts.
    '''
    return np.sqrt(delta_r2(v1, v2))

def delta_phi(v1, v2):
    return (v1.phi - v2.phi + np.pi) % (2 * np.pi) - np.pi

def delta_eta(v1, v2):
    return (v1.eta - v2.eta + np.pi) % (2 * np.pi) - np.pi

def cos_opening_angle(v, pv, p1, p2):
    '''Calculates cosine of opening angle between passed
    vertex v (with methods x, y, z), primary vertex pv,
    and the four-vector sum of two passed particles
    p1 and p2 (with methods pt/eta/phi)
    '''
    x = v.vtx_x - pv.x
    y = v.vtx_y - pv.y
    z = v.vtx_z - pv.z
    px = (p1.pt * np.cos(p1.phi) + p2.pt * np.cos(p2.phi))
    py = (p1.pt * np.sin(p1.phi) + p2.pt * np.sin(p2.phi))
    pz = (p1.pt * np.sinh(p1.eta) + p2.pt * np.sinh(p2.eta))
    
    num = x*px + y*py + z*pz
    den = np.sqrt(x**2 + y**2 + z**2) * np.sqrt(px**2 + py**2 + pz**2)

    return num/den

def dxy_significance(v, pv):
    dxy2 = (v.vtx_x - pv.x)**2 + (v.vtx_y - pv.y)**2
    edxy2 = v.vtx_ex**2 + v.vtx_ey**2 #PV has no but anyway negligible error
    return np.sqrt(dxy2/edxy2)

def inv_mass(p1, p2):
    '''Calculates invariant mass from 
    two input particles in pt/eta/phi
    representation, assuming zero mass 
    (so it works for track input)'''
    px1 = p1.pt * np.cos(p1.phi)
    px2 = p2.pt * np.cos(p2.phi)
    px = px1 + px2

    py1 = p1.pt * np.sin(p1.phi) 
    py2 = p2.pt * np.sin(p2.phi)
    py = py1+ py2
    pz1 = p1.pt * np.sinh(p1.eta)
    pz2 = p2.pt * np.sinh(p2.eta)
    pz = pz1 + pz2
    #e = np.sqrt(p1.pt**2 + pz1**2) + np.sqrt(p2.pt**2 + pz2**2)
    e = np.sqrt(px1**2 + py1**2 + pz1**2) + np.sqrt(px2**2 + py2**2 + pz2**2)
    e2 = e**2
    m2 = e2 - px**2 - py**2 - pz**2
    
    if(m2 < 0):
        return 0
    else:
        return np.sqrt(m2)

def inv_mass_3p(p1, p2, p3):
    '''Calculates invariant mass from 
    three input particles in pt/eta/phi
    representation, assuming zero mass 
    (so it works for track input)'''
    x = p1.pt * np.cos(p1.phi) + p2.pt * np.cos(p2.phi) + p3.pt * np.cos(p3.phi)
    y = p1.pt * np.sin(p1.phi) + p2.pt * np.sin(p2.phi) + p3.pt * np.sin(p3.phi)
    z = p1.pt * np.sinh(p1.eta) + p2.pt * np.sinh(p2.eta) + p3.pt * np.sinh(p3.eta)
    e = np.sqrt(p1.pt**2 + (p1.pt * np.sinh(p1.eta))**2) + np.sqrt(p2.pt**2 + (p2.pt * np.sinh(p2.eta))**2) + np.sqrt(p3.pt**2 + (p3.pt * np.sinh(p3.eta))**2)
    m2 = e**2 - x**2 - y**2 - z**2
    return np.sqrt(m2)

def files_from_dir(d):
    '''Returns a list of all ROOT files in the passed directory.
    '''
    files = os.listdir(d)
    return ['/'.join([d, f]) for f in files if f.endswith('.root')]


def HNL_from_dir(d, mass):
    '''Returns a list of all ROOT files in the passed directory.
    '''
    files = os.listdir(d)
    return ['/'.join([d, f]) for f in files if f.endswith(str(mass)+'.root')]

def one_file_from_dir(d):
    '''Returns a list of all ROOT files in the passed directory.
    '''
    files = os.listdir(d)
    files = files[:1]
    return ['/'.join([d, f]) for f in files if f.endswith('.root')]



def files_from_dirs(dirs):
    '''Returns a list of all ROOT files from the directories in the passed list.
    '''
    files = []
    for d in dirs:
        files += files_from_dir(d)
    return files

