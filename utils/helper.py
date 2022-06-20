from enum import Enum

class PlaceType(Enum):
  P_HOMO = 'Paracetamol Homogenate'
  G_HOMO = 'Paracetamol-glucuronide Homogenate'
  S_HOMO = 'Paracetamol-sulfate Homogenate'
  P_EXC = 'Paracetamol Excreted'
  G_EXC = 'Paracetamol-glucuronide Excreted'
  S_EXC = 'Paracetamol-sulfate Excreted'
  CYP_HOMO = 'Paracetamol-oxidation'
  CYP_EXC = 'Paracetamol-oxidation excreted'

"""
Mapping from types to places
"""
places = {
    PlaceType.P_HOMO: 'P',
    PlaceType.S_HOMO: 'S',
    PlaceType.G_HOMO: 'G',
    PlaceType.CYP_HOMO: 'CYP',
    PlaceType.P_EXC: 'P excreted',
    PlaceType.S_EXC: 'S excreted',
    PlaceType.G_EXC: 'G excreted',
    PlaceType.CYP_EXC: 'CYP excreted',
}
"""
Mapping from types to tokens
"""
tokens = {
    PlaceType.P_HOMO: [],
    PlaceType.S_HOMO: [],
    PlaceType.G_HOMO: [],
    PlaceType.P_EXC: [],
    PlaceType.S_EXC: [],
    PlaceType.G_EXC: [],
    PlaceType.CYP_HOMO: [],
    PlaceType.CYP_EXC: [],
}
"""
Mapping from type to output file
"""
out = {
    PlaceType.P_HOMO: 'p_val.npy',
    PlaceType.S_HOMO: 's_val.npy',
    PlaceType.G_HOMO: 'g_val.npy',
    PlaceType.P_EXC: 'pe_val.npy',
    PlaceType.S_EXC: 'se_val.npy',
    PlaceType.G_EXC: 'ge_val.npy',
    PlaceType.CYP_HOMO: 'cyp_val.npy',
    PlaceType.CYP_EXC: 'cype_val.npy',
}

metabolic = [PlaceType.G_HOMO, PlaceType.P_EXC, PlaceType.S_HOMO]


def get_tokens(net, place):
    return  sum(list(net.place(place).tokens))


def get_clean_tokens():
  for place in PlaceType:
    tokens[place] = []
  return tokens