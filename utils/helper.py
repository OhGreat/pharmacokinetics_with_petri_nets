from enum import Enum
from tkinter import Place

from prometheus_client import PlatformCollector

class PlaceType(Enum):
  P_HOMO = 'Paracetamol Homogenate'
  G_HOMO = 'Paracetamol-glucuronide Homogenate'
  S_HOMO = 'Paracetamol-sulfate Homogenate'
  P_EXC = 'Paracetamol Excreted'
  G_EXC = 'Paracetamol-glucuronide Excreted'
  S_EXC = 'Paracetamol-sulfate Excreted'

"""
Mapping from types to places
"""
places = {
    PlaceType.P_HOMO: 'P',
    PlaceType.S_HOMO: 'S',
    PlaceType.G_HOMO: 'G',
    PlaceType.P_EXC: 'P excreted',
    PlaceType.S_EXC: 'S excreted',
    PlaceType.G_EXC: 'G excreted',
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
}

metabolic = [PlaceType.G_HOMO, PlaceType.P_EXC, PlaceType.S_HOMO]


def get_tokens(net, place):
    return  sum(list(net.place(place).tokens))


def get_clean_tokens():
  for place in PlaceType:
    tokens[place] = []
  return tokens