import numpy as np
import snakes.plugins
snakes.plugins.load('gv', 'snakes.nets', 'my_nets')
from snakes.nets import *
from my_nets import *
from Zebrafish_model import *
from enum import Enum


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

def get_tokens(net, place):
    return  sum(list(net.place(place).tokens))

def main(verbose=False):
    zebra_model = ZebraMol()

    net = zebra_model.net
    #net.add_marking(Marking(P=([1])))
    print("initial_marking:",net.get_marking())
    net.add_marking(Marking(input=([10])))

    for t in range(100): # 25
        zebra_model.k_PS_f = zebra_model.k_PS_f_0*(1- (t/(zebra_model.t_50 + t)))
        zebra_model.kappas['k_PS,f'] = zebra_model.k_PS_f
        zebra_model.net.remove_transition('Sulfation')
        add_sequence(net=net,
            name="Sulfation",
            from_place='P',
            to_place='S',
            in_var="x",
            expr=f"x * {zebra_model.kappas['k_PS,f']} ")

        fire_continuous(net, ['P absorption'], verbose=False)
        fire_continuous(net, ['P excretion', 'Glucuronidation', 'Sulfation'], verbose=False)
        fire_continuous(net, ['G excretion'], verbose=False)
        fire_continuous(net, ['S excretion'], verbose=False)
        # Paracetamol Marking
        for pl in PlaceType:
            if verbose:
                print(f'Adding tokens for: {pl.value}')
            tokens[pl].append(get_tokens(net, places[pl]))



    for pl in PlaceType:
        if verbose:
            print(pl.value)
            print(tokens[pl])
        np.save(out[pl], tokens[pl])

    print(net.get_marking())

if __name__ == '__main__':
    main()
