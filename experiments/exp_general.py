from models.Zebrafish_model import *
from exp_general import *
from utils.helper import *
from utils.SNAKES_extensions import update_transition, fire_continuous
import numpy as np

class Experimenter():
    def __init__(
            self, 
            total_timesteps:int, 
            washout:int, 
            zebra_model:ZebraMol,
            exp_name: str,
            tokens:dict
        ):
      self.zebra_model = zebra_model
      self.total_timesteps = total_timesteps
      self.washout = washout
      self.exp_name = exp_name
      self.tokens = tokens
      self.net = zebra_model.net

    def __del__(self):
        print("Object deleted")
    
    def run_exp(self, verbose=False):
        for t in range(self.total_timesteps):
            if self.washout is not None and t < self.washout:
                self.net.place('input').empty()
                self.net.add_marking(Marking(input=([1])))
            if self.washout is not None and t == self.washout:
                self.net.place('input').empty()
            # Sulfation is time dependent and needs to be fixed at every timestep
            update_transition(  self.net, "Sulfation",
                                self.zebra_model.k_PS_f, t)
            # washout + no head model rules
            if self.washout is not None:
                fire_continuous(self.net, ['P absorption'], verbose=False)
            # no head model
            elif self.washout is None and t == 0:
                self.net.add_marking(Marking(P=([1])))
            # fire our transitions sequentially
            self.tokens[PlaceType.P_HOMO].append(get_tokens(self.net, places[PlaceType.P_HOMO]))
            fire_continuous(self.net, ['P excretion', 'Glucuronidation', 'Sulfation'], verbose=False)
            for m in metabolic:
                self.tokens[m].append(get_tokens(self.net, places[m]))
            fire_continuous(self.net, ['G excretion'], verbose=False)
            self.tokens[PlaceType.G_EXC].append(get_tokens(self.net, places[PlaceType.G_EXC]))
            fire_continuous(self.net, ['S excretion'], verbose=False)
            self.tokens[PlaceType.S_EXC].append(get_tokens(self.net, places[PlaceType.S_EXC]))

        for pl in PlaceType:
            if verbose:
                print(pl.value)
                print(self.tokens[pl])
        
            np.save(f'../results/{self.exp_name}/{out[pl]}', self.tokens[pl])
        print(self.net.get_marking())
