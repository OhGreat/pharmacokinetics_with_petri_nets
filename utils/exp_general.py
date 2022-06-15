from Zebrafish_model import *
from utils.helper import *
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
        print(self.washout)
        for t in range(self.total_timesteps):
            if t < self.washout:
                self.net.place('input').empty()
                self.net.add_marking(Marking(input=([1])))
            if t == self.washout:
                self.net.place('input').empty()

            self.zebra_model.k_PS_f = self.zebra_model.k_PS_f_0*(1-(t/(self.zebra_model.t_50 + t)))
            self.zebra_model.kappas['k_PS,f'] = self.zebra_model.k_PS_f
            self.net.remove_transition('Sulfation')
            add_sequence(net=self.net,
                name="Sulfation",
                from_place='P',
                to_place='S',
                in_var="x",
                expr=f"x * {self.zebra_model.kappas['k_PS,f']} ")

            fire_continuous(self.net, ['P absorption'], verbose=False)
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
        
            np.save(f'results/{self.exp_name}/{out[pl]}', self.tokens[pl])
        print(self.net.get_marking())
