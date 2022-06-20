from Zebrafish_model import *
from utils.helper import *
from utils.SNAKES_extensions import update_transition, fire_continuous
import numpy as np
from model_adult import AdultModel


class ExperimenterAdult():
    def __init__(
            self, 
            total_timesteps:int, 
            zebra_model:ZebraMol,
            exp_name: str,
            washout: int,
            tokens:dict
        ):
      self.washout = washout
      self.zebra_model = zebra_model
      self.total_timesteps = total_timesteps
      self.exp_name = exp_name
      self.tokens = tokens
      self.net = zebra_model.net

    
    def run_exp(self, verbose=False):
        for t in range(self.total_timesteps):
          if t < self.washout:
              self.net.place('input').empty()
              self.net.add_marking(Marking(input=([500])))
          # if t == self.washout:
              # self.net.place('input').empty()
          # Sulfation is time dependent and needs to be fixed at every timestep
          update_transition(  
                              self.net, "Sulfation",
                              self.zebra_model.k_PS_f, 
                              t
                            )
          # fire our transitions sequentially
          fire_continuous(self.net, ['P absorption'], verbose=False)
          self.tokens[PlaceType.P_HOMO].append(get_tokens(self.net, places[PlaceType.P_HOMO]))
          fire_continuous(self.net, ['P excretion', 'Glucuronidation', 'Sulfation', 'Oxidation'], verbose=False)
          for m in metabolic:
              self.tokens[m].append(get_tokens(self.net, places[m]))
          self.tokens[PlaceType.CYP_HOMO].append(get_tokens(self.net, places[PlaceType.CYP_HOMO]))
          fire_continuous(self.net, ['G excretion'], verbose=False)
          self.tokens[PlaceType.G_EXC].append(get_tokens(self.net, places[PlaceType.G_EXC]))
          fire_continuous(self.net, ['S excretion'], verbose=False)
          self.tokens[PlaceType.S_EXC].append(get_tokens(self.net, places[PlaceType.S_EXC]))
          fire_continuous(self.net, ['CYP excretion'], verbose=False)
          self.tokens[PlaceType.CYP_EXC].append(get_tokens(self.net, places[PlaceType.CYP_EXC]))

        for pl in PlaceType:
          if verbose:
            print(pl.value)
            print(self.tokens[pl])
        
          np.save(f'results/{self.exp_name}/{out[pl]}', self.tokens[pl])
        print(self.net.get_marking())

if __name__ == '__main__':
  # Adult experimenter
  zebra_model = AdultModel()
  exp = ExperimenterAdult(
              total_timesteps=300,
              washout=60,
              zebra_model=zebra_model,
              exp_name='exp3',
              tokens=get_clean_tokens()
          )
  exp.run_exp()