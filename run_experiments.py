import sys
from Zebrafish_model import *
from utils.exp_general import *
from utils.helper import *

experiments = [
  {
    'total_timesteps': 180,
    'washout': sys.maxsize,
    'exp_name': 'exp1'
  },
  {
    'total_timesteps': 300,
    'washout': 60,
    'exp_name': 'exp2'
  }
]

if __name__ == '__main__':
  for experiment in experiments:
    zebra_model = ZebraMol()
    exp = Experimenter(
                total_timesteps=experiment['total_timesteps'],
                washout=experiment['washout'],
                zebra_model=zebra_model,
                exp_name=experiment['exp_name'],
                tokens=get_clean_tokens()
            )
    exp.run_exp()
    del exp, zebra_model
