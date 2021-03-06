import sys
sys.path.append('../')
from models.Zebrafish_model import *
from models.Zebrafish_model_no_head import ZebraMolNoHead
from exp_general import *
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
  # paper reproduction experiments
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

  # no head experiment
  zebra_model = ZebraMolNoHead()
  exp = Experimenter(
            total_timesteps=300,
            washout=None,
            zebra_model=zebra_model,
            exp_name='exp_no_head',
            tokens=get_clean_tokens()
        )
  exp.run_exp()
