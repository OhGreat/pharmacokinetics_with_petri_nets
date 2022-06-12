import numpy as np
import snakes.plugins
snakes.plugins.load('gv', 'snakes.nets', 'my_nets')
from snakes.nets import *
from my_nets import *
from Zebrafish_model_2 import *


def main():
    zebra_model = ZebraMol()

    net = zebra_model.net
    #net.add_marking(Marking(P=([1])))
    print("initial_marking:",net.get_marking())
    net.add_marking(Marking(input=([1])))

    p_values = []
    for i in range(300):
        fire_continuous(net, ['P absorption'], verbose=False)
        fire_continuous(net, ['P excretion', 'Glucuronidation', 'Sulfation'], verbose=False)
        fire_continuous(net, ['G excretion'], verbose=False)
        fire_continuous(net, ['S excretion'], verbose=False)

        P_markings = net.place('P')
        p_values.append(list(P_markings.tokens)[0])

    np.save("test.npy",p_values)
    print(net.get_marking())

if __name__ == '__main__':
    main()
