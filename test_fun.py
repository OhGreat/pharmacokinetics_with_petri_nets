import snakes.plugins
snakes.plugins.load('gv', 'snakes.nets', 'my_nets')
from snakes.nets import *
from my_nets import *
from Zebrafish_model import *


def main():
    zebra_model = ZebraMol()

    net = zebra_model.net
    #net.add_marking(Marking(P=([1])))
    print("initial_marking:",net.get_marking())
    #fire_continuous(net,['P excretion', 'Glucuronidation', 'Sulfation'])
    net.add_marking(Marking(input=([1])))
    fire_continuous(net, ['P absorption'], verbose=True)
    #net.add_marking(Marking(input=([1])))
    #fire_continuous(net, ['P absorption'], verbose=True)



if __name__ == '__main__':
    main()
