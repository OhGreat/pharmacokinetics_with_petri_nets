import snakes.plugins
snakes.plugins.load('gv', 'snakes.nets', 'nets')
from nets import *

def addSequence(net, input_name, output_name, transition_name, cons, prod, expression=None):
    if expression is None:
        net.add_transition(Transition(transition_name))
    else:
        net.add_transition(Transition(transition_name, expression))
    net.add_input(input_name, transition_name, cons)
    net.add_output(output_name, transition_name, prod)

def createNet(cons, prod, init=[]):
    n = PetriNet('Paracetamol Zebrafish')
    n.add_place(Place('i', init))
    n.add_place(Place('p', []))
    addSequence(n, 'i', 'p', 'Ka', cons, prod)

    n.add_place(Place('Pe', []))
    addSequence(n, 'p', 'Pe', 'KPe', cons, prod)

    n.add_place(Place('G', []))
    addSequence(n, 'p', 'G', 'KPG, f', cons, prod)

    n.add_place(Place('S', []))
    addSequence(n, 'p', 'S', 'KPS, f', cons, prod)

    n.add_place(Place('Ge', []))
    addSequence(n, 'G', 'Ge', 'KG,e', cons, prod)

    n.add_place(Place('Se', []))
    addSequence(n, 'S', 'Se', 'KS,e', cons, prod)

    return n

if __name__ == '__main__':
    net = createNet(Variable('x'), Variable('x'))
    net.draw('zebrafish.png', engine='dot')
