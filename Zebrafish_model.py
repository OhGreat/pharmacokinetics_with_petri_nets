import snakes.plugins
snakes.plugins.load('gv', 'snakes.nets', 'my_nets')
from snakes.nets import *
from my_nets import *  # required to draw networks


class ZebraMol():
    def __init__(self, kappas=None, t=0):
        if kappas is None:
            # model time dependancy
            k_PS_f_0 = 0.422
            t_50 = 1.42  # time in minutes at which the formation
                         # rate for the sulfate metabolite is at 50%
                         # of its value at time 0.
            k_PS_f = k_PS_f_0*(1- (t/(t_50 - t)))

            self.kappas = { 'k_a': 0.760,
                            'k_PG,f': 0.00327,
                            'k_PS,f': k_PS_f,
                            'k_P,e': 0.0185,
                            'k_G,e': 0.00743,
                            'k_S,e': 0.000664, }
        else: self.kappas = kappas
        self.net = self.create_model()

    def save_img(self, path):
        self.net.draw(path, engine='dot')

    def create_model(self):
        net = PetriNet('Zebrafish Paracetamol net')
        
        net.add_place(Place('input')) # paracetamol in water
        net.add_place(Place('P')) # paracetamol in zebrafish
        add_sequence(net=net,
                    name="P absorption",
                    from_place='input',
                    to_place='P',
                    in_var="x",
                    expr=f"x * {self.kappas['k_a']} ")

        net.add_place(Place('P excreted'))
        add_sequence(net=net,
                    name="P excretion",
                    from_place='P',
                    to_place='P excreted',
                    in_var="x",
                    expr=f"x * {self.kappas['k_P,e']} ")
    
        net.add_place(Place('G'))
        add_sequence(net=net,
                    name="Glucuronidation ",
                    from_place='P',
                    to_place='G',
                    in_var="x",
                    expr=f"x * {self.kappas['k_PG,f']} ")
        
        net.add_place(Place('G excreted'))
        add_sequence(net=net,
                    name="G excretion",
                    from_place='G',
                    to_place='G excreted',
                    in_var="x",
                    expr=f"x *{self.kappas['k_G,e']} ")
        
        net.add_place(Place('S'))
        add_sequence(net=net,
                    name="Sulfation",
                    from_place='P',
                    to_place='S',
                    in_var="x",
                    expr=f"x * {self.kappas['k_PS,f']} ")

        net.add_place(Place('S excreted'))
        add_sequence(net=net,
                    name="S excretion",
                    from_place='S',
                    to_place='S excreted',
                    in_var="x",
                    expr=f"x * {self.kappas['k_S,e']} ")
        return net

def add_sequence(net, name, from_place, to_place, in_var, expr, t_exp=None):
    """ Params:
        - name: anme of transition
        - from_place: start point of transition
        - to_place: end point of transition
        - in_var: input variable for transition, e.g. 'x'
        - expr: expression for the output of the transition, e.g. 'x+1'
    """
    if t_exp is not None:
        net.add_transition(Transition(name), Expression(t_exp))
    else:    
        net.add_transition(Transition(name))
    # create arcs
    net.add_input(from_place, name, Variable(in_var))
    net.add_output(to_place, name, Expression(expr))
