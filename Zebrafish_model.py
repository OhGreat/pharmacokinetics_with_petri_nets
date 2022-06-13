import snakes.plugins
snakes.plugins.load(['gv', 'ops'], 'snakes.nets', 'my_nets')
from snakes.nets import *
from my_nets import *  # required to draw networks

def draw_place (place, attr) :
    attr['label'] = place.name.upper()
    attr['color'] = '#FF0000'

def draw_transition (trans, attr) :
    if str(trans.guard) == 'True' :
        attr['label'] = trans.name
    else :
        attr['label'] = '%s\n%s' % (trans.name, trans.guard)


class ZebraMol():
    def __init__(self, kappas=None, t=0):
        if kappas is None:
            # model time dependancy
            self.k_PS_f_0 = 0.422
            self.t_50 = 2  # time in minutes at which the formation
                         # rate for the sulfate metabolite is at 50%
                         # of its value at time 0.  (1.42)
            self.k_PS_f = self.k_PS_f_0*(1- (t/(self.t_50 + t)))

            self.kappas = { 'k_a': 0.760,
                            'k_PG,f': 0.00327,
                            'k_PS,f': self.k_PS_f,
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
                    name="Glucuronidation",
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

        # add_sequence(net=net,
        #             name="S reabsorption",
        #             from_place='S excreted',
        #             to_place='input',
        #             in_var="x",
        #             expr=f"x")
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

def fire_continuous(net: PetriNet, transitions, verbose=False):
    
    # used to update 
    init_mark_final = 0

    # take starting place     
    start_place = [curr_in for curr_in in net.transition(transitions[0]).__dict__['_input']][0]
    # all initial markings
    init_markings = list(start_place.tokens)
    # if no markings, return
    if len(init_markings) == 0: 
        # print('here')
        return

    init_mark = sum(init_markings)

    if verbose:
        print("initial marking:",init_mark)

    for transition in transitions:
        curr_trans = net.transition(transition) 

        # take the output place
        output_place = [curr_in for curr_in in curr_trans.__dict__['_output']][0]
        # intial marking of output place
        output_mark_init = sum(list(output_place.tokens))

        # fire transition
        modes = net.transition(transition).modes()
        net.transition(transition).fire(modes[0])

        # Get output mark final
        output_mark_final = sum(list(output_place.tokens))
        init_mark_final += (output_mark_final - output_mark_init)
        # reset start place marking for next transitions
        start_place.add([init_mark])
        output_place.tokens = MultiSet([sum(list(output_place.tokens))])

        
    # remove the last added transition from the above loop
    start_place.remove([init_mark])    
    # fix start place marking value
    start_place.add([init_mark - init_mark_final])
    # Add all tokens of output

    if verbose:
        print(net.get_marking())


if __name__ == '__main__':
    net = ZebraMol()
    net.save_img('temp/model.png')
