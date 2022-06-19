import snakes.plugins
snakes.plugins.load(['ops', 'gv'], 'snakes.nets', 'my_nets')
from snakes.nets import *
from my_nets import *

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
        return
    # we assume we are starting a fresh experiment.
    init_mark = sum(init_markings)
    # debug print
    if verbose:
        print("initial marking:",init_mark)
    # loop through transitions
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
    if verbose:
        print(net.get_marking())

def update_transition(net, trans_name,  
                    update_fun, t):
    """ Used to model time dependent transitions
        Params:
            - net: SNAKES PetriNet
            - trans_name: name of the transition (string)
            - update_fun: lambda function used to update transition parameter
            - t: current timestep used to update 
    """
    # calculate new transition value multiplier
    new_mul = update_fun(t)
    # current transition
    trans = net.transition(trans_name)
    # get input place
    input_p = list(trans.__dict__['_input'].keys())[0].name
    # get output place
    out_p = list(trans.__dict__['_output'].keys())[0].name
    # create new expression
    in_var = "x"
    expr_new = in_var+f" * {new_mul}"
    # remove old transition
    net.remove_transition(trans_name)
    # create new transition
    add_sequence(net=net,
            name=trans_name,
            from_place=input_p,
            to_place=out_p,
            in_var=in_var,
            expr=expr_new)

def draw_place (place, attr) :
    attr['label'] = place.name.upper()
    attr['color'] = '#FF0000'

def draw_transition (trans, attr) :
    if str(trans.guard) == 'True' :
        attr['label'] = trans.name
    else :
        attr['label'] = '%s\n%s' % (trans.name, trans.guard)