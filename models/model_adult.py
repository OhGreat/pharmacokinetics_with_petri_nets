from utils.SNAKES_extensions import add_sequence
import snakes.plugins
snakes.plugins.load(['gv', 'ops'], 'snakes.nets', 'my_nets')
from snakes.nets import *
from my_nets import *  # required to draw networks


class AdultModel():
    def __init__(self, kappas=None, t=0):
        if kappas is None:
            # model time dependancy
            self.k_PS_f_0 = 0.422
            self.t_50 = 120

            # k_PS_f is time dependant
            self.k_PS_f = lambda t: self.k_PS_f_0 * (1 - (t/(self.t_50 + t)))

            self.kappas = { 'k_a': 0.80,
                            'k_PG,f': 0.565,
                            'k_PS,f': self.k_PS_f(t),
                            'k_PC,f': 0.0925,
                            'k_PC,e': 0.132,
                            'k_P,e': 0.0154,
                            'k_G,e': 0.00743,
                            'k_S,e': 0.000664, }
            print(self.kappas)
        else: self.kappas = kappas
        self.net = self.create_model()

    def create_model(self):
        net = PetriNet('Zebrafish Paracetamol net')
        
        net.add_place(Place('input')) # paracetamol in water
        net.add_place(Place('P')) # paracetamol in zebrafish (homogenate)
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

        net.add_place(Place('CYP'))
        add_sequence(net=net,
                    name="Oxidation",
                    from_place='P',
                    to_place='CYP',
                    in_var="x",
                    expr=f"x * {self.kappas['k_PC,f']} ")
        
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

        net.add_place(Place('CYP excreted'))
        add_sequence(net=net,
            name="CYP excretion",
            from_place='CYP',
            to_place='CYP excreted',
            in_var="x",
            expr=f"x * {self.kappas['k_PC,e']} ")

        return net

    def save_img(self, path="zebrafish_model.png"):
        self.net.draw(path)

if __name__ == '__main__':
    zebra_model = AdultModel()
    zebra_model.save_img('temp/adult_model.png')
