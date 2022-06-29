# Bio-modeling: Usage of PetriNets for Pharmacokinetics modeling

## Intoduction
This project aims to model pharmacokinetics of paracetamol in Zebrafish larvae with Petri Nets.

Following the original publication[[1]](#1), translation of the proposed NONMEM model using PetriNets has been achieved, in an attempt to reproduce the results present in the paper.

___

## Experiments

There are two main types of experimental procedures. In order to run the respective experiments you can use the following commands from the main directory of your  cloned repository:
```
cd experiments
python run_experiments.py
```
The above also runs the adapted model with no head (input place). 


In order to run the extended model for adult humans and represent the enzymatic maturation process, the following instructions should be used: 
```
cd experiments
python run_adult.py
```

## Experiment visualization and statistics
Visualization can be achieved using the provided notebooks, which will create all of the graphs using the retrived experimental results. The notebooks can be found under the `visualizations` directory.

___
## References
<a id="1">[1]</a> 
Mechanistic and Quantitative Understanding of
Pharmacokinetics in Zebrafish Larvae through Nanoscale Blood
Sampling and Metabolite Modeling of Paracetamols - 
Rob C. Van Wijk (2019). 
