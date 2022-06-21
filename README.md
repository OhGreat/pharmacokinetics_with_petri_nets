# Bio-modeling: Usage of PetriNets for Pharmacokinetics modeling

## Intoduction
Following the original paper[[1]](#1), implementation of the proposed model using PetriNets has been achieved, in an attempt to reproduce the results (initially using NONMEM).

___

## Instructions

There are two main types of experimental procedures. In order to run the respective experiments:
```
python run_experiments.py
```
The above also runs the adapted model with no head (input place). 


In order to run extended model for adult humans and represent enzymatic maturity, the following should be used. 
```
python run_adult.py
```

Visualization can be achieved using the provided notebook, which will create all of the graphs using the retrived experimental results. The notebook is 'visualization.ipynb', 'model_no_head_res.ipynb', 'adult_model_visualization.iypnb'

___
## References
<a id="1">[1]</a> 
Mechanistic and Quantitative Understanding of
Pharmacokinetics in Zebrafish Larvae through Nanoscale Blood
Sampling and Metabolite Modeling of Paracetamols - 
Rob C. Van Wijk (2019). 
