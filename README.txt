NRamps:
This repository containds a Python 3 implementation of # NRamps, a multi-module tool for Antimicrobial Peptides prediction. Module AmpPro predicts C-terminus cleavage site, and returns a list of predicted mature sequences. Module NaiveBayesNRamps performs a sequence-dependent prediction of AMPs. Module ChemAtrClassNR predicts AMPs calculating diferent physicochemical attributes.

Dependencies:
    linux/OSX 64X
    Anaconda/Miniconda (https://conda.io/docs/user-guide/install/index.html), for enviroment instalation


INSTALLATION:
    from this folder:
        conda env create -f NRampsEnv.yml
        
USAGE:
    activate the NRamps enviroment:
        $ source activate NRamps
    from this folder execute:
        $ python modules/master.py -m prueba.fa

    
