# Dynamic learning rates for continual unsupervised learning

This is the source code used to get the results in the paper "Dynamic learning rates for continual unsupervised learning", acepted for publication in the journal "Integrated Computer-Aided Engineering".

The code is written by Jose David Fernandez Rodriguez, building upon previous code for competitive clustering techniques (in a classic learning setting, i.e. before devising techniques for continual unsupervised learning) by Esteban Palomo Ferrer and Ezequiel Lopez Rubio.

The code is written as a set of Matlab scripts; it was executed in `Matlab 2022a`. The main entry points to run experiments are `DemoCORA.m` and `DemoCORACrossValidation.m`, while `clusteringsCORA.m` is used to run classic clustering techniques. The main entry points to crunch data from the experiments are `getBestCORAResults.m`, `showExperimentsBatches.m`, `showCORAClusters.m` and `showCORAAndNetwork.m`, `writeExperimentResults.m`, `writeExperimentResultsBatches.m`, `writeExperimentResultsNoBatches.m`, and `searchGoodEnough.m`. A good way to get a sense of the dependencies between the scripts is to use the project analyzer in the Matlab GUI. The file `run.py` was part of an effort to run minibatch K-means on the dataset used in the paper, as a comparison with our continual learning method.

The dataset used in the paper is available in the `cora` folder. The `images` folder contains images used as a dataset to perform an earlier, preliminary set of experiments on continual techniques learning.
