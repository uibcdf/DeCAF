# README

## About

This repository contains data and scripts that were used to test a new method for representing and comparing small molecules.
In this approach a molecule is represented as an undirected graph, in which nodes correspond to atoms with pharmacophoric properties, and edges of the graph represent distances between features.
Two models are compared by finding their maximal common subgraph.
Consequently, this approach combines the benefits of a conformation-free representation of a molecule with the additional spatial and physicochemical information.


The method was implemented as a Python module called [DeCAF](http://bitbucket.org/marta-sd/decaf).
Here we present its performance in ligand-based drug discovery.
We used DeCAF, [USRCAT](https://bitbucket.org/aschreyer/usrcat)[1], and all fingerprints implemented in [Open Babel](http://openbabel.org/)[2] in a methodology based on the Similarity Ensemble Approach (SEA)[3],[4] and compared their behaviour.

You can view the notebooks as HTML files or run them by yourself with Jupyter.



**References:**

[1] Adrian M Schreyer and Tom Blundell. "USRCAT: real-time ultrafast shape recognition with pharmacophoric constraints." *J Cheminf*, 4:27, 2012.

[2] Noel M OLBoyle, Michael Banck, Craig A James, Chris Morley, Tim Vandermeersch, and Geoffrey R Hutchison. "Open babel: An open chemcal toolbox." *J Cheminf*, 3:33, 2011.

[3] Michael J Keiser, Bryan L Roth, Blaine N Armbruster, Paul Ernsberger, John J Irwin, and Brian K Shoichet. "Relating protein pharmacology by ligand chemistry." *Nat. Biotechnol.*, 25(2):197–206, 2007.

[4] Eugen Lounkine, Michael J Keiser, Steven Whitebread, Dmitri Mikhailov, Jacques Hamon, Jeremy L Jenkins, Paul Lavan, Eckhard Weber, Allison K Doak, Serge Côté, et al. "Large-scale prediction and testing of drug activity on side-effect targets." *Nature*, 486(7403):361–367, 2012.

------


## Requirements

To generate and analyze the data we used:

* python 2.7.12
* numpy 1.11.1
* openbabel development version (conda build 2.3.90.4git.9dd91b from [mwojcikowski](http://anaconda.org/mwojcikowski/openbabel) channel)
* pandas 0.18.1
* scipy 0.17.1
* seaborn 0.7.0
* decaf 2.0.0 (available through [PyPI](http://pypi.python.org/pypi/DeCAF) and [conda](http://anaconda.org/marta-sd/decaf), or directly from [Bitbucket repository](http://bitbucket.org/marta-sd/decaf))
* jupyter notebook 4.2.1



## Prepare environment

You can re-create the environment that we used in our study with [miniconda](http://conda.pydata.org/miniconda.html).
After installing miniconda, you should add it to your PATH:

```
export PATH=/path/to/miniconda/bin:$PATH
```

Then, you can create python environment with all packages needed:

```
conda config --add channels marta-sd
conda config --add channels mwojcikowski
conda config --add channels rdkit
conda create -n "decaf_tests" python=2.7.12 notebook=4.2.1 numpy=1.11.1 openbabel=2.3.90.4git.9dd91b pandas=0.18.1 scipy=0.17.1 seaborn=0.7.0 decaf=2.0.0
```

To use the environment you need to activate it:

```
source activate decaf_tests
```
If you want to use USRCAT, install it from [Bitbucket repository](https://bitbucket.org/aschreyer/usrcat).

Now you are ready to run our scripts and notebooks.


## Run scripts and notebooks

To get the scripts, notebooks, and data used in our study you need to clone this repository:

```
git clone http://bitbucket.org/marta-sd/decaf-supplementary.git
```


### Scripts
* **scores_decaf.py**, **scores_usrcat**, and **scores_tc.py** - computes similarity scores between drugs and target (i.e. its active ligands) with DeCAF, USRCAT or Tanimoto coefficient, respectively.
* **random_set.py** - creates random datasets used in the statistical analysis.
* **stats_decaf.py**, **stats_usrcat.py**, and **stats_tc.py** - computes similarity scores for random dataset with DeCAF, USRCAT, and Tanimoto coefficient, respectively. Data were used to find the background distributions for raw scores. The distribution provided parameters for calculating Z-scores, p-values and E-values for predictions computed with the *4_results.ipynb* and *5_usrcat.ipynb* notebooks.


### Notebooks
* **1_dataset.ipynb** - finds active ligands for all targets.
* **2_drug_target_data.ipynb** - combines activity data from different sources and finds ChEMBL IDs for drugs used in the study and their extremely close analogues. It also finds trivial hits in the dataset and low activity hits in ChEMBL data (not used in dataset describing a target).
* **3_stats.ipynb** - finds background distribution that is used in statistical analysis of the results. It uses similarity scores for random molecules. Results for random samples were computed with *stats_decaf.py* and *stats_tc.py*, and samples were generated with *random_set.py*.
* **4_results.ipynb** - computes raw scores, Z-scores, p-values and E-values for the data generated with *scores_decaf.py* and *scores_tc.py*.
* **5_usrcat.ipynb** - finds background distribution and generate the results for USRCAT (similar to notebooks *3_stats.ipynb* and *4_results.ipynb*)
* **6_analysis.ipynb** - analyzes and explains the results.

------

## Contact

Please contact us with comments, suggestions, questions or bug reports:

* Marta M. Stepniewska ( martasd[at]ibb.waw.pl )
