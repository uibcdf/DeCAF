# ! This project has been moved to [gitlab](https://gitlab.com/cheminfIBB/decaf) !

---

# README #

## DeCAF - Discrimination, Comparison, Alignment tool for small molecules ##

DeCAF is a method for describing molecules' pharmacophoric properties and a fast and efective tool for comparing and combining ligands.

DeCAF is written as a Python module and can be easily combined with OpenBabel, RDKit and other chemoinformatic tools.


## Examples ##

### Basic example: ###
```
#!python

#use RDKit
from rdkit.Chem import MolFromSmiles
from decaf.toolkits import rd
from decaf.utils import similarity, combine_pharmacophores, draw

#create models
mol1 = MolFromSmiles("c1cc(cc(c1)N)c1cc(cc(c1O)c1[nH]c2ccc(cc2n1)C(=N)N)Cl")
phar1 = rd.phar_from_mol(mol1)
mol2 = MolFromSmiles("c1cc(cc(c1)N(=O)=O)c1cc(cc(c1O)c1cc2cc(ccc2[nH]1)C(=N)N)CC(=O)O")
phar2 = rd.phar_from_mol(mol2)

#compare and combine models
print similarity(phar1, phar2)
phar = combine_pharmacophores(phar1, phar2)

#draw pharmacophore
draw(phar)
```

### Demos ###
There are also three demo scripts, showing prossible applications of DeCAF:

* compare_demo.py - compare two sets of ligands to refrence set
* filter_demo.py - screen database for molecules similar to given pharmacophore model
* model_demo.py - create pharmacophore models from set of molecules

## Requirements ##
* numpy (for basic functionalities)
* Pybel (OpenBabel) and/or RDKit (for creating models from molecules)
* matplotlib (for pharmacophore drawing)
* scipy (for spring layout)

## Documentation ##
Automatic documentation for DeCAF can be build with Sphinx (http://sphinx-doc.org/)

```
cd docs
make html
make latexpdf
```

## Installation ##
After installing OpenBabel or RDKit, you can install DeCAF with setuptools:
```
python setup.py install
```

...or use pip:

```
pip install decaf
```

You can also install DeCAF and all its dependencies with [conda](http://conda.pydata.org/):

```
#add conda channels for DeCAF, OpenBabel and RDKit
conda config --add channels marta-sd
conda config --add channels mwojcikowski
conda config --add rdkit

#install DeCAF with all dependencies
conda install decaf
```

## Contact ##

Please contact us with comments, suggestions, questions or bug reports:

* Marta M. Stepniewska-Dziubinska ( martasd[at]ibb.waw.pl )

## Citing DeCAF ##

If you use DeCAF in your research, please cite:

Stepniewska-Dziubinska, M. M., Zielenkiewicz, P., & Siedlecki, P. (2017). DeCAF—Discrimination, Comparison, Alignment Tool for 2D PHarmacophores. *Molecules*, 22(7), 1128; [doi:10.3390/molecules22071128](http://dx.doi.org/10.3390/molecules22071128)
