# DeCAF

This repository is a fork of [DeCAF](https://bitbucket.org/marta-sd/) by [Marta Stepniewska-Dziubinska et alt.](https://doi.org/10.3390/molecules22071128)

DeCAF (Discrimination, Comparison, Alignment tool for 2D PHarmacophores) is an open-source Python to perform ligand-based virtual screenings based on 2D pharmacophoric models where chemical features -pharmacophoric elements- and distances among them take the shape of a graph. These pharmacophoric graphs are in the core of a set of algorithms to compare, align, and screen large sets of small molecules. DeCAF is freely available at http://bitbucket.org/marta-sd/decaf.

## Citation

[Stepniewska-Dziubinska, Marta M., Piotr Zielenkiewicz, and Pawel Siedlecki. “DeCAF—Discrimination, Comparison, Alignment Tool for 2D PHarmacophores.” Molecules : A Journal of Synthetic Chemistry and Natural Product Chemistry 22, no. 7 (2017): 1128.](https://doi.org/10.3390/molecules22071128)

## Why this fork?

The original DeCAF repository was forked here with an academic purpose only. **If you are planning to use DeCAF for production in the context of a scientific project, please DO NOT USE THIS REPOSITORY**. Use [the original DeCAF's repository](https://bitbucket.org/marta-sd/) and give credit only to the [legitime authors of DeCAF](https://doi.org/10.3390/molecules22071128).

### Changes in this forked version of DeCAF

DeCAF's workflow and strategy were implemented by [Marta M. Stepniewska-Dziubinska, Piotr Zielenkiewicz and Pawel Siedlecki](https://doi.org/10.3390/molecules22071128). The code shown here has some variations in its form, but not in the form. As such, the contributors to this repository are not the intelectual parents of the ideas behind DeCAFs. 

The following is a brief description of the changes made by the contributors to this forked
repository of DeCAF:

- New directory 'devtools' with instructions to create a conda environment to work with DeCAF.
- Every test script was placed in the directory 'tests'.
- The scripts, notebooks and data included as supplementary material of DeCAF's main paper were
  included here in the directory 'supplementary'.

## Installation

### Conda environment with required packages

```
cd DeCAF/devtools/conda-envs
python create_conda_env.py development_env.yaml -n DeCAF -p 3.7
```

### Installing the decaf module

```
conda activate DeCAF
cd DeCAF
python setup.py develop
```

## DeCAF License

The original [DeCAF License](./LICENSE_DeCAF) is included in this repository: 

Copyright (c) 2015, Marta M. Stepniewska
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the {organization} nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



