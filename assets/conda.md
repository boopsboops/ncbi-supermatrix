## Conda Installation

This script should install the software on a conda environment on a remote server or machine that you don't have sudo for.

```bash
# clone the repository 
git clone https://github.com/boopsboops/ncbi-supermatrix.git
cd ncbi-supermatrix

# recreate the conda environment from the yaml
conda env create -f conda.yaml

# activate conda env
conda activate ncbi-supermatrix

# important! need to spoof some symlinks for shared libs to keep R and conda happy
ln -s $CONDA_PREFIX/include/freetype2/freetype $CONDA_PREFIX/include/freetype
ln -s $CONDA_PREFIX/include/cairo/* $CONDA_PREFIX/include/

# install the R packages
# you are perhaps better off opening R and running 'renv::restore()'
# can also restore individually to diagnose problems better, e.g. renv::restore(packages="ape")
Rscript -e "renv::restore()"

# once you have finished analysis leave conda env
conda deactivate

# useful to nuke the whole env if you need to start fresh
conda remove --name ncbi-supermatrix --all
```
