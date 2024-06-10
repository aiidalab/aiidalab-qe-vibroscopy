# aiidalab-qe-vibroscopy
Plugin to compute vibrational properties of materials via the aiida-vibroscopy AiiDA plugin

## Installation

In order to load this plugin into QeApp, you need to switch to the `features/phonopy` branch for your `aiidalab-qe`.

Then, install this plugin by:

```shell
pip install -e .
pip install phonopy --user # if phonopy is not installed in your machine
ln -s <phonopy path> /opt/conda/bin/phonopy # if it doesn't already exist, create a symbolic link for phonopy
phonopy_setup.py # setup phonopy in AiiDA
```

### h5py installation for arm64 architectures

In case the installation of the `aiidalab-qe-vibroscopy` fails due to `h5py` installation problem, you may try to first install
`h5py` by:

```shell
conda install h5py==3.11.0
```

this will install also the `hdf5` library as dependency.


## License

MIT

## Contact

miki.bonacci@psi.ch
andres.ortega-guerrero@empa.ch
