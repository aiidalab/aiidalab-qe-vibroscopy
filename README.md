# aiidalab-qe-vibroscopy
Plugin to compute vibrational properties of materials via the aiida-vibroscopy AiiDA plugin

## Installation

Once cloned the repository, `cd` into it and:

```shell
pip install --user .
```

If you want to easily set up phonopy, use the CLI of this package (inspect it via `aiidalab-qe-vibroscopy --help`):

```shell.
pip install phonopy --user # if phonopy is not installed in your machine; it should be already installed as it is a dependency of the package.
aiidalab-qe-vibroscopy setup-phonopy # setup phonopy@localhost in AiiDA; this post-install command is automatically triggered if you install the plugin from the aiidalab-qe interface.
```

### Specific details for arm64 architectures

#### Installation of scipy from conda is required

In case of `arm64` architecture, please run `aiidalab-qe-vibroscopy setup-phonopy`.
This will install the correct version of scipy to work with the package.

#### h5py issues
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
