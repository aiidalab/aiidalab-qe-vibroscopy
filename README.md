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

## Selected-atom mode participation

The Raman, IR and phonon DOS result panels accept one-based atom indices in
input-structure order, for example `1 3 6..8` or `(1,3; 5. 8 10 .. 30 44)`.
Ranges are inclusive; duplicates are counted once. Invalid or out-of-range
indices are reported before changing the plot. Leave the selection empty to
restore the original curves.

For Raman and IR, a thick line shows the spectrum weighted by the selected
atoms' **mode participation**. For a mode with Cartesian displacements
$u_{a\nu}$, the weight is

$$
w_\nu(S) =
\frac{\sum_{a\in S} M_a |u_{a\nu}|^2}
     {\sum_a M_a |u_{a\nu}|^2}.
$$

This recovers the mass-weighted eigenvector convention used for phonon PDOS.
Weights are averaged over frequency-degenerate modes (using the backend's
$10^{-5}$ THz tolerance) to avoid dependence on arbitrary rotations within
a degenerate subspace. Each optical mode intensity is multiplied by this
weight **before broadening**, using the total curve's normalization.
Selecting every atom therefore reproduces the total; complementary selections
add back to it. This describes the character of the modes carrying the
optical signal, not an additive decomposition into atomic IR/Raman
intensities, whose amplitudes can interfere.

Phonon DOS selections sum the saved atomic PDOS and report it per input unit
cell. This requires automatic per-atom PDOS and verifiable primitive-to-unit
cell mapping; if these are unavailable, the panel explains why atom selection
cannot be applied. Phonon dispersion curves are not atom projected.

Both panels offer a selected-only view and include the selection and its
definition in their JSON downloads. Existing results can be used without
rerunning the original calculation.

## License

MIT

## Contact

miki.bonacci@psi.ch
andres.ortega-guerrero@empa.ch
