import pathlib
import tempfile


from IPython.display import display

import ipywidgets as ipw

# from ..euphonic.bands_pdos import *


from aiidalab_qe_vibroscopy.utils.euphonic.data.phonopy_interface import (
    generate_force_constant_from_phonopy,
)

###### START for detached app:


# Upload buttons
class UploadPhonopyYamlWidget(ipw.FileUpload):
    def __init__(self, **kwargs):
        super().__init__(
            description="upload phonopy YAML file",
            multiple=False,
            layout={"width": "initial"},
        )


class UploadForceConstantsHdf5Widget(ipw.FileUpload):
    def __init__(self, **kwargs):
        super().__init__(
            description="upload force constants HDF5 file",
            multiple=False,
            layout={"width": "initial"},
        )


class UploadPhonopyWidget(ipw.HBox):
    def __init__(self, **kwargs):
        self.upload_phonopy_yaml = UploadPhonopyYamlWidget(**kwargs)
        self.upload_phonopy_hdf5 = UploadForceConstantsHdf5Widget(**kwargs)

        self.reset_uploads = ipw.Button(
            description="Discard uploaded files",
            icon="pencil",
            button_style="warning",
            disabled=False,
            layout=ipw.Layout(width="auto"),
        )

        super().__init__(
            children=[
                self.upload_phonopy_yaml,
                self.upload_phonopy_hdf5,
                self.reset_uploads,
            ],
            **kwargs,
        )

    def _read_phonopy_files(self, fname, phonopy_yaml_content, fc_hdf5_content=None):
        suffix = "".join(pathlib.Path(fname).suffixes)

        with tempfile.NamedTemporaryFile(suffix=suffix) as temp_yaml:
            temp_yaml.write(phonopy_yaml_content)
            temp_yaml.flush()

            if fc_hdf5_content:
                with tempfile.NamedTemporaryFile(suffix=".hdf5") as temp_file:
                    temp_file.write(fc_hdf5_content)
                    temp_file.flush()
                    temp_hdf5_name = temp_file.name

                    try:
                        fc = generate_force_constant_from_phonopy(
                            path=pathlib.Path(fname),
                            summary_name=temp_yaml.name,
                            fc_name=temp_hdf5_name,
                        )
                    except ValueError:
                        return None

                    return fc
            else:
                temp_hdf5_name = None

                try:
                    fc = generate_force_constant_from_phonopy(
                        path=pathlib.Path(fname),
                        summary_name=temp_yaml.name,
                        # fc_name=temp_hdf5_name,
                    )
                except ValueError:
                    return None

                return fc


#### END for detached app


class DownloadYamlHdf5Widget(ipw.HBox):
    def __init__(self, model):
        self._model = model

        self.download_button = ipw.Button(
            description="Download phonopy data",
            icon="pencil",
            button_style="primary",
            disabled=False,
            layout=ipw.Layout(width="auto"),
        )
        self.download_button.on_click(self._download_data)

        super().__init__(
            children=[
                self.download_button,
            ],
        )

    def _download_data(self, _=None):
        """
        Download both the phonopy.yaml and fc.hdf5 files.
        """
        phonopy_yaml, fc_hdf5 = self._model.produce_phonopy_files()
        self._download(payload=phonopy_yaml, filename="phonopy" + ".yaml")
        self._download(payload=fc_hdf5, filename="fc" + ".hdf5")

    @staticmethod
    def _download(payload, filename):
        from IPython.display import Javascript

        javas = Javascript(
            """
            var link = document.createElement('a');
            link.href = 'data:text/json;charset=utf-8;base64,{payload}'
            link.download = "{filename}"
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
            """.format(payload=payload, filename=filename)
        )
        display(javas)
