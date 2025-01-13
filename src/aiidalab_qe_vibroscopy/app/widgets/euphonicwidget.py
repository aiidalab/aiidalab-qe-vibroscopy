import ipywidgets as ipw
from IPython.display import display

from aiidalab_qe.common.widgets import LoadingWidget

from aiidalab_qe_vibroscopy.app.widgets.structurefactorwidget import (
    EuphonicStructureFactorWidget,
)
from aiidalab_qe_vibroscopy.app.widgets.euphonicmodel import EuphonicResultsModel


##### START OVERALL WIDGET TO DISPLAY EVERYTHING:


class EuphonicWidget(ipw.VBox):
    """
    Widget that will include everything,
    from the upload widget to the tabs with single crystal and powder predictions.
    In between, we trigger the initialization of plots via a button.
    """

    def __init__(
        self,
        model: EuphonicResultsModel,
        node=None,
        detached_app=False,
        **kwargs,
    ):
        """
        Initialize the Euphonic utility class.
        Parameters:
        -----------
        mode : str, optional
            The mode of operation, default is "aiidalab-qe app plugin".
        fc : optional
            Force constants, default is None.
        q_path : optional
            Q-path for phonon calculations, default is None. If Low-D system, this can be provided.
            It is the same path obtained for the PhonopyCalculation of the phonopy_bands.
        Attributes:
        -----------
        mode : str
            The mode of operation.
        upload_widget : UploadPhonopyWidget
            Widget for uploading phonopy files.
        fc_hdf5_content : None
            Content of the force constants HDF5 file.
        tab_widget : ipw.Tab
            Tab widget for different views.
        plot_button : ipw.Button
            Button to initialize INS data.
        fc : optional
            Force constants if provided.
        """

        super().__init__()

        self._model = model  # this is the single crystal model.
        self._model.vibro = node
        self._model.detached_app = detached_app
        self._model.fc_hdf5_content = None

        self.rendered = False

    def render(self):
        if self.rendered:
            return

        self.tab_widget = ipw.Tab()
        self.tab_widget.layout.display = "none"
        self.tab_widget.set_title(0, "Single crystal")
        self.tab_widget.set_title(1, "Powder sample")
        self.tab_widget.set_title(2, "Q-plane view")
        self.tab_widget.children = ()

        self.plot_button = ipw.Button(
            description="Initialise INS data",
            icon="pencil",
            button_style="primary",
            disabled=True,
            layout=ipw.Layout(width="auto"),
        )
        self.plot_button.on_click(self._render_for_real)

        self.loading_widget = LoadingWidget("Loading INS data")
        self.loading_widget.layout.display = "none"

        if not self._model.detached_app:
            self.plot_button.disabled = False
        else:
            from aiidalab_qe_vibroscopy.utils.euphonic.detached_app.uploadwidgets import (
                UploadPhonopyWidget,
            )

            self.upload_widget = UploadPhonopyWidget()
            self.upload_widget.reset_uploads.on_click(
                self._on_reset_uploads_button_clicked
            )
            self.upload_widget.children[0].observe(self._on_upload_yaml, "_counter")
            self.upload_widget.children[1].observe(self._on_upload_hdf5, "_counter")
            self._model.upload_widget = self.upload_widget
            self.children += (self.upload_widget,)

        self.download_widget = DownloadYamlHdf5Widget(model=self._model)
        self.download_widget.layout.display = "none"

        self.children += (
            self.plot_button,
            self.tab_widget,
            self.download_widget,
            self.loading_widget,
        )

        self.rendered = True

    def _render_for_real(self, change=None):
        # It creates the widgets
        self.plot_button.layout.display = "none"
        self.loading_widget.layout.display = "block"

        self._model.fetch_data()  # should be in the model, but I can do it here once for all and then clone the model.
        powder_model = EuphonicResultsModel(spectrum_type="powder")
        qsection_model = EuphonicResultsModel(spectrum_type="q_planes")

        for data in ["fc", "q_path"]:
            setattr(powder_model, data, getattr(self._model, data))
            setattr(qsection_model, data, getattr(self._model, data))

        # I first initialise this widget, to then have the 0K ref for the other two.
        # the model is passed to the widget. For the other two, I need to generate the model.
        singlecrystalwidget = EuphonicStructureFactorWidget(
            node=self._model.vibro, model=self._model, spectrum_type="single_crystal"
        )

        self.tab_widget.children = (
            singlecrystalwidget,
            EuphonicStructureFactorWidget(
                node=self._model.vibro, model=powder_model, spectrum_type="powder"
            ),
            EuphonicStructureFactorWidget(
                node=self._model.vibro, model=qsection_model, spectrum_type="q_planes"
            ),
        )

        for widget in self.tab_widget.children:
            widget.render()  # this is the render method of the widget.

        self.loading_widget.layout.display = "none"
        self.tab_widget.layout.display = "block"
        self.download_widget.layout.display = "block"

    def _on_reset_uploads_button_clicked(self, change):
        self.upload_widget.upload_phonopy_yaml.value.clear()
        self.upload_widget.upload_phonopy_yaml._counter = 0
        self.upload_widget.upload_phonopy_hdf5.value.clear()
        self.upload_widget.upload_phonopy_hdf5._counter = 0

        self.plot_button.layout.display = "block"
        self.plot_button.disabled = True

        self.tab_widget.children = ()

        self.tab_widget.layout.display = "none"

    def _on_upload_yaml(self, change):
        if change["new"] != change["old"]:
            for fname in self.upload_widget.children[
                0
            ].value.keys():  # always one key because I allow only one file at the time.
                self._model.fname = fname
                self._model.phonopy_yaml_content = self.upload_widget.children[0].value[
                    fname
                ]["content"]

        if self.plot_button.disabled:
            self.plot_button.disabled = False

    def _on_upload_hdf5(self, change):
        if change["new"] != change["old"]:
            for fname in self.upload_widget.children[1].value.keys():
                self._model.fc_hdf5_content = self.upload_widget.children[1].value[
                    fname
                ]["content"]


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
