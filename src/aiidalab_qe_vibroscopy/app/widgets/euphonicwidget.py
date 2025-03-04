import ipywidgets as ipw
from IPython.display import display

from aiidalab_qe.common.widgets import LoadingWidget

from aiidalab_qe_vibroscopy.app.widgets.structurefactorwidget import (
    EuphonicStructureFactorWidget,
)
from aiidalab_qe_vibroscopy.app.widgets.euphonicmodel import EuphonicResultsModel

from aiidalab_qe.common.infobox import InAppGuide

##### EUPHONIC WIDGET TO DISPLAY EVERYTHING: UPLOAD, PLOTS, DOWNLOAD #####


class EuphonicWidget(ipw.VBox):
    """
    Widget that will include everything,
    from the upload widget to the tabs with single crystal, powder and Q-plane results.

    The first render() is not the real rendering, is just the rendering of the initialize analysis button.
    The real rendering is done by the _render_for_real method. The reason is that it can take a while,
    so we don't want to block the entire app. In the future, this could happen in a thread.

    PLEASE NOTE: the EuphonicResultsModel which are initialized are actually three, and are the models for the three
    different types of calculations: single crystal, powder, and Q-plane, each of them corresponding to a different
    EuphonicStructureFactorWidget, the true visualization widget for the INS results.
    """

    def __init__(
        self,
        model: EuphonicResultsModel,
        node=None,
        detached_app=False,
        **kwargs,
    ):
        """
        Initialize the Euphonic widget class.

        model : EuphonicResultsModel
            The model containing the results for Euphonic calculations.
        node : optional
            The node associated with the Euphonic calculations, default is None.
        detached_app : bool, optional
            Flag indicating if the app is running in detached mode, default is False.
        **kwargs : dict
            Additional keyword arguments.

        _model : EuphonicResultsModel
            The model containing the results for Euphonic calculations.
        rendered : bool
            Flag indicating if the widget has been rendered.
        loading_widget : LoadingWidget
            Widget indicating loading of INS data.
            Widget for uploading phonopy files (only if detached_app is True).
        download_widget : DownloadYamlHdf5Widget
            Widget for downloading YAML and HDF5 files.

        Methods:
        --------
        render():
            Render just the *plot_button* if it has not been rendered yet.
        _render_for_real(change=None):
            Initialize and render the widgets for single crystal, powder, and Q-plane views.
        _on_reset_uploads_button_clicked(change):
            Reset the upload widgets and disable the plot button.
        _on_upload_yaml(change):
            Handle the upload of a phonopy YAML file.
        _on_upload_hdf5(change):
            Handle the upload of a phonopy HDF5 file.
        """

        super().__init__()

        self._model = model  # this is the single crystal model.
        if not hasattr(self._model, "vibro"):
            self._model.vibro = node

        # For the detached app (i.e. when the widget is used outside the QeApp),
        self._model.detached_app = detached_app
        self._model.fc_hdf5_content = None

        self.rendered = False

    def render(self):
        if self.rendered:
            return

        # tabs creation
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
            # we are in QeApp, for sure data are already there.
            self.plot_button.disabled = False
        else:
            # we are in detached app, we need to upload the files.
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
            InAppGuide(identifier="ins-results"),
            self.plot_button,
            self.tab_widget,
            self.download_widget,
            self.loading_widget,
        )

        self.rendered = True

    def _render_for_real(self, change=None):
        """True render method for the widget(s)

        It fetches the data, creates the models for the second and third different types of
        results (i.e. powder, and Q-plane), and then creates the three widgets for the
        single crystal, powder, and Q-plane, respectively.
        """

        # rendering transition
        self.plot_button.layout.display = "none"
        self.loading_widget.layout.display = "block"

        # fetch the data
        self._model.fetch_data()

        # create the models for the other two types of results.
        powder_model = EuphonicResultsModel(spectrum_type="powder")
        qsection_model = EuphonicResultsModel(spectrum_type="q_planes")

        # setting the data for the other two models, because these are exactly the same.
        # the difference is in the post processiong routines.
        for data in ["fc", "q_path"]:
            setattr(powder_model, data, getattr(self._model, data))
            setattr(qsection_model, data, getattr(self._model, data))

        # NOTE:
        # (1) the vibro node is the same for all the models.
        # (2) the self._model is the one for the single crystal.
        # (3) the powder_model and qsection_model are the ones for the powder and Q-plane views.
        self.tab_widget.children = (
            EuphonicStructureFactorWidget(
                node=self._model.vibro,
                model=self._model,
                spectrum_type="single_crystal",
            ),
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
        # method employed in the detached app to reset the upload widgets.
        self.upload_widget.upload_phonopy_yaml.value.clear()
        self.upload_widget.upload_phonopy_yaml._counter = 0
        self.upload_widget.upload_phonopy_hdf5.value.clear()
        self.upload_widget.upload_phonopy_hdf5._counter = 0

        self.plot_button.layout.display = "block"
        self.plot_button.disabled = True

        self.tab_widget.children = ()

        self.tab_widget.layout.display = "none"

    def _on_upload_yaml(self, change):
        # detached app method to handle the upload of a phonopy YAML file.
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
        # detached app method to handle the upload of a phonopy HDF5 file.
        if change["new"] != change["old"]:
            for fname in self.upload_widget.children[1].value.keys():
                self._model.fc_hdf5_content = self.upload_widget.children[1].value[
                    fname
                ]["content"]


class DownloadYamlHdf5Widget(ipw.HBox):
    """
    Widget to download the phonopy.yaml and fc.hdf5 files.

    The download button will trigger the download of the phonopy.yaml and fc.hdf5 files
    via the _download_data method of the widget and the produce_phonopy_files method of the model.
    In the future, this could be done in a thread, and in a way to upload the files in the phonon
    visualizer of MaterialsCloud.
    """

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
