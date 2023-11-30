"""Implementation of the VibroWorkchain for managing the aiida-vibroscopy workchains."""
from aiida.common import AttributeDict
from aiida.engine import ToContext, WorkChain, calcfunction
from aiida.orm import AbstractCode, Int, Float, Dict, Code, StructureData, load_code
from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.utils.mapping import prepare_process_inputs
from aiida_quantumespresso.common.types import ElectronicType, SpinType
from aiida.engine import WorkChain, calcfunction, if_
from aiida_vibroscopy.common.properties import PhononProperty


IRamanSpectraWorkChain = WorkflowFactory("vibroscopy.spectra.iraman")
HarmonicWorkChain = WorkflowFactory("vibroscopy.phonons.harmonic")
DielectricWorkChain = WorkflowFactory("vibroscopy.dielectric")
PhononWorkChain = WorkflowFactory("vibroscopy.phonons.phonon")

class VibroWorkChain(WorkChain):
    "WorkChain to compute vibrational property of a crystal."
    label = "vibro"

    @classmethod
    def define(cls, spec):
        """Specify inputs and outputs."""
        super().define(spec)
        spec.input('structure', valid_type=StructureData) #Maybe not needed as input... just in the protocols. but in this way it is not easy to automate it in the app, after the relaxation. So let's keep it for now. 

        spec.expose_inputs(
            PhononWorkChain,
            namespace='phonon',
            exclude=('clean_workdir'), #AAA check this... maybe not needed.
            namespace_options={
                'required': False,
                'populate_defaults': False,
                'help': 'Inputs for the `PhononWorkChain`, triggered selecting "trigger"="phonon".',
            }
        )
        spec.expose_inputs(
            DielectricWorkChain, 
            namespace='dielectric',
            exclude=('clean_workdir'), #AAA check this... maybe not needed.
            namespace_options={
                'required': False,
                'populate_defaults': False,
                'help': (
                    'Inputs for the `DielectricWorkChain` that will be'
                    'used to calculate the mixed derivatives with electric field.'
                )
            },
            #exclude=('symmetry')
        )
        spec.expose_inputs(
            HarmonicWorkChain,
            namespace='harmonic',
            exclude=('clean_workdir'), #AAA check this... maybe not needed.
            namespace_options={
                'required': False,
                'populate_defaults': False,
                'help': 'Inputs for the `HarmonicWorkChain`.',
            }
        )
        spec.expose_inputs(
            IRamanSpectraWorkChain,
            namespace='iraman',
            exclude=('clean_workdir'), #AAA check this... maybe not needed.
            namespace_options={
                'required': False,
                'populate_defaults': False,
                'help': 'Inputs for the `IRamanSpectraWorkChain`.',
            }
        )
        ###
        spec.outline(
            cls.setup,
            cls.vibrate,
            cls.results,
        )
        ###
        spec.expose_outputs(
            PhononWorkChain, namespace='phonon',
            namespace_options={'required': False, 'help':'Outputs of the `PhononWorkChain`.'},
        )
        spec.expose_outputs(
            DielectricWorkChain, namespace='dielectric',
            namespace_options={'required': False, 'help':'Outputs of the `DielectricWorkChain`.'},
        )
        spec.expose_outputs(
            HarmonicWorkChain, namespace='harmonic',
            namespace_options={'required': False, 'help':'Outputs of the `HarmonicWorkChain`.'},
        )
        spec.expose_outputs(
            IRamanSpectraWorkChain, namespace='iraman',
            namespace_options={'required': False, 'help':'Outputs of the `IRamanSpectraWorkChain`.'},
        )
        ###
        spec.exit_code(400, 'ERROR_WORKCHAIN_FAILED', message='The workchain failed.')
    
    @classmethod
    def get_builder_from_protocol(
        cls,
        pw_code,
        structure,
        protocol=None,
        phonopy_code=None,
        overrides=None,
        options=None,
        trigger=None,
        phonon_property=PhononProperty.NONE,
        dielectric_property=None,
        **kwargs
    ):
        """Return a builder prepopulated with inputs selected according to the chosen protocol.

        :param pw_code: the ``Code`` instance configured for the ``quantumespresso.pw`` plugin.
        :param structure: the ``StructureData`` instance to use.
        :param phonopy_code: the ``Code`` instance configured for the ``phonopy.phonopy`` plugin.
        :param protocol: protocol to use, if not specified, the default will be used.
        :param overrides: optional dictionary of inputs to override the defaults of the protocol.
        :param options: A dictionary of options that will be recursively set for the ``metadata.options`` input of all
            the ``CalcJobs`` that are nested in this work chain.
        :param kwargs: additional keyword arguments that will be passed to the ``get_builder_from_protocol`` of all the
            sub processes that are called by this workchain.
        :return: a process builder instance with all inputs defined ready for launch.
        """
        from aiida_quantumespresso.workflows.protocols.utils import recursive_merge
        
        if trigger not in ["phonon","dielectric","harmonic","iraman",]:
            raise ValueError('trigger not in "phonon","dielectric","harmonic","iraman"') 
        
        builder = cls.get_builder()
        
        
        if trigger == "phonon":
            
            builder_phonon = PhononWorkChain.get_builder_from_protocol(
                pw_code=pw_code,
                phonopy_code=phonopy_code,
                structure=structure,
                protocol=protocol,
                overrides=overrides["phonon"],
                phonon_property=phonon_property,
                **kwargs
            )
            
            builder_phonon.phonopy.metadata.options.resources = {
                "num_machines": 1,
                "num_mpiprocs_per_machine": 1,
            }

            #MBO: I do not understand why I have to do this, but it works
            symmetry = builder_phonon.pop('symmetry')
            builder.phonon = builder_phonon
            builder.phonon.symmetry = symmetry
        
        elif trigger == "dielectric":
            builder_dielectric = DielectricWorkChain.get_builder_from_protocol(
                code=pw_code,
                structure=structure,
                protocol=protocol,
                overrides=overrides["dielectric"],
                **kwargs
            )

            #MBO: I do not understand why I have to do this, but it works. maybe related with excludes.
            symmetry = builder_dielectric.pop('symmetry')
            builder.dielectric = builder_dielectric
            builder.dielectric.symmetry = symmetry
            builder.dielectric.property = dielectric_property
        
        elif trigger == "harmonic":
            builder_harmonic = HarmonicWorkChain.get_builder_from_protocol(
                pw_code=pw_code,
                phonopy_code=phonopy_code,
                structure=structure,
                protocol=protocol,
                overrides=overrides,
                phonon_property=phonon_property,
                **kwargs
            )
        
            
            # MB supposes phonopy will always run serially, otherwise choose phono3py 
            # also this is needed to be set here.
            builder_harmonic.phonopy.metadata.options.resources = {
                    "num_machines": 1,
                    "num_mpiprocs_per_machine": 1,
                }

            builder_harmonic.phonon.phonopy.metadata.options.resources = builder_harmonic.phonopy.metadata.options.resources

            #should be automatic inside HarmonicWorkchain.
            builder_harmonic.phonon.phonopy.parameters = Dict(dict={})
            builder_harmonic.phonopy.parameters = builder_harmonic.phonon.phonopy.parameters
            builder_harmonic.phonon.phonopy.code = builder_harmonic.phonopy.code
            
            builder_harmonic.phonopy.parameters = Dict(dict=phonon_property.value)

            builder.harmonic = builder_harmonic
            
        elif trigger == "iraman":
            builder_iraman = IRamanSpectraWorkChain.get_builder_from_protocol(
                code=pw_code,
                structure=structure,
                protocol=protocol,
                overrides=overrides,
                **kwargs
            )
            
            builder_iraman.dielectric.property = dielectric_property
    
            builder_iraman.phonon.phonopy.code = phonopy_code
            builder_iraman.phonon.phonopy.metadata.options.resources = {
                "num_machines": 1,
                "num_mpiprocs_per_machine": 1,
            }
            builder_iraman.phonon.phonopy.parameters = Dict({})
            
            #Remove kpoints_parallel_distance if for 1D and 2D materials
            if structure.pbc != (True, True, True):
                builder_iraman.dielectric.pop("kpoints_parallel_distance", None)

            builder.iraman = builder_iraman


        
        for wchain in ["phonon","dielectric","harmonic","iraman",]:
            if trigger != wchain: builder.pop(wchain,None)
        
        builder.structure = structure
        
        return builder
        
    def setup(self):
        #setup general contest variables... see in HarmonicWorkChain.
        #also see self.ctx.key = 'phonon'... maybe you can initialise here to simplify next functions.
        #key, class, outputs namespace.
        if "phonon" in self.inputs:
            self.ctx.key = "phonon"
            self.ctx.workchain = PhononWorkChain
        elif "dielectric" in self.inputs:
            self.ctx.key = "dielectric"
            self.ctx.workchain = DielectricWorkChain
        elif "harmonic" in self.inputs:
            self.ctx.key = "harmonic"
            self.ctx.workchain = HarmonicWorkChain
        elif "iraman" in self.inputs:
            self.ctx.key = "iraman"
            self.ctx.workchain = IRamanSpectraWorkChain
            
        self.ctx.structure = self.inputs.structure

    
    def vibrate(self):
        """Run a WorkChain for vibrational properties."""
        #maybe we can unify this, thanks to a wise setup.
        inputs = AttributeDict(self.exposed_inputs(self.ctx.workchain, namespace=self.ctx.key))
        inputs.metadata.call_link_label = self.ctx.key
        if self.ctx.key in ["phonon","dielectric"]:
            inputs.scf.pw.structure = self.ctx.structure
        else:
            inputs.structure = self.ctx.structure

        future = self.submit(self.ctx.workchain, **inputs)
        self.report(f'submitting `WorkChain` <PK={future.pk}>')
        self.to_context(**{self.ctx.key: future})

    def results(self):
        """Inspect all sub-processes."""
        workchain = self.ctx[self.ctx.key]

        if not workchain.is_finished_ok:
            self.report(f'the child WorkChain with <PK={workchain.pk}> failed')
            return self.exit_codes.ERROR_WORKCHAIN_FAILED

        self.out_many(self.exposed_outputs(self.ctx[self.ctx.key], self.ctx.workchain, namespace=self.ctx.key))
