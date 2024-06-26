from aiidalab_widgets_base import (
    ComputationalResourcesWidget,
)


PhonopyCode = ComputationalResourcesWidget(
    description="phonopy:",
    default_calc_job_plugin="phonopy.phonopy",
)
