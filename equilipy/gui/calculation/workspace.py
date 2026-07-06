"""Calculation workspace behavior mixin."""

from __future__ import annotations

from equilipy.gui.calculation.module_input import CalculationModuleInputMixin
from equilipy.gui.calculation.runner import CalculationRunnerMixin
from equilipy.gui.calculation.session_files import CalculationSessionFilesMixin
from equilipy.gui.calculation.session_overview import CalculationSessionOverviewMixin


class CalculationWorkspaceLogicMixin(
    CalculationSessionOverviewMixin,
    CalculationModuleInputMixin,
    CalculationSessionFilesMixin,
    CalculationRunnerMixin,
):
    """Combined calculation workspace behavior."""
