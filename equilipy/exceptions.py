"""Custom exception types used by Equilipy."""

from __future__ import annotations


class EquilibError(Exception):
    """Raise when equilibrium calculation fails."""

    def __init__(self, message="Equilibrium calculation failed."):
        self.message = message
        super().__init__(self.message)


class DatabaseParsingError(Exception):
    """Raise when database parsing fails."""

    def __init__(self, message="Failed to parse the database"):
        self.message = message
        super().__init__(self.message)


class DatabaseLoadError(Exception):
    """Raise when a parsed database cannot be loaded into global state."""

    def __init__(self, message="Failed to load the database"):
        self.message = message
        super().__init__(self.message)


class InputConditionError(Exception):
    """Raise when input conditions cannot be interpreted."""

    def __init__(self, message="Input condition is not recognized properly"):
        self.message = message
        super().__init__(self.message)


class PostProcessError(Exception):
    """Raise when post-processing fails."""

    def __init__(self, message="Failed processing calcualtion result."):
        self.message = message
        super().__init__(self.message)


class TransitionError(Exception):
    """Raise when transition search cannot produce a valid transition."""

    def __init__(self, message="Failed to find a phase transition."):
        self.message = message
        super().__init__(self.message)
