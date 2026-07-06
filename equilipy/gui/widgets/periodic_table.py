"""Periodic-table element selector widget.

Renders the standard 18-column periodic table as checkable element cells.
Elements available in the current context (for example a loaded database)
are clickable; selection is shown as an accent-colored edge box.  Elements
outside the context are disabled and rendered 70% transparent.  Symbols that
are not real elements (pseudo-elements) fall back to an extra row below the
table so nothing silently disappears.
"""

from __future__ import annotations

from collections.abc import Iterable

from PySide6.QtCore import Signal
from PySide6.QtGui import QPalette
from PySide6.QtWidgets import QGridLayout, QSizePolicy, QToolButton, QWidget

# symbol -> (row, column) in the 18-column layout; rows 0-6 are periods 1-7,
# row 7 is a spacer, rows 8/9 hold the lanthanide/actinide blocks.
_PERIODIC_LAYOUT: dict[str, tuple[int, int]] = {
    "H": (0, 0), "HE": (0, 17),
    "LI": (1, 0), "BE": (1, 1), "B": (1, 12), "C": (1, 13), "N": (1, 14),
    "O": (1, 15), "F": (1, 16), "NE": (1, 17),
    "NA": (2, 0), "MG": (2, 1), "AL": (2, 12), "SI": (2, 13), "P": (2, 14),
    "S": (2, 15), "CL": (2, 16), "AR": (2, 17),
    "K": (3, 0), "CA": (3, 1), "SC": (3, 2), "TI": (3, 3), "V": (3, 4),
    "CR": (3, 5), "MN": (3, 6), "FE": (3, 7), "CO": (3, 8), "NI": (3, 9),
    "CU": (3, 10), "ZN": (3, 11), "GA": (3, 12), "GE": (3, 13),
    "AS": (3, 14), "SE": (3, 15), "BR": (3, 16), "KR": (3, 17),
    "RB": (4, 0), "SR": (4, 1), "Y": (4, 2), "ZR": (4, 3), "NB": (4, 4),
    "MO": (4, 5), "TC": (4, 6), "RU": (4, 7), "RH": (4, 8), "PD": (4, 9),
    "AG": (4, 10), "CD": (4, 11), "IN": (4, 12), "SN": (4, 13),
    "SB": (4, 14), "TE": (4, 15), "I": (4, 16), "XE": (4, 17),
    "CS": (5, 0), "BA": (5, 1), "HF": (5, 3), "TA": (5, 4), "W": (5, 5),
    "RE": (5, 6), "OS": (5, 7), "IR": (5, 8), "PT": (5, 9), "AU": (5, 10),
    "HG": (5, 11), "TL": (5, 12), "PB": (5, 13), "BI": (5, 14),
    "PO": (5, 15), "AT": (5, 16), "RN": (5, 17),
    "FR": (6, 0), "RA": (6, 1), "RF": (6, 3), "DB": (6, 4), "SG": (6, 5),
    "BH": (6, 6), "HS": (6, 7), "MT": (6, 8), "DS": (6, 9), "RG": (6, 10),
    "CN": (6, 11), "NH": (6, 12), "FL": (6, 13), "MC": (6, 14),
    "LV": (6, 15), "TS": (6, 16), "OG": (6, 17),
    "LA": (8, 2), "CE": (8, 3), "PR": (8, 4), "ND": (8, 5), "PM": (8, 6),
    "SM": (8, 7), "EU": (8, 8), "GD": (8, 9), "TB": (8, 10), "DY": (8, 11),
    "HO": (8, 12), "ER": (8, 13), "TM": (8, 14), "YB": (8, 15),
    "LU": (8, 16),
    "AC": (9, 2), "TH": (9, 3), "PA": (9, 4), "U": (9, 5), "NP": (9, 6),
    "PU": (9, 7), "AM": (9, 8), "CM": (9, 9), "BK": (9, 10), "CF": (9, 11),
    "ES": (9, 12), "FM": (9, 13), "MD": (9, 14), "NO": (9, 15),
    "LR": (9, 16),
}

_FALLBACK_ROW = 11
_CELL_SIZE = 40
_UNAVAILABLE_ALPHA = 77  # 30% opacity == 70% transparent


class PeriodicTableSelector(QWidget):
    """Checkable periodic table restricted to an available element set."""

    selectionChanged = Signal()

    def __init__(
        self,
        available: Iterable[str],
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        # canonical (upper) -> original casing, alphabetical for stable output.
        self._available = {
            symbol.upper(): symbol
            for symbol in sorted(available, key=str.upper)
        }
        self.buttons: dict[str, QToolButton] = {}

        grid = QGridLayout(self)
        grid.setSpacing(4)
        grid.setRowMinimumHeight(7, 12)
        grid.setRowMinimumHeight(10, 12)

        fallback_column = 0
        for canonical, position in _PERIODIC_LAYOUT.items():
            grid.addWidget(self._make_cell(canonical), *position)
        for canonical in self._available:
            if canonical not in _PERIODIC_LAYOUT:
                grid.addWidget(
                    self._make_cell(canonical),
                    _FALLBACK_ROW,
                    fallback_column,
                )
                fallback_column += 1

        # Cells are added in table-position order; expose the button map in
        # alphabetical order so callers get a deterministic contract.
        self.buttons = {
            symbol: self.buttons[symbol]
            for symbol in sorted(self.buttons, key=str.upper)
        }

        accent = self.palette().color(QPalette.ColorRole.Highlight)
        accent_name = accent.name()
        self.setStyleSheet(
            "QToolButton {"
            " border: 1px solid rgba(128, 128, 128, 110);"
            " border-radius: 4px;"
            " background: rgba(128, 128, 128, 20);"
            "}"
            "QToolButton:checked {"
            f" border: 2px solid {accent_name};"
            f" background: rgba({accent.red()}, {accent.green()},"
            f" {accent.blue()}, 45);"
            " font-weight: 600;"
            "}"
            "QToolButton:disabled {"
            f" border: 1px solid rgba(128, 128, 128, {_UNAVAILABLE_ALPHA});"
            f" color: rgba(128, 128, 128, {_UNAVAILABLE_ALPHA});"
            " background: rgba(128, 128, 128, 8);"
            "}"
        )

    def _make_cell(self, canonical: str) -> QToolButton:
        original = self._available.get(canonical)
        button = QToolButton(self)
        button.setText(original if original is not None else canonical.title())
        button.setFixedSize(_CELL_SIZE, _CELL_SIZE)
        button.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        font = button.font()
        if font.pointSizeF() > 0:
            font.setPointSizeF(max(font.pointSizeF(), 10.0))
            button.setFont(font)
        if original is None:
            button.setEnabled(False)
        else:
            button.setCheckable(True)
            button.setChecked(True)
            button.toggled.connect(self.selectionChanged)
            self.buttons[original] = button
        return button

    def selected_symbols(self) -> list[str]:
        """Return the checked element symbols in alphabetical order."""
        return [
            symbol
            for symbol, button in self.buttons.items()
            if button.isChecked()
        ]

    def set_all(self, checked: bool) -> None:
        """Check or uncheck every available element."""
        for button in self.buttons.values():
            button.setChecked(checked)
