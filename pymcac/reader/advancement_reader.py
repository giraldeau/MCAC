#!/usr/bin/env python3

# MCAC
# Copyright (C) 2020 CORIA
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Read the advancement file of MCAC output files."""

from pathlib import Path
from typing import Union

import pandas as pd


class AdvancementReader:
    """All functions necessary to read an advancement file."""

    __slots__ = ("filename",)

    def __init__(self, filename: Union[str, Path]) -> None:
        """Init."""
        self.filename = filename

    @classmethod
    def read_advancement(cls, folder: Union[str, Path]) -> pd.DataFrame:
        """Read the advancement file."""
        return cls(Path(folder) / "advancement.dat").advancement

    @property
    def advancement(self) -> pd.DataFrame:
        """Read the advancement file."""
        df = (
            pd.read_csv(
                self.filename,
                sep=" ",
                names=[
                    "Time",
                    "concentration",
                    "volume_fraction",
                    "avg_npp",
                    "temperature",
                    "BoxVolume",
                    "monomer_concentration",
                    "u_sg",
                    "flux_nucleation",
                ],
            )
            .groupby(by="Time", sort=True)
            .last()
        )

        return df
