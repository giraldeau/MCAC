#!/usr/bin/env python3
# coding: utf-8

# MCAC
# Copyright (C) 2020 CORIA
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Read the advancement file of MCAC output files
"""

from pathlib import Path
from typing import Union

import dask.dataframe as dd


class AdvancementReader:
    """
    Object containing all functions necessary to read an advancement file
    """

    __slots__ = ("filename",)

    def __init__(self, filename: Union[str, Path]) -> None:
        self.filename = filename

    @classmethod
    def read_advancement(cls, dir: Union[str, Path]) -> dd.DataFrame:
        """
        Read the advancement file
        """
        return cls(Path(dir) / "advancement.dat").advancement

    @property
    def advancement(self) -> dd.DataFrame:
        """
        Read the advancement file
        """
        return dd.read_csv(
            self.filename,
            sep=" ",
            names=["Time", "concentration", "volume_fraction", "avg_npp", "temperature", "BoxSize"],
        ).set_index("Time")
