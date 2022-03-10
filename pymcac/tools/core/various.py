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

"""Various little functions."""


def get_idx_name(ds):
    """Infer the name of the index."""
    potentials = set(ds.dims) - ({"k", "Time"} & set(ds.dims))
    if "Num" in potentials:
        return "Num"
    if "Label" in potentials:
        return "Label"
    if len(potentials) == 1:
        return potentials.pop()
    return None
