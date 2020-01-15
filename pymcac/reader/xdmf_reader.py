#!/usr/bin/env python3
# coding: utf-8
"""
Read the xdmf part of the MCAC output files

TODO do not use lxml
"""

from pathlib import Path
from typing import Union, Dict, Optional, Tuple

from lxml import etree


class XdmfReader:
    """
    Object containing all functions necessary to read an xdmf file
    """
    __slots__ = ("filename",
                 "xml",
                 "metadata",
                 "h5_groups")

    def __init__(self, filename: Union[str, Path]) -> None:
        self.filename = Path(filename).with_suffix(".xmf")
        self.xml = None
        self.metadata: Optional[Union[bool, float]] = None
        self.h5_groups: Optional[Dict[float, Dict[str, str]]] = None

    # noinspection PyUnusedFunction
    @classmethod
    def read_file(cls,
                  filename: Union[str, Path]) -> Tuple[Dict[str, Union[bool, float]],
                                                       Dict[float, Dict[str, str]]]:
        """
        Read the xdmf file and return metadata and h5_groups
        """
        return cls(filename).read()

    def parse_xml(self) -> None:
        """
        Parse xml file
        """
        parser = etree.XMLParser(remove_blank_text=True)
        for data in open(str(self.filename)):
            parser.feed(data)
        self.xml = parser.close()

    @staticmethod
    def bool_from_any(s: str) -> bool:
        """convert printable to boolean"""
        return str(s).lower() in ['true', '1', 't', 'y', 'yes', 'oui']

    def extract_metadata(self) -> Dict[str, Union[bool, float]]:
        """
        Extract metadata from xml
        """
        if self.xml is None:
            self.parse_xml()

        metadata = {}
        for element in self.xml.iter("Information"):
            key = element.get("Name")
            if key in ["Copyright", "Physics"]:
                continue
            value = element.get("Value")
            if "Active" in key:
                value = self.bool_from_any(value)
            else:
                value = float(value)
            metadata[key] = value
        return metadata

    def extract_h5_groups(self) -> Dict[float, Dict[str, str]]:
        """
        Extract which h5_group will correspond to which data
        """
        if self.xml is None:
            self.parse_xml()

        h5_groups = {}
        for step in self.xml.iter("Grid"):
            if step.get("Name") == "Collection":
                continue
            time = float(step.find('Time').get("Value"))

            h5_groups.setdefault(time, dict())
            h5_groups[time]["Positions"] = step.find('Geometry').getchildren()[0].text.split(":")[-1]
            for attrib in step.findall("Attribute"):
                key = attrib.get("Name")
                h5_groups[time][key] = attrib.getchildren()[0].text.split(":")[-1]
        return h5_groups

    def extract_sizes(self) -> Dict[float, int]:
        """
        Extract which h5_group will correspond to which data
        """
        if self.xml is None:
            self.parse_xml()

        sizes = {}
        for step in self.xml.iter("Grid"):
            if step.get("Name") == "Collection":
                continue
            time = float(step.find('Time').get("Value"))

            for attrib in step.findall("Attribute"):
                if attrib.get("Name") in ("BoxSize", "Time"):
                    continue
                size = attrib.getchildren()[0].get("Dimensions")
                sizes[time] = int(size)
                break
        return sizes

    def read(self) -> Tuple[Dict[str, Union[bool, float]],
                            Dict[float, Dict[str, str]]]:
        """
        Read the xdmf file and return metadata and h5_groups
        """
        self.metadata = self.extract_metadata()
        self.h5_groups = self.extract_h5_groups()
        return self.metadata, self.h5_groups
