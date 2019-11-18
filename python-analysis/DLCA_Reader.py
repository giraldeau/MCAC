#!/usr/bin/python3
# coding: utf-8

import numpy as np
import pandas as pd
#from dask import dataframe as dd

from h5py import File as h5file

from pathlib import Path
from lxml import etree
from threading import Thread, RLock
from multiprocessing import Pool, cpu_count

class XdmfReader(object):
    """
    Object containing all functions necessary to read an xdmf file
    """
    def __init__(self, filename):
        self.filename = filename.with_suffix(".xmf")
        self.xml = None
        self.metadata = None
        self.h5_groups = None

    def parse_xml(self):
        """
        Parse xml file
        """
        parser = etree.XMLParser(remove_blank_text=True)
        for data in open(str(self.filename)):
            parser.feed(data)
        self.xml = parser.close()

    @staticmethod
    def bool_from_any(s):
        """convert printable to boolean"""
        return str(s).lower() in ['true', '1', 't', 'y', 'yes', 'oui']

    def extract_metadata(self):
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
                value=float(value)
            metadata[key] = value
        return metadata

    def extract_h5_groups(self):
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

    def read(self):
        """
        Read the xdmf file and return metadata and h5_groups
        """
        self.metadata = self.extract_metadata()
        self.h5_groups = self.extract_h5_groups()
        return self.metadata, self.h5_groups

class H5Reader(object):
    """
    Object containing all functions necessary to read an xdmf file
    """
    def __init__(self, filename):
        self.filename = filename.with_suffix(".h5")

    @classmethod
    def read_file(cls, filename, index):
        return cls(filename).read(index=index)

    def read(self, h5_groups=None, index=None):
        """
        Read the h5 file

        We need info from the xmf file, so this file is also read if info is not given
        """
        if h5_groups is None:
            h5_groups = XdmfReader(self.filename).extract_h5_groups()

        data = {}
        with h5file(str(self.filename), 'r') as file_h5:
            for time, step in h5_groups.items():
                datat = {}
                for attrib, group in step.items():
                    npdata = np.asarray(file_h5[group])
                    if attrib == "Positions":
                        pos = npdata.reshape((-1,3))
                        datat['Posx'] = pos[:,0]
                        datat['Posy'] = pos[:,1]
                        datat['Posz'] = pos[:,2]
                    else:
                        datat[attrib] = npdata
#                datat["Aggkey"] = list(map(lambda x: hash((time,x)),datat["Label"]))

                n = datat['Posx'].size
                datat["Time"] = np.repeat(datat["Time"], n)
                datat["BoxSize"] = np.repeat(datat["BoxSize"], n)

                del datat["Time"]

                data[time] = pd.DataFrame.from_dict(datat)

                if index is not None:
                    data[time].set_index(index, inplace=True)
        return data

class DLCA(object):
    """
    This object read and contains the simulation result from DLCA
    """
    def __init__(self, datadir, seq=False):
        self.dir = datadir
        self.metadata = None
        self.Aggregates = None
        self.Spheres = None
        self.times = None

        self.seq = seq
        if seq:
            self.read_data = self.read_data_seq
        else:
            self.read_data = self.read_data_par
        self.have_data = False

    def read_metadata(self):
        """
        Read the metadata of the simulation from one of the files
        """
        # usually spheres are smaller file to read
        files = list(self.dir.glob("Spheres*.xmf"))

        # usually the last one is smaller
        filename = files[-1]

        self.metadata = XdmfReader(filename).extract_metadata()

    def extract_time(self, data, lock=None):
        """
        Extract time steps from data.
        Save them sorted and return them unsorted (for later reordering)
        """
        times = np.fromiter(data.keys(), dtype=float)

        # save sorted time steps
        if lock is not None:
            with lock:
                if self.times is None:
                    self.times = np.sort(times)
        else:
            if self.times is None:
                self.times = np.sort(times)
        return times

    #@profile
    def read_data_seq(self, files, index=None):
        """
        Read all the data into a dict(dict(ndarray))

        the result is in the form data[time][attribut][label]
        """
        args = [(file, index) for file in files]

        datas = [H5Reader.read_file(arg) for arg in args]

        data = {k: v for d in datas for k, v in d.items() }

        return data

    def read_data_par(self, files, index=None):
        """
        Read all the data into a dict(dict(ndarray))

        the result is in the form data[time][attribut][label]
        """

        args = [ (file, index) for file in files]

        with Pool(processes=cpu_count()) as pool :
            datas = pool.starmap(H5Reader.read_file, args)

            data = { k: v for d in datas for k, v in d.items() }

        return data

    #@profile
    def read_all_aggregates(self, lock=None):
        """
        Read all the data from the aggregates files

        If not already stored, store the time steps

        The result is a large panda multiindex dataframe in the form data[attribut][time][label]
        """
        files = list(self.dir.glob("Aggregats*.xmf"))

        # Read data
        data = self.read_data(files, index='Label')
        times = self.extract_time(data, lock)

        # turn data into a large panda multiindex dataframe
        self.Aggregates = pd.concat(data, names=['Time', 'Label']).sort_index(0)

    #@profile
    def read_all_spheres(self, lock=None):
        """
        Read all the data from the aggregates files

        If not already stored, store the time steps

        The result is a large panda multiindex dataframe in the form data[attribut][time][label]
        """
        files = list(self.dir.glob("Spheres*.xmf"))

        # Read data
        data = self.read_data(files)
        times = self.extract_time(data, lock)

        # turn data into a large panda multiindex dataframe
        self.Spheres = pd.concat(data, names=['Time', 'Num']).sort_index(0)

    def read_seq(self):
        """
        Seqential version of the full reader
        """
        self.read_metadata()
        self.read_all_spheres()
        self.read_all_aggregates()

        return self.metadata, self.times, self.Spheres, self.Aggregates


    def read_par(self):
        """
        Parallel version of the full reader
        """
        self.read_metadata()

        lock = RLock()
        dlca = self
        class AggregatesJob(Thread):
            def run(self):
                dlca.read_all_aggregates(lock)

        aggregates_job = AggregatesJob()
        aggregates_job.start()

        class SphereJob(Thread):
            def run(self):
                dlca.read_all_spheres(lock)

        sphere_job = SphereJob()
        sphere_job.start()
        sphere_job.join()
        aggregates_job.join()

        return self.metadata, self.times, self.Spheres, self.Aggregates


    def read(self, dask=False):
        """
        Read paraview output or preprocessed files if availables
        """
        data = self.dir / "data.h5"
        if not data.exists():
            if self.seq:
                self.read_seq()
            else:
                self.read_par()

            #self.Spheres.to_hdf(data, key="spheres", mode='w', format="table")
            #self.Aggregates.to_hdf(data, key="aggregats", mode='a', format="table")
            self.have_data = True

        if dask:
            #self.Spheres = pd.utils.from_pandas(self.Spheres, len(self.Spheres) / 2**16 + 1)
            #self.Aggregates = pd.utils.from_pandas(self.Aggregates, len(self.Aggregates) / 2**16 + 1)


            #self.Spheres = dd.from_pandas(self.Spheres.reset_index(0), len(self.Spheres) / 2**16 + 1)
            #self.Aggregates = dd.from_pandas(self.Aggregates.reset_index(0), len(self.Aggregates) / 2**16 + 1)

            self.Spheres = dd.read_hdf(str(data), key="spheres")
            self.Aggregates = dd.read_hdf(str(data), key="aggregats")
            self.have_data = True

        if not self.have_data:
            self.Spheres = pd.read_hdf(data, key="spheres")
            self.Aggregates = pd.read_hdf(data, key="aggregats")

        return self.Spheres, self.Aggregates

