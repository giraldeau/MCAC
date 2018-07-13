
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import ctypes
import os
from pathlib import Path
from time import time
from queue import Empty
from scipy.special import erf
from mpl_toolkits.mplot3d import Axes3D
from skimage import measure

from DLCA_Reader import DLCA
from Analyze_cython import coverages_cython, label_argsort

WithSBL = True


def SBLError():
    print("SBL not available, autocorrelation no possible")


try:
    from sbl_wrapper import compute_volume_surface_sbl
except ImportError:
    WithSBL = False
    SBLError()


def random_vector():
    """ Return a random vector on the unit sphere"""
    theta = np.random.uniform(0, 2 * np.pi)
    phi = np.arccos(1 - 2 * np.random.uniform(0, 1))

    return np.array([np.sin(phi) * np.cos(theta),
                     np.sin(phi) * np.sin(theta),
                     np.cos(phi)])


def translated_union(vec, spheres):
    """
        Compute the union of an aggegate and a translated copy of itself
    """
    nspheres = spheres.shape[0]

    # duplication
    duplicate = np.vstack((spheres, spheres))

    # translation of one version
    duplicate[:nspheres, :3] += vec

    return duplicate


def single_volume_autoco(i, radius, spheres):
    """
        Compute the volume needed for the autocorrelation

        Except for the initial volume (i == 0),
        compute the union of an aggegate and a translated copy of itself
        in a random direction

    """
    if i == 0:
        aggregate = spheres
    else:
        vec = radius * random_vector()
        aggregate = translated_union(vec, spheres)

    # compute volume with sbl
    return compute_volume_surface_sbl(aggregate)[0]


def surface_volume_sbl(Agg, Spheres):
    """
        Compute surface and volume using discretisation

        alphagangue is a parameter allowing some gangue around the aggregate
    """
    t, label = Agg.name

    cols = ["Posx", "Posy", "Posz", "Radius"]

    # get the spheres of the given aggregate
    spheres = Spheres[cols][Spheres["Label"] == label].loc[t].values.copy()

    return compute_volume_surface_sbl(spheres)


class Advancement(object):
    """
        A context to print the advancement during computation
    """
    def __init__(self, njobs):
        self.njobs = njobs
        self.last_percent = 0
        self.start = 0.
        self.last_time = self.start
        self.current = 0

    def __enter__(self):
        """
            Start the timer
        """
        self.start = time()
        self.last_time = self.start
        print("{:4.0f}%".format(0), end="")
        return self

    def print(self, it):
        """
            Print the current advancement
        """
        self.current = it
        current_percent = (100 * it) // self.njobs
        current_time = time()
        if (current_percent > self.last_percent or
            current_time - self.last_time > 1):

            if it > 0:
                eta = (self.njobs - it) * (current_time - self.start) / it
            else:
                eta = np.inf

            print("\r{:4.0f}% - ETA:{:4.0f}s".format(current_percent, eta),
                  end="")

            self.last_percent = current_percent
            self.last_time = current_time

    def __exit__(self, type, value, traceback):
        """
            Print the total time
        """
        end = time()
        current_percent = (100 * self.current) // self.njobs
        print("\r{:4.0f}% - Time:{:4.0f}s".format(current_percent,
                                                  end - self.start))


def par_volumes_autoco(q, r, spheres):
    """
        This code will be executed by each process

        Pick one job and do it
    """
    first = True

    # Each process must have a different random seed
    np.random.seed()

    while True:
        try:
            # get work
            i, radius = q.get(block=first)

            # work
            volume = single_volume_autoco(i, radius, spheres)

            # put result
            r.put((i, volume))

            first = False
        except Empty:
            if q.empty() and q.qsize() == 0:
                # nothing else to do, shutting
                break


class Parallel(object):
    """
        A context to do parallel computation
    """
    def __init__(self, target, array, nprocs=None):
        self.target = target
        self.nprocs = nprocs

        if self.nprocs is None:
            self.nprocs = len(os.sched_getaffinity(0))

        # input data to shared numpy array
        self._shared_array = mp.Array(ctypes.c_double, int(array.size))
        self.shared_array = np.frombuffer(self._shared_array.get_obj())
        self.shared_array = self.shared_array.reshape(array.shape)
        self.shared_array[:] = array

        self.qin = mp.Queue()
        self.qout = mp.Queue()

        self.processes = [mp.Process(target=target,
                                     args=(self.qin, self.qout,
                                           self.shared_array))
                          for i in range(self.nprocs)]

    def __enter__(self):
        """
            Start the processes
        """
        for p in self.processes:
            p.start()
        return self

    def __exit__(self, type, value, traceback):
        """
            End the processes
        """
        # first kill the, just in case
        for p in self.processes:
            p.terminate()

        for p in self.processes:
            p.join()


def parallel_volumes_autoco(spheres, lradius, nsamples, nprocs=None):
    """ Compute all the volumes needed for the autocorrelation

        including the volume of the initial aggregate
    """
    # prepare output array
    volumes = np.zeros(len(lradius))

    # useless to make a contiguous copy, Process will do it
    with Parallel(par_volumes_autoco, spheres, nprocs) as p:

        # Put work inside the queue
        njobs = 0
        for i, radius in enumerate(lradius):
            for j in range(nsamples):
                p.qin.put((i, radius))
                njobs += 1

                # Only compute the volume of the initial aggregate once
                if i == 0:
                    break

        with Advancement(njobs) as advancement:
            # Recieve result
            for k in range(njobs):
                i, volume = p.qout.get()
                volumes[i] += volume

                advancement.print(k)

    return volumes


def seq_volumes_autoco(spheres, lradius, nsamples, nprocs=None):
    """ Compute all the volumes needed for the autocorrelation

        including the volume of the initial aggregate
    """
    # prepare output array
    volumes = np.zeros(len(lradius))

    # make a contiguous copy
    spheres = spheres.copy()

    # Put work inside the queue
    njobs = (len(lradius) - 1) * nsamples + 1

    it = 0
    with Advancement(njobs) as advancement:
        for i, radius in enumerate(lradius):
            for j in range(nsamples):

                volumes[i] += single_volume_autoco(i, radius, spheres)
                it += 1
                advancement.print(it)

                # Only compute the volume of the initial aggregate once
                if i == 0:
                    break

    return volumes


def autocorrelation3d(agg, Spheres, nradius, nsamples, repart=np.geomspace,
                      start=None, end=None, nprocs=None):
    """
        Compute the normalized 3d autocorrelation graph of an aggregate
    """

    if not WithSBL:
        SBLError()
        return

    # Get the spheres of the given aggregate
    t, label = agg.name
    cols = ["Posx", "Posy", "Posz", "Radius"]
    spheres = Spheres[cols][Spheres["Label"] == label].loc[t]

    # Define range of radius
    if start is None:
        start = spheres.Radius.min() / 2
    if end is None:
        end = agg.Rmax * 2
    radius = np.array([0] + list(repart(start, end, nradius)))

    # convert pandas to numpy
    spheres = spheres.values

    # computye volumes
    if nprocs is None or nprocs > 0:
        volumes = parallel_volumes_autoco(spheres, radius, nsamples)
    else:
        volumes = seq_volumes_autoco(spheres, radius, nsamples)

    # intersection volume
    # using A∩B = A + B - A∪b
    # and normalization
    volumes = 2 - volumes / (nsamples * volumes[0])

    return radius[1:], volumes[1:]


def discretize(Agg, Spheres, resolution=128, alphagangue=0.4):
    """
        Discretize the aggregate in a grid

        alphagangue is a parameter allowing some gangue around the aggregate
    """
    t, label = Agg.name

    cols = ["Posx", "Posy", "Posz", "Radius"]

    # get the spheres of the given aggregate
    spheres = Spheres[cols][Spheres["Label"] == label].loc[t].values.copy()

    return discretize_spherelist(spheres, resolution, alphagangue)


def surrounding_box(spheres):
    """
        compute dimensions of the surrounding box
    """
    rmax = spheres[:, 3].max()
    xmin, xmax = spheres[:, 0].min() - rmax, spheres[:, 0].max() + rmax
    ymin, ymax = spheres[:, 1].min() - rmax, spheres[:, 1].max() + rmax
    zmin, zmax = spheres[:, 2].min() - rmax, spheres[:, 2].max() + rmax

    return (xmin, xmax), (ymin, ymax), (zmin, zmax)


def mkgrid(xbounds, ybounds, zbounds, resolution):
    """
        compute dimensions of the surrounding box
    """
    xmin, xmax = xbounds
    ymin, ymax = ybounds
    zmin, zmax = zbounds

    # create a grid of the requested resolution
    x, y, z = np.mgrid[xmin:xmax:resolution*1j,
                       ymin:ymax:resolution*1j,
                       zmin:zmax:resolution*1j]

    return (x, y, z)

def discretize_spherelist(spheres, resolution=128, alphagangue=0.4):
    """
        Discretize the aggregate in a grid

        alphagangue is a parameter allowing some gangue around the aggregate
    """

    xbounds, ybounds, zbounds = surrounding_box(spheres)

    (x, y, z) = mkgrid(xbounds, ybounds, zbounds, resolution)

    # fill the domain : positive value if inside at least a sphere
    if alphagangue > 0:
        data = np.zeros_like(x)
        for posx, posy, posz, radius in spheres:
            d = np.sqrt((posx - x) ** 2 + (posy - y) ** 2 + (posz - z) ** 2)
            data += 0.5 * (1 + erf(-(d / radius - 1) / alphagangue))
        # centering on 0
        data -= 0.5
    else:
        data = - np.ones_like(x)*np.inf
        for posx, posy, posz, radius in spheres:
            d = (posx - x) ** 2 + (posy - y) ** 2 + (posz - z) ** 2
            data = np.maximum(data, radius ** 2 - d)

    return data, (x, y, z)


def view_proj_agg(Agg, Spheres, resolution=128, projdir='z',
                  reduce=np.sum, alphagangue=0.4):
    """
        Plot a 2d projection of the aggregate in the direction x, y or z

        if reduce is np.sum (or np.mean) the aggregate will be "transparent"
        if reduce is np.max  the aggregate will be fully opaque

        alphagangue is a parameter allowing some gangue around the aggregate
    """
    data, grid = discretize(Agg, Spheres, resolution, alphagangue)

    view_proj_disc(data, grid, projdir, reduce)


def view_proj_disc(data, grid, projdir='z', reduce=np.sum):
    """
        Plot a 2d projection of the aggregate in the direction x, y or z

        if reduce is np.sum (or np.mean) the aggregate will be "transparent"
        if reduce is np.max  the aggregate will be fully opaque

        alphagangue is a parameter allowing some gangue around the aggregate
    """
    x, y, z = grid

    # binarization
    data = (data > 0.).astype(float)

    # start the figure
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # projection
    if projdir == 'z':
        data = - reduce(data, axis=2)
        x = x[:, :, 0]
        y = y[:, :, 0]
        ax.set_xlabel('x')
        ax.set_ylabel('y')
    if projdir == 'y':
        data = - reduce(data, axis=1)
        x = z[:, 0, :]
        y = x[:, 0, :]
        ax.set_xlabel('z')
        ax.set_ylabel('x')
    if projdir == 'x':
        data = - reduce(data, axis=0)
        x = y[0, :, :]
        y = z[0, :, :]
        ax.set_xlabel('y')
        ax.set_ylabel('z')

    ax.contourf(x, y, data, cmap='gray')
    plt.show()


def view_agg(Agg, Spheres, resolution=32, alphagangue=0.):
    """
        Plot a 3d view of the aggregate in the direction x, y or z

        if reduce is np.sum (or np.mean) the aggregate will be "transparent"
        if reduce is np.max  the aggregate will be fully opaque

        alphagangue is a parameter allowing some gangue around the aggregate
    """
    data, grid = discretize(Agg, Spheres, resolution, alphagangue)

    view_agg_disc(data, grid)


def view_agg_disc(data, grid):
    """
        Plot a 3d view of the aggregate in the direction x, y or z

        if reduce is np.sum (or np.mean) the aggregate will be "transparent"
        if reduce is np.max  the aggregate will be fully opaque

        alphagangue is a parameter allowing some gangue around the aggregate
    """
    x, y, z = grid

    # extract space steps
    dx = abs(x[0, 0, 0] - x[1, 0, 0])
    dy = abs(y[0, 0, 0] - y[0, 1, 0])
    dz = abs(z[0, 0, 0] - z[0, 0, 1])

    # compute the 3d surface of the aggregate
    verts, faces, *_ = measure.marching_cubes_lewiner(data, 0.,
                                                      spacing=(dx, dy, dz))

    verts[:, 0] += x[0, 0, 0]
    verts[:, 1] += y[0, 0, 0]
    verts[:, 2] += z[0, 0, 0]

    # Plot the result
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2])

    plt.show()


def coverages(Spheres, Aggregates, nprocs=None):
    """
        Compute the overlapping coefficient COV

        (which may happen with surface growth)
    """
    if nprocs is None:
        nprocs = len(os.sched_getaffinity(0))

    # Extracting
    data = Spheres[['Label', 'Posx', 'Posy', 'Posz', 'Radius']]
    data = data.copy()

    # Counting spheres
    Nspheres = data.index.get_level_values("Num")[-1] + 1

    # Reset index
    data.reset_index(level="Num", drop=True, inplace=True)
    data.reset_index(inplace=True)

    # New index
    data.set_index(["Time", "Label"], inplace=True)

    # Sorting
    Label_idx = data.index.get_level_values("Label").values
    idx = label_argsort(Nspheres, Label_idx, nprocs)

    # extracting
    X = data["Posx"].values[idx]
    Y = data["Posy"].values[idx]
    Z = data["Posz"].values[idx]
    R = data["Radius"].values[idx]

    # Get NP (already sorted)
    Np = Aggregates["Np"].values

    # Compute
    return coverages_cython(Np, X, Y, Z, R, nprocs)


def inertia(Agg, Spheres, resolution=128, alphagangue=0.4):
    """
        Compute gyration radius using discretisation

        alphagangue is a parameter allowing some gangue around the aggregate
    """
    data, grid = discretize(Agg, Spheres, resolution, alphagangue)

    return inertia_disc(data, grid)


def inertia_disc(data, grid):
    """
        Compute gyration radius using discretisation
    """
    X, Y, Z = grid

    coords = np.where(data > 0)

    # centering
    x = X[coords] - X[coords].mean()
    y = Y[coords] - Y[coords].mean()
    z = Z[coords] - Z[coords].mean()

    tensor = np.zeros((3, 3))
    tensor[0, 0] = (y ** 2 + z ** 2).mean()
    tensor[1, 1] = (x ** 2 + z ** 2).mean()
    tensor[2, 2] = (x ** 2 + y ** 2).mean()

    tensor[0, 1] = tensor[1, 0] = -(x * y).mean()
    tensor[0, 2] = tensor[2, 0] = -(x * z).mean()
    tensor[1, 2] = tensor[2, 1] = -(y * z).mean()

    eig = np.linalg.eigvals(tensor)

    return eig


def export_ddscat(Agg, Spheres, filename="ddscat_shape.dat", resolution=64,
                  alphagangue=0.4, Aggregates=None, type_limits=[0.]):
    """
        Ecriture du fichier shape.dat pour DDSCAT
    """
    if "cov" not in Agg:
        if Aggregates is None:
            print("""In export_ddscat
I need the cov attribute
You can compute it before using the `coverages` function
Or give me the Aggregates array so I can compute it for you""")
            return

        Aggregates["cov"] = coverages(Spheres, Aggregates)
        Agg = Aggregates.loc[Agg.name]

    data, grid = discretize(Agg, Spheres, resolution, alphagangue)
    x, y, z = grid
    nx, ny, nz = data.shape

    inertia_agg = inertia_disc(data, grid)
    anisotropy = inertia_agg.max()/inertia_agg.min()

    typedData = np.zeros_like(data, dtype=np.int)
    for limit in type_limits:
        typedData += data > limit

    coords = np.where(typedData > 0)
    NbDipoles = len(coords[0])

    if True:
        volume, _ = surface_volume_sbl(Agg, Spheres)

    else:
        # extract space steps
        dx = abs(x[0, 0, 0] - x[1, 0, 0])
        dy = abs(y[0, 0, 0] - y[0, 1, 0])
        dz = abs(z[0, 0, 0] - z[0, 0, 1])

        volume = NbDipoles * dx * dy * dz

    aeff = (3 * volume/(4 * np.pi)) ** (1 / 3)

    with open(filename, 'w') as f:
        f.write("Param alpha : {:9.2e},".format(alphagangue))
        f.write(" Cov : {:9.2e},".format(Agg["cov"]))
        f.write(" Rg (nm): {:9.2e},".format(Agg.Rg))
        f.write(" Anisotropie : {:9.2e},".format(anisotropy))
        f.write(" aeff (nm) : {:9.2e}\n".format(aeff))

        f.write("{}\t= NAT\n".format(NbDipoles))
        f.write("1.000000 0.000000 0.000000 = targer vector a1 (in TF)\n")
        f.write("0.000000 1.000000 0.000000 = targer vector a2 (in TF)\n")

        f.write("{}\t{}\t{}\t".format(1.0, 1.0, 1.0))
        f.write(" = lattice spacings (d_x,d_y,d_z)/d\n")

        f.write("{}\t{}\t{}\t".format(nx/2, ny/2, nz/2))
        f.write(" = lattice offset x0(1-3)")
        f.write(" = (x_TF,y_TF,z_TF)/d for dipole 0 0 0\n")

        f.write("JA IX IY IZ ICOMP(x,y,z)\n")
        for i, coord in enumerate(zip(*coords)):
            f.write("{}".format(i))
            f.write("\t{}".format(x[coord]))
            f.write("\t{}".format(y[coord]))
            f.write("\t{}".format(z[coord]))
            f.write("\t{}".format(typedData[coord]) * 3)
            f.write("\n")


if __name__ == '__main__':

    # The folder with all .h5 and .xmf files
    data_dir = Path("output_dir/")

    # Read all data
    Spheres, Aggregates = DLCA(data_dir).read()

    if False:
        # Get the one aggregate (the last)
        agg = Aggregates.iloc[-1]

        # Compute autocorrelation
        radius, volumes = autocorrelation3d(agg, Spheres, 32, 16)

        # plot result
        plt.figure()
        plt.loglog(radius, volumes, "-o")
        plt.show()

    if False:
        # select aggregates from the first half time period
        times = Aggregates.index.get_level_values(0)

        # Compute the overlapping coefficient for all the aggregates
        Aggregates["cov"] = coverages(Spheres, Aggregates)

        # Filter out aggregates that have only one spheres
        cov = Aggregates["cov"][Aggregates["Np"] > 1]

        # Compute the mean overlapping coefficient for each time step
        time_cov = cov.groupby(by="Time").agg((np.mean))

        # Plot the evolution of the overlapping coefficient in time.
        fig, ax = plt.subplots()
        ax.loglog(time_cov)
        ax.set_xlabel('Time (s)', fontsize=9)
        ax.set_ylabel('Overlapping coefficient', fontsize=9)
        plt.suptitle('Time evolution of the overlapping coefficient', fontsize=11)
        plt.show()

    if False:
        # Get the one aggregate (the last)
        agg = Aggregates.iloc[-1]

        volsbl, surfsbl = surface_volume_sbl(agg, Spheres)
        print("    Volume:  ", volsbl)
        print("    Surface: ", surfsbl)

