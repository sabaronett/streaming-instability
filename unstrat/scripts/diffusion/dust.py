#!/usr/bin/env python3
#==============================================================================
# dust.py
#
# Defines classes and functions for processing particle data from Athena++.
#
# Author:  Stanley A. Baronett, Chao-Chin Yang
# Created: 2021-01-20
# Last Modified: 2022-06-22
#==============================================================================
class DensityDistribution:
    """Contains methods to analyze the density distribution of dust
    particles.

    Instance Variables
        bins
            Boundaries of bins in density, including both endpoints.
            The density is normalized by the mean gas density.
        datadir
            Data directory.
        logscale
            Whether or not the bins are in logscale.
        taus
            List of dimensionless stopping times.
    """
    # Author: Chao-Chin Yang
    # Created: 2021-03-03
    # Last Modified: 2021-03-18

    def __init__(self, datadir="./data", max_workers=None, nproc=None):
        """Initiates an object.

        Keyword Arguments
            datadir
                Data directory.
            max_workders
                Maximum number of workers for concurrent processing.
            nproc
                Number of processes for method process() to use.
        """
        # Author: Chao-Chin Yang
        # Created: 2021-03-03
        # Last Modified: 2021-04-01
        import math
        from multiprocessing import set_start_method
        import numpy as np
        import os
        import PencilCode as pc

        # Get the normalization.
        par = pc.read.parameters(datadir=datadir)
        rho0 = par.rho0
        period = par.cs0 / par.omega

        # Determine the density range.
        time, slices = pc.read.slices("rhop", datadir=datadir)
        rhopmin, rhopmax = [], []
        for plane in slices.keys():
            s = getattr(slices, plane)
            rhopmin.append(s.min())
            rhopmax.append(s.max())
        ms = MultiSpecies(datadir=datadir, par=par, pvarfile="PVAR0")
        rhopmin = max(min(rhopmin), 0.375 * ms.get_rhopswarm().min()) / rho0
        rhopmax = max(rhopmax) / rho0

        # Arrange the density bins.
        self.logscale = rhopmax / rhopmin > 4
        if self.logscale:
            self.bins = np.logspace(math.log10(rhopmin), math.log10(rhopmax))
        else:
            self.bins = np.linspace(rhopmin, rhopmax)

        # Record the stopping times.
        self.taus = par.omega * np.array(ms.tausp)

        # Compute optimal number of workers for concurrent processing.
        if max_workers is None:
            max_workers = os.cpu_count()
        if nproc is None:
            nproc = min(ms.nsp, math.floor(math.sqrt(max_workers)))
            for n in range(nproc,0,-1):
                if max_workers % n == 0:
                    nproc = n
                    break
        self._nproc = nproc
        self._nts = max_workers // nproc

        # Save other repeating constant variables.
        self._nsp = ms.nsp
        self._par = par
        self._period = period
        self._rho0 = rho0
        self.datadir = datadir

        # Initiate the caches.
        self._cached = False

        # Initiate multiprocessing for process() and time_series().
        set_start_method("fork", force=True)
#-----------------------------------------------------------------------
    def process(self, varfile="var.dat"):
        """Finds the density distribution of a given snapshot.

        Returned Values
            time
                Time stamp of the snapshot, normalized by the local
                orbital period.
            pdf
                Probability density of each density bin; pdf[i] is the
                probability density of bins[i].
            pdfc
                Contribution to the probability density by each dust
                species; pdfc[j][i] is the part of pdf[i] contributed to
                by species j.

        Keyword Arguments
            varfile
                File name of the snapshot.
        """
        # Author: Chao-Chin Yang
        # Created: 2021-03-03
        # Last Modified: 2021-03-20
        from concurrent import futures
        import numpy as np
        from numpy import ma
        import PencilCode as pc

        # Find the distribution of the total density.
        f = pc.read.var(datadir=self.datadir, varfile=varfile, verbose=False)
        rhop0 = f.rhop / self._rho0
        pdf, bins = np.histogram(rhop0, bins=self.bins, density=True)
        counts, bins = np.histogram(rhop0, bins=self.bins)

        # Feed constant variables and (small) arrays to _one_species().
        self._pvarfile = ('p' if varfile[0] == 'v' else 'P') + varfile
        self._rhop = ma.array(f.rhop, mask = (f.rhop <= 0))
        self._counts = ma.array(counts, mask = (counts <= 0))

        # Process each dust species.
        nbins = len(self.bins) - 1
        frac = np.empty((self._nsp, nbins))
        with futures.ProcessPoolExecutor(max_workers=self._nproc) as executor:
            f2i = {}
            for i in range(self._nsp):
                f2i[executor.submit(self._one_species, i)] = i
            for future in futures.as_completed(f2i):
                frac[f2i[future]] = future.result()

        return f.t / self._period, pdf, pdf * frac
#-----------------------------------------------------------------------
    def time_average(self, taus_bins=None, tmin=None):
        """Finds the time average and variability of the distribution.

        Keyword Arguments
            taus_bins
                Bin boundaries in dimensionless stopping times,
                excluding the two ends.
            tmin
                Minimum time used for the average; if None, read from
                data file "turbulence.json"; 0 if the file does not
                exist.

        Returned Values
            (avg, var)
                Time average and variability of the probability density,
                respectively; avg[i] and var[i] are the average and
                variability of density bins[i].
            (avgc, varc)
                Part of (avg, var) contributed to by each dust species;
                avgc[j] and varc[j] are the part of (avg, var)
                contributed to by species j.
        """
        # Author: Chao-Chin Yang
        # Created: 2021-03-06
        # Last Modified: 2021-03-10
        import json
        import numpy as np
        from pathlib import Path

        # Get the time series.
        times, pdfs, pdfcs = self.time_series(taus_bins=taus_bins)

        if tmin is None:
            # Look up the saturation time, if any.
            p = Path(self.datadir) / "turbulence.json"
            if p.exists():
                with open(p) as f:
                    x = json.load(f)
                tmin = x["tsat"]
            else:
                print("Data file '{}' does not exist. ".format(p.as_posix()))
                tmin = 0
            print("Set tmin = {}. ".format(tmin))

        # Define function to compute average and variability.
        def stats(t, v):
            duration = t[-1] - t[0]
            avg = np.trapz(v, t, axis=0) / duration
            var = np.sqrt(np.trapz((v - avg)**2, t, axis=0) / duration)
            return avg, var

        # Cut off the series below tmin.
        indices = np.nonzero(times >= tmin)
        t = times[indices]
        p = pdfs[indices]
        q = pdfcs.swapaxes(0,1)[indices]

        # Process the PDF.
        return stats(t, p), stats(t, q)
#-----------------------------------------------------------------------
    def time_series(self, taus_bins=None):
        """Finds the evolution of the density distribution.

        Keyword Arguments
            taus_bins
                Bin boundaries in dimensionless stopping times,
                excluding the two ends.

        Returned Values
            times
                Time series, normalized by the local orbital period.
            pdfs
                Probability density of each density bin at each instant
                of time; pdfs[j][i] is the probability density of
                bins[i] at times[j].
            pdfcs
                Contribution to the probability density by each dust
                species at each instant of time; pdfcs[j][k][i] is the
                part of pdfs[j][i] contributed to by species k.
        """
        # Author: Chao-Chin Yang
        # Created: 2021-03-05
        # Last Modified: 2021-03-16
        from concurrent import futures
        import numpy as np
        from pathlib import Path

        # Use the caches if they are filled.
        p = Path(self.datadir)
        flist = p / "allprocs" / "varN.list"
        dfile = p / "rhop-pdf.npz"
        if not self._cached and dfile.exists():
            if dfile.stat().st_mtime > flist.stat().st_mtime:
                print("Reading from data file {}......".format(dfile))
                a = np.load(dfile)
                self._times = a["times"]
                self._pdfs = a["pdfs"]
                self._pdfcs = a["pdfcs"]
                self._cached = True
            else:
                print("Data file '{}' is outdated. ".format(dfile))
        if self._cached:
            pdfcs = self._bin_taus(self._pdfcs, bins=taus_bins)
            return self._times, self._pdfs, pdfcs

        # Get the list of snapshots.
        files = []
        with open(flist) as f:
            for line in f:
                varfile, time = line.split()
                files.append(varfile)
        nfiles = len(files)

        # Assemble the time series.
        progress = "\rProcessing ({:4.0%})......"
        print(progress.format(0), end='', flush=True)
        nbins = len(self.bins) - 1
        times = np.empty((nfiles,))
        pdfs = np.empty((nfiles,nbins))
        pdfcs = np.empty((nfiles,self._nsp,nbins))
        with futures.ProcessPoolExecutor(max_workers=self._nts) as executor:
            f2i = {}
            for i, varfile in enumerate(files):
                f2i[executor.submit(self.process, varfile=varfile)] = i
            for n, f in enumerate(futures.as_completed(f2i)):
                print(progress.format((n + 1) / nfiles), end='', flush=True)
                time, pdf, pdfc = f.result()
                i = f2i[f]
                times[i] = time
                pdfs[i] = pdf
                pdfcs[i] = pdfc
        pdfcs = pdfcs.swapaxes(0,1)
        print("Done. ")

        # Save, cache and return the results.
        np.savez(dfile, times=times, pdfs=pdfs, pdfcs=pdfcs)
        print("Updated data file {}. ".format(dfile))
        self._times = times
        self._pdfs = pdfs
        self._pdfcs = pdfcs
        self._cached = True
        pdfcs = self._bin_taus(self._pdfcs, bins=taus_bins)
        return times, pdfs, pdfcs
#-----------------------------------------------------------------------
    def _bin_taus(self, a, bins=None):
        """Sums an array by bins of dust species.

        Positional Arguments
            a
                An array, with the first dimension being the species.

        Keyword Arguments
            bins
                Bin boundaries in dimensionless stopping times,
                excluding the two ends.

        Returned Value
            A list of the summed arrays by bins.
        """
        # Author: Chao-Chin Yang
        # Created: 2021-03-09
        # Last Modified: 2021-03-10
        import numpy as np

        # Do nothing if no bins are specified.
        if bins is None:
            return a

        # Define function to select species in a bin.
        all_true = np.full_like(self.taus, True, dtype=bool)
        def masked_sum(tausmin=None, tausmax=None):
            u = all_true if tausmin is None else (self.taus >= tausmin)
            v = all_true if tausmax is None else (self.taus < tausmax)
            return np.expand_dims(a[u & v].sum(axis=0), 0)

        # Process each bin.
        s = masked_sum(tausmax=bins[0])
        for t1, t2 in zip(bins[:-1], bins[1:]):
            s = np.vstack((s, masked_sum(t1, t2)))
        s = np.vstack((s, masked_sum(tausmin=bins[-1])))

        return s
#-----------------------------------------------------------------------
    def _one_species(self, s):
        """Helps method process() as a child process to process one dust
        species.

        Positional Arguments
            s
                Index of the dust species.

        Returned Values
            Fractional contribution of the species to the probability
            density function.
        """
        # Author: Chao-Chin Yang
        # Created: 2021-03-11
        # Last Modified: 2021-03-20
        import numpy as np

        # Get the particle data.
        ms = MultiSpecies(datadir=self.datadir, pvarfile=self._pvarfile,
                par=self._par)
        fp = ms.get_species(s)

        # Find the fractional contribution in each cell.
        pm = ParticleMesh(datadir=self.datadir, par=self._par)
        weights = pm.assign(fp.rhopswarm, fp.xp, fp.yp, fp.zp)
        weights = weights / self._rhop

        # Compute the histogram.
        rhop0 = self._rhop / self._rho0
        hist, bins = np.histogram(rhop0, bins=self.bins, weights=weights)
        return hist / self._counts
#==============================================================================
class Diffusion:
    """Contains methods to analyze the diffusion of dust particles.

    Instance Variables
        time
            Time series, normalized by local orbital period.
        dxp
            A 2D numpy array, where dxp[i,j] is the radial component of the
            displacement of particle j at time[i], normalized by the local
            scale height of the gas.  It is None if this dimension is neither
            active nor periodic.
        dxp_mean, dxp_std
            2D numpy arrays, where dxp_mean[k,i] and dxp_std[k,i] is the mean
            and standard deviation of the radial displacement for species k at
            time[i], respectively.  They are None if this dimension is neither
            active nor periodic.
        dyp, dyp_mean, dyp_std
            Similar to dxp, dxp_mean and dxp_std, except for the azimuthal
            component.
        dzp, dzp_mean, dzp_std
            Similar to dzp, dxp_mean and dxp_std, except for the vertical
            component.
    """
    # Author: Stanley A. Baronett, Chao-Chin Yang
    # Created: 2021-05-07
    # Last Modified: 2022-06-28

    def __init__(self, datadir="./dat", athinput=None, nout=3):
        """Finds the displacement of each dust particle as a function of
        time.

        Keyword Arugments
            datadir
                Data directory.
            athinput
                If not None, a dictionary of dictionaries for athinput file.
            nout
                Number of particle outputs (.dat) to include, ending with the
                last available (highest numbered) file.
        """
        # Author: Stanley A. Baronett, Chao-Chin Yang
        # Created: 2021-05-07
        # Last Modified: 2022-06-28
        import sys
        sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
        import athena_read
        import math
        import numpy as np
        # import PencilCode as pc
        import find
        # import read

        # Find the displacement.
        time, (dxp, dyp, dzp) = find.par_disp(datadir=datadir, save_to='disp',
                                              athinput=athinput, nout=nout)

        # Normalize the time and displacement.
        if athinput is None:
            athinput = athena_read.athinput('athinput.si')
        omega = athinput['problem']['omega']
        time /= 2 * math.pi / omega
        hgas = athinput['hydro']['iso_sound_speed'] / omega
        if dxp is not None: dxp /= hgas
        if dyp is not None: dyp /= hgas
        if dzp is not None: dzp /= hgas

        # Save the results.
        self.time = time
        self.dxp, self.dyp, self.dzp = dxp, dyp, dzp

        # Define function to process one dimension [for each species.]
        # ms = MultiSpecies(datadir=datadir, par=par)
        def process_1d(dxp):
            if dxp is None:
                return None, None

            # Find the mean and standard deviation of the displacement.
            mean, std = [], []
            # for s in range(ms.nsp):
            # dxp1 = dxp[:,ms.get_indices_of_species(s)]
            mean.append(dxp.mean(axis=1))
            std.append(dxp.std(axis=1))

            return np.vstack(mean), np.vstack(std)

        # Process each dimension.
        self.dxp_mean, self.dxp_std = process_1d(dxp)
        self.dyp_mean, self.dyp_std = process_1d(dyp)
        self.dzp_mean, self.dzp_std = process_1d(dzp)
#-----------------------------------------------------------------------
    def coefficient(self):
        """Finds the diffusion coefficient of each dust species.

        Returned Values:
            dpx
                Diffusion coefficient for the radial dimension, where
                dpx[i] is the coefficient of species i, in terms of
                Hg^2/P with Hg being gas scale height and P local
                orbital period.  It is None if the dimension is neither
                active nor periodic.
            dpy
                Similar to dpx, except for the azimuthal dimension.
            dpz
                Similar to dpx, except for the vertical dimension.
        """
        # Author: Chao-Chin Yang
        # Created: 2021-06-04
        # Last Modified: 2021-06-04
        import numpy as np

        # Define function to compute the diffusion coefficient.
        def compute(time, disp):
            if disp is None:
                return None

            tsqrt = np.sqrt(time)
            tsum = np.sum(time)
            sigma = []
            for d in disp:
                sigma.append(np.sum(d * tsqrt) / tsum)
            return 0.5 * np.array(sigma)**2

        # Compute and return the coefficient.
        dpx = compute(self.time, self.dxp_std)
        dpy = compute(self.time, self.dyp_std)
        dpz = compute(self.time, self.dzp_std)
        return dpx, dpy, dpz
#=======================================================================
class Kinematics:
    """Contains the analysis of the kinematics of each dust species
    from a given Pencil Code model.

    Instance Variables:
        nsp
            Number of dust species.
        taus
            Dimensionless stopping time of each species.
        epsilon
            Solid-to-gas density ratio of each species.
        time
            Time series, normalized by orbital period.
        vpxm, dvpx
            vpxm[i][j] and dvpx[i][j] are the mean and the standard
            deviation of the radial velocity of species i at time[j],
            respectively, normalized by the speed of sound.
        vpym, dvpy
            Similar to vpxm and dvpx, except for the azimuthal
            direction.
        vpzm, dvpz
            Similar to vpxm and dvpx, except for the vertical direction.
        vpx0, vpy0
            vpx0[i] and vpy0[i] are the radial and azimuthal components
            of the initial equilibrium velocity of species i, normalized
            by the speed of sound.
        ux0, uy0
            Radial and azimuthal components of the initial equilibrium
            velocity of gas, normalized by the speed of sound.
    """
    # Author: Chao-Chin Yang
    # Created: 2020-12-17
    # Last Modified: 2020-04-25

    def __init__(self, datadir="./data"):
        """Reads and processes the particle data and saves the results.

        Keyword Arugments
            datadir
                Data directory.
        """
        # Author: Chao-Chin Yang
        # Created: 2020-12-17
        # Last Modified: 2021-04-25
        import math
        import numpy as np
        import PencilCode as pc
        from streaming_instability import SIMode

        # Get a list of PVAR files.
        pvarlist = pc.read.varlist(datadir=datadir, listname="pvarN.list")[0]

        # Process each PVAR file.
        rhopswarm = None
        time, mean, std = [], [], []
        par = pc.read.parameters(datadir=datadir)
        for pvarfile in pvarlist:
            print("Processing {}......".format(pvarfile))

            # Read the particle data.
            ms = MultiSpecies(datadir=datadir, par=par, pvarfile=pvarfile)
            time.append(ms.fp.t)

            # Record and check the mass of each species.
            if rhopswarm is None:
                rhopswarm = ms.get_rhopswarm()
            elif any(rhopswarm != ms.get_rhopswarm()):
                raise RuntimeError("inconsistent rhopswarm over time")

            # Find mean velocity and velocity dispersion.
            for i in range(ms.nsp):
                fp = ms.get_species(i)
                for v in (fp.vpx, fp.vpy, fp.vpz):
                    mean.append(v.mean())
                    std.append(v.std())

        # Convert rhopswarm to epsilon.
        rhopswarm *= ms.nppsc / par.rho0

        # Vectorize the data and normalize.
        time = np.array(time) * (par.omega / (2 * math.pi))
        shape = len(time), ms.nsp, 3
        mean = (np.array(mean) / par.cs0).reshape(shape).transpose()
        std = (np.array(std) / par.cs0).reshape(shape).transpose()

        # Record the results.
        self.nsp = ms.nsp
        self.taus = par.omega * np.array(ms.tausp)
        self.epsilon = rhopswarm
        self.time = time
        self.vpxm, self.vpym, self.vpzm = mean
        self.dvpx, self.dvpy, self.dvpz = std

        # Find the equilibrium velocities.
        self.pi = -0.5 * par.beta_glnrho_global[0]
        eq = SIMode.equilibrium(self.pi, self.taus, self.epsilon)
        self.vpx0 = eq.vpx
        self.vpy0 = eq.vpy
        self.ux0 = eq.ux
        self.uy0 = eq.uy
#-----------------------------------------------------------------------
    def time_average(self, fpath=None, tmin=None):
        """Finds the time average of the mean and standard deviation of
        the velocity for each species.

        Keyword Arguments:
            fpath
                If not None, the path to a file where the results and
                some other data are saved.
            tmin
                Minimum time to begin the average, if not None.

        Returned Values:
            (vpxm, vpym, vpzm)
                Components of the mean velocity of each species,
                normalized by the speed of sound.
            (dvpx, dvpy, dvpz)
                Components of the velocity dispersion of each species,
                normalized by the speed of sound.
        """
        # Author: Chao-Chin Yang
        # Created: 2020-12-28
        # Last Modified: 2021-04-25
        import numpy as np
        from utilities import time_average

        # Integrate.
        mean, std = [], []
        for v, dv in zip((self.vpxm, self.vpym, self.vpzm),
                         (self.dvpx, self.dvpy, self.dvpz)):
            for vi, dvi in zip(v, dv):
                mean.append(time_average(self.time, vi, tmin=tmin))
                std.append(time_average(self.time, dvi, tmin=tmin))

        # Vectorize the results.
        vpxm, vpym, vpzm = np.array(mean).reshape(3,self.nsp)
        dvpx, dvpy, dvpz = np.array(std).reshape(3,self.nsp)

        if fpath is not None:
            # Save the results.
            if tmin is None: tmin = self.time[0]
            np.savez(fpath, tmin=tmin,
                    pi=self.pi, taus=self.taus, epsilon=self.epsilon,
                    vpx0=self.vpx0, vpy0=self.vpy0, ux0=self.ux0, uy0=self.uy0,
                    vpxm=vpxm, dvpx=dvpx,
                    vpym=vpym, dvpy=dvpy,
                    vpzm=vpzm, dvpz=dvpz)
            print("Saved the time averages in file " + fpath + ".npz. ")

        return (vpxm, vpym, vpzm), (dvpx, dvpy, dvpz)
#=======================================================================
class MultiSpecies:
    """A class that holds a snapshot of particle data, which can be
    decomposed into individual dust species.

    Instance Variables:
        fp
            Snapshot of the full particle data.
        npps
            Number of particles per species.
        nppsc
            Number of particles per species per cell.
        nsp
            Number of dust species.
        tausp
            List of (physical) stopping times in code units.
    """
    # Author: Chao-Chin Yang
    # Created: 2021-01-20
    # Last Modified: 2021-02-19

    def __init__(self, datadir="./data", par=None, pvarfile="pvar.dat"):
        """Initiates an object, holding a snapshot of particle data.

        Keyword Arguments
            datadir
                Data directory.
            par
                If not None, init parameters provided by module
                PencilCode.
            pvarfile
                Name of the particle data file.
        """
        # Author: Chao-Chin Yang
        # Created: 2021-01-20
        # Last Modified: 2021-02-19
        import PencilCode as pc

        # Get the number of species.
        if par is None:
            par = pc.read.parameters(datadir=datadir)
        if type(par.tausp_species) is float:
            self.nsp = 1
        else:
            self.nsp = len(par.tausp_species)

        # Get the stopping times.
        self.tausp = par.tausp_species if self.nsp > 1 else (par.tausp,)

        # Find the number of particles per species.
        pdim = pc.read.pardim(datadir=datadir)
        self.npps = pdim.npar // self.nsp

        # Find the number of particles per species per cell.
        dim = pc.read.dimensions(datadir=datadir)
        self.nppsc = self.npps // (dim.nxgrid * dim.nygrid * dim.nzgrid)

        # Read the snapshot.
        self.fp = pc.read.pvar(datadir=datadir, pvarfile=pvarfile)
#-----------------------------------------------------------------------
    def get_indices_of_species(self, species):
        """Lists the (ipar) indices of the particles of a given species.

        Positional Arguments
            species
                Index of the species (starting from 0).

        Returned Values:
            The (ipar) indices.
        """
        # Author: Chao-Chin Yang
        # Created: 2021-05-07
        # Last Modified: 2021-05-07

        return range(species * self.npps, (species + 1) * self.npps)
#-----------------------------------------------------------------------
    def get_rhopswarm(self):
        """Finds the mass (rhopswarm) of each dust species.

        Returned Values:
            A list of rhopswarm of individual species as a numpy array.
        """
        # Author: Chao-Chin Yang
        # Created: 2020-12-17
        # Last Modified: 2021-05-07
        import numpy as np

        rhopswarm = []
        for i in range(self.nsp):
            rhop = self.fp.rhopswarm[self.get_indices_of_species(i)]
            rhop0 = rhop.min()
            if rhop0 != rhop.max():
                raise RuntimeError("inconsistent rhopswarm for " +
                        "species {}: ".format(i+1) +
                        "min, max = {}, {}".format(rhop0, rhop.max()))
            rhopswarm.append(rhop0)

        return np.array(rhopswarm)
#-----------------------------------------------------------------------
    def get_species(self, species):
        """Gets the data of an individual species.

        Positional Arguments
            species
                Index of the species (starting from 0).

        Returned Values:
            Named tuple of the particle data of the specified species.
        """
        # Author: Chao-Chin Yang
        # Created: 2020-01-20
        # Last Modified: 2021-09-29
        from PencilCode.read import Dict

        # Sanity check.
        if type(species) is not int:
            raise TypeError("Argument species is not an integer. ")
        elif species < 0 or species >= self.nsp:
            raise IndexError(f"species = {species} is out of range. ")

        # Extract the subset of particle data.
        npar = self.npps
        pvar = dict(npar=npar, t=self.fp.t)
        indices = self.get_indices_of_species(species)
        for v in self.fp.keys():
            if v == "npar" or v == "t": continue
            pvar[v] = getattr(self.fp, v)[indices]
        return Dict(**pvar)
#=======================================================================
class ParticleMesh:
    """A class that handles particle-mesh operations. """
    # Author: Chao-Chin Yang
    # Created: 2021-01-25
    # Last Modified: 2021-03-17

    def __init__(self, datadir="./data", par=None):
        """Initiates an object, registering mesh info.

        Keyword Arguments
            datadir
                Data directory.
            par
                If not None, init parameters provided by module
                PencilCode.
        """
        # Author: Chao-Chin Yang
        # Created: 2021-01-25
        # Last Modified: 2021-03-17
        import numpy as np
        import PencilCode as pc

        # Read the dimensions.
        dim = pc.read.dimensions(datadir=datadir)
        nxyz = dim.nxgrid, dim.nygrid, dim.nzgrid
        active_dim = np.array(nxyz) > 1
        self._nghost = 1

        # Compute the dimensions of the mesh array.
        if par is None: par = pc.read.parameters(datadir=datadir)
        mshape, wshape, xyz0, dxyz = [], [], [], []
        for i, active in enumerate(active_dim):
            if active:
                mshape.append(nxyz[i] + 2 * self._nghost)
                wshape.append(1 + 2 * self._nghost)
                xyz0.append(par.xyz0[i])
                dxyz.append(par.lxyz[i] / nxyz[i])
        self._active_dim = np.nonzero(active_dim)
        self._mshape = tuple(mshape)
        self._wshape = tuple(wshape)
        self._xyz0 = tuple(xyz0)
        self._dxyz = tuple(dxyz)

        # Sanity check.
        if not any(par.lequidist):
            raise NotImplemented("not working with non-equidistant grid. ")
#-----------------------------------------------------------------------
    def assign(self, vp, xp, yp, zp):
        """Assign a property of the particles to the grid.

        Positional Arguments:
            vp
                The property to be assigned.
            xp, yp, zp
                Position of each particle.

        Returned Values:
            The property of the particles on the grid.
        """
        # Author: Chao-Chin Yang
        # Created: 2021-01-25
        # Last Modified: 2021-01-26
        import numpy as np

        # Define function to convert a coordinate to the index space.
        def index(x, x0, dx):
            return self._nghost + (x - x0) / dx

        # Define function to find the weight.
        @np.vectorize
        def weight(dxi):
            x = abs(dxi)
            if x < 0.5:
                return 0.75 - x**2
            elif x < 1.5:
                return 0.5 * (1.5 - x)**2
            else:
                return 0

        # Arrange the coordinates.
        rp = np.vstack((xp, yp, zp))[self._active_dim].transpose()

        # Assign each particle to the mesh.
        vg = np.zeros(self._mshape)
        nghost = self._nghost
        for v, r in zip(vp, rp):
            w, s = np.array([1]), []
            for i, x in enumerate(r):
                xi = index(x, self._xyz0[i], self._dxyz[i])
                arg = int(xi)
                arg = arg - nghost, arg + nghost + 1
                w = np.outer(w, weight(xi - 0.5 - np.arange(*arg)))
                s.append(slice(*arg))
            vg[tuple(s)] += v * w.reshape(self._wshape)

        # Fold the boundary values, assuming periodic.
        for i in range(len(self._mshape)):
            vg = np.moveaxis(vg, 0, -1)
            vg[nghost:2*nghost] += vg[-nghost:]
            vg[-2*nghost:-nghost] += vg[:nghost]
            vg = vg[nghost:-nghost]

        return vg
