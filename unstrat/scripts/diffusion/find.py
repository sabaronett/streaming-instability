#!/usr/bin/env python3
#==============================================================================
# find.py
#
# Facilities for analyzing Athena++ data.
#
# Author: Stanley A. Baronett, Chao-Chin Yang
# Created: 2013-10-21
# Last Modified: 2022-06-28
#==============================================================================
def par_disp(datadir='./dat', athinput=None, save_to=None):
    """Finds the displacement of each particle as a function of time.

    Keyword Arguments:
        datadir
            Path to the data directory.
        athinput
            If not None, a dictionary of dictionaries for athinput file.
        save_to
            If not None, a string of the filename (without .npz
            extension) to save the results (t and dxp) to a numpy data
            file under datadir.

    Returned Values:
        t
            A numpy array of times, starting from zero.
        dxp, dyp, dzp
            A tuple of three 2D numpy arrays, where dxp[i,j], dyp[i,j],
            and dzp[i,j] are the components of the displacement of
            particle j at t[i].  A component is returned None if the
            dimension is neither active nor periodic.
    """
    # Author: Stanley A. Baronett, Chao-Chin Yang
    # Created: 2018-01-16
    # Last Modified: 2022-06-28
    import sys
    sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
    import athena_read
    import numpy as np
    from pathlib import Path

    # Retrieve athinput file
    if athinput is None:
        athinput = athena_read.athinput('athinput.si')

    # Determine numbers of particles and snapshots.
    outputs = sorted(list(Path(datadir).glob(athinput['job']['problem_id'] +
                                             '.pout.*.dat')))
    del outputs[0]  # exclude t=0 snapshot
    time, pdata = athena_read.particles(str(outputs[0]))
    npar = pdata.size
    nt = len(outputs)
    print(f'Total number of particles: {npar}')
    print(f'Total number of snapshots: {nt}')

    # Determine which dimension(s) should be considered.
    mesh = athinput['mesh']
    lx = mesh['x1max'] - mesh['x1min']
    ly = mesh['x3max'] - mesh['x3min']
    lz = mesh['x2max'] - mesh['x2min']
    active = (mesh[f'nx1'] > 1 and mesh[f'ix1_bc'] == 'periodic'\
                and mesh[f'ox1_bc'] == 'periodic',
              mesh[f'nx2'] > 1 and mesh[f'ix2_bc'] == 'periodic'\
                and mesh[f'ox2_bc'] == 'periodic',
              mesh[f'nx3'] > 1 and mesh[f'ix3_bc'] == 'periodic'\
                and mesh[f'ox3_bc'] == 'periodic')
    print(f'Active dimensions: {active}')

    # Define function to detect boundary jumping.
    @np.vectorize
    def jump(xp1, xp2, lx):
        if 2 * abs(xp2 - xp1) > lx:
            return -1 if xp2 > xp1 else +1
        else:
            return 0

    # Define function to record positions.
    def record(xp, xp_old, xshift, lx):
        xshift += jump(xp_old, xp, lx)
        return xp + xshift * lx, xp, xshift

    # Allocate memory.
    t = np.empty(nt,)
    def alloc(switch):
        dxp, xshift = None, None
        if switch:
            dxp = np.empty((nt,npar))
            xshift = np.zeros((npar,), dtype=int)
        return dxp, xshift
    dxp, xshift = alloc(active[0])
    dyp, yshift = alloc(active[1])
    dzp, zshift = alloc(active[2])
    vpxmin, vpxmax = float('inf'), -float('inf')
    vpymin, vpymax = float('inf'), -float('inf')
    vpzmin, vpzmax = float('inf'), -float('inf')

    # Process each snapshot of particles.
    for i, output in enumerate(outputs):
        print('\rProcessing pout files ({:6.1%})......'.format((i+1)/nt),
              end='', flush=True)
        time, pdata = athena_read.particles(str(output))

        # Monitor velocity extrema.
        vpxmin = min(vpxmin, min(pdata['vpx']))
        vpxmax = max(vpxmax, max(pdata['vpx']))
        vpymin = min(vpymin, min(pdata['vpy']))
        vpymax = max(vpymax, max(pdata['vpy']))
        vpzmin = min(vpzmin, min(pdata['vpz']))
        vpzmax = max(vpzmax, max(pdata['vpz']))

        # Record time and positions.
        t[i] = time
        if i == 0:
            if active[0]: dxp[0] = xp_old = pdata['xp']
            if active[1]: dyp[0] = yp_old = pdata['yp']
            if active[2]: dzp[0] = zp_old = pdata['zp']
        else:
            if active[0]:
                dxp[i], xp_old, xshift = record(pdata['xp'], xp_old, xshift,
                                                lx)
            if active[1]:
                dyp[i], yp_old, yshift = record(pdata['yp'], yp_old, yshift,
                                                ly)
            if active[2]:
                dzp[i], zp_old, zshift = record(pdata['zp'], zp_old, zshift,
                                                lz)

    print('Done. ')

    # Check the dimensions.
    dtmax = max(t[1:] - t[:-1])
    if (active[0] and 2 * max(abs(vpxmax), abs(vpxmin)) * dtmax > lx or
        active[1] and 2 * max(abs(vpymax), abs(vpymin)) * dtmax > ly or
        active[2] and 2 * max(abs(vpzmax), abs(vpzmin)) * dtmax > lz):
        print(f'dtmax = {dtmax}')
        if active[0]: print(f'vpxmin, vpxmax = {vpxmin}, {vpxmax}')
        if active[1]: print(f'vpymin, vpymax = {vpymin}, {vpymax}')
        if active[2]: print(f'vpzmin, vpzmax = {vpzmin}, {vpzmax}')
        print('Warning: Boundary jumping may not be properly detected. ')

    # Remove duplicate data.
    print('Removing duplicate data......', end='', flush=True)
    t, indices = np.unique(t, return_index=True)
    if active[0]: dxp = dxp[indices]
    if active[1]: dyp = dyp[indices]
    if active[2]: dzp = dzp[indices]
    nt = len(t)
    print('Done. ')

    # Find the displacement.
    print('Computing the displacement......', end='', flush=True)
    t -= t[0]
    if active[0]: dxp -= dxp[0]
    if active[1]: dyp -= dyp[0]
    if active[2]: dzp -= dzp[0]
    print('Done. ')

    # Save the results if requested.
    if save_to is not None:
        file = f'{datadir}/' + save_to
        print('Saving the results to ' + file + '.npz......',
              end='', flush=True)
        kw = dict(t=t)
        if active[0]: kw['dxp'] = dxp
        if active[1]: kw['dyp'] = dyp
        if active[2]: kw['dzp'] = dzp
        np.savez(file, **kw)
        print('Done. ')

    return t, (dxp, dyp, dzp)
