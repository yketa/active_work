import numpy as np
from collections import OrderedDict

def squared_displacement(dat, frame, dt):
    """
    Returns squared displacement between two frames.
    """

    displacements = (dat.getPositions(frame + dt, wrapped=False)
        - dat.getPositions(frame, wrapped=False))

    if dat.N > 1:
        sqdisp = np.sum((displacements - np.mean(displacements, axis=0))**2,
            axis=-1)
    else:
        sqdisp = np.sum(displacements**2, axis=-1)

    return sqdisp

def msd(dat, dt_max=1000, n_max=1000):
    """
    Parameters
    ----------
    dat : Dat
        Data class.
    dt_max : int
        Maximum number of computed lag times.
    n_max : int
        Maximum number of intervals for computation of mean square displacement.

    Returns
    -------
    msdstderr : (*, 3) Numpy array
        (lag time, stderr, mean square displacement)

    NOTE: ASSUMES TIME STEPS ARE ALWAYS THE SAME.
    """

    timeStep = dat.getTimeStep(1)
    msdstderr = []

    if dt_max > dat.frames - 1:
        lag_times = range(1, dat.frames)
    else:
        lag_times = list(OrderedDict.fromkeys(map(
            int,
            np.exp(np.linspace(np.log(1), np.log(dat.frames - 1), dt_max)))))

    for dt in lag_times:

        lag_time = timeStep*dt

        frames = list(OrderedDict.fromkeys(
            np.linspace(0, dat.frames - dt - 1,
                min(n_max, dat.frames - dt), dtype=int)))
        sqdisp = np.array(list(map(
            lambda frame: squared_displacement(dat, frame, dt),
            frames))).flatten()

        msd, stderr = (np.mean(sqdisp),
            np.std(sqdisp)/np.sqrt(np.prod(sqdisp.shape)))

        msdstderr += [[lag_time, msd, stderr]]

    return np.array(msdstderr)

def msd_th(dt, lp):
    """
    Returns theoretical mean square displacement for ABP at dimensionless time
    according to dimensionless persistence length.
    """

    return 4*(1/(3*lp) + lp/2)*(dt + lp*(np.exp(-dt/lp) - 1))
