#!/usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import healpy as hp
from rubin_sim.selfcal import generate_catalog, OffsetSys, OffsetSNR, LsqrSolver
import sqlite3
import os
from rubin_sim.utils import healbin
import argparse


def run_selfcal(pointings_db, stars_db="msrgb_1e6.sqlite", filtername="r", night_max=365.25,
                error_sys=0.03, n_patches=16):
    """Run a selfcalibration simulation
    """

    con = sqlite3.connect(pointings_db)
    # read with pandas
    visits_df = pd.read_sql("select observationId, fieldra as ra,fielddec as "
                            "dec,fiveSigmaDepth,rotSkyPos "
                            "from observations where filter='%s' and night < %f;" % (filtername, night_max), con)
    con.close()

    # Dump the dataframe into a simple numpy array. The to_numpy() method doesn't handle mixed dtypes well
    names = ['observationId', 'ra', 'dec', 'fiveSigmaDepth', 'rotSkyPos']
    types = [int, float, float, float, float]
    visits = np.empty(visits_df.shape[0], dtype=list(zip(names, types)))
    for key in names:
        visits[key] = visits_df[key].values

    # Read the stars in
    con = sqlite3.connect(stars_db)
    stars_df = pd.read_sql("select id,ra,decl,%smag from stars;" % filtername, con)
    con.close()
    # Repackage as numpy array
    names = ['id', 'ra', 'decl', '%smag' % filtername]
    types = [int, float, float, float]
    stars = np.empty(stars_df.shape[0], dtype=list(zip(names, types)))
    for key in names:
        stars[key] = stars_df[key].values

    offsets = [OffsetSys(error_sys=error_sys), OffsetSNR(lsst_filter=filtername)]

    observed_stars = generate_catalog(visits, stars, offsets=offsets,
                                      lsst_filter=filtername, n_patches=n_patches,
                                      verbose=False)

    solver = LsqrSolver(observed_stars, patch_out=None, star_out=None)

    solver.run()
    fit_patches, fit_stars = solver.return_solution()

    return fit_patches, fit_stars


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--filtername", type=str, default="r")
    parser.add_argument("--night_max", type=float, default=365.25)
    parser.add_argument("--n_patches", type=int, default=16)
    parser.add_argument("--pointings_db", type=str, default='')

    args = parser.parse_args()

    run_name = os.path.basename(args.pointings_db).replace(".db", "")
    years = args.night_max/365.25

    out_name = '%s_%s_%.1f_%i.npz' % (run_name, args.filtername, years, args.n_patches)

    fit_patches, fit_stars = run_selfcal(args.pointings_db, filtername=args.filtername,
                                         night_max=args.night_max, n_patches=args.n_patches)

    np.savez(out_name, fit_patches=fit_patches, fit_stars=fit_stars)
