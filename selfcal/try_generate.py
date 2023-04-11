import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from rubin_sim.selfcal import generate_catalog, OffsetSys, OffsetSNR, LsqrSolver
import sqlite3
from rubin_sim.data import get_baseline


if __name__ == "__main__":

    filter_name = "u"

    # Load up pointings
    baseline_survey = get_baseline()
    baseline_survey
    con = sqlite3.connect(baseline_survey)
    visits_df = pd.read_sql("select observationId, fieldra as ra,fielddec as "
                            "dec,fiveSigmaDepth,rotSkyPos from observations where filter='%s' and night < 365;" % filter_name, con)
    con.close()
    # what is the right way to do this?
    names = ['observationId', 'ra', 'dec', 'fiveSigmaDepth', 'rotSkyPos']
    types = [int, float, float, float, float]
    visits = np.empty(visits_df.shape[0], dtype=list(zip(names, types)))
    for key in names:
        visits[key] = visits_df[key].values

    # Load up stars
    con = sqlite3.connect("msrgb_1e6.sqlite")
    stars_df = pd.read_sql("select id,ra,decl,%smag from stars;" % filter_name, con)
    con.close()
    # what is the right way to do this? There's gotta be a better way!
    names = ['id', 'ra', 'decl', '%smag' % filter_name]
    types = [int, float, float, float]
    stars = np.empty(stars_df.shape[0], dtype=list(zip(names, types)))
    for key in names:
        stars[key] = stars_df[key].values

    # Set offsets to apply. Doing a 3% systematic to make sure there is some good noise
    offsets = [OffsetSys(error_sys=0.03), OffsetSNR(lsst_filter=filter_name)]
    observed_stars = generate_catalog(visits, stars, offsets=offsets, lsst_filter=filter_name)

    # Save results
    np.savez('test_generate.npz', observed_stars=observed_stars)

    # now to try running the solution
    solver = LsqrSolver(iter_lim=10000)
    solver.run()
