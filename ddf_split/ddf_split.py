#!/usr/bin/env python

import numpy as np
import matplotlib.pylab as plt
import healpy as hp
from rubin_sim.scheduler.model_observatory import ModelObservatory
from rubin_sim.scheduler.schedulers import CoreScheduler, FilterSchedUzy
from rubin_sim.scheduler.utils import (
    EuclidOverlapFootprint,
    ConstantFootprint,
    make_rolling_footprints,
)
import rubin_sim.scheduler.basis_functions as bf
from rubin_sim.scheduler.surveys import (
    GreedySurvey,
    BlobSurvey,
    ScriptedSurvey,
    LongGapSurvey,
)
from rubin_sim.scheduler import sim_runner
import rubin_sim.scheduler.detailers as detailers
import sys
import subprocess
import os
import argparse
from astropy.coordinates import SkyCoord
from astropy import units as u
from rubin_sim.utils import _hpid2_ra_dec
import rubin_sim

from rubin_sim.utils import calc_season, ddf_locations
from rubin_sim.scheduler.utils import scheduled_observation
from rubin_sim.data import get_data_dir

# So things don't fail on hyak
from astropy.utils import iers

iers.conf.auto_download = False


def ddf_slopes(ddf_name, raw_obs, night_season):
    """
    Let's make custom slopes for each DDF

    Parameters
    ----------
    ddf_name : str
       The DDF name to use
    raw_obs : np.array
        An array with values of 1 or zero. One element per night, value of 1 indicates the night
        is during an active observing season.
    night_season : np.array
        An array of floats with the fractional season value (e.g., 0.5 would be half way through the first season)
    """

    # OK, so 258 sequences is ~1% of the survey
    # so a 25.8 sequences is a 0.1% season
    # COSMOS is going to be 0.7% for 3 years, then 0.175 for the rest.

    ss = 30  # standard season, was 45

    if (ddf_name == "ELAISS1") | (ddf_name == "XMM_LSS") | (ddf_name == "ECDFS"):
        # Dict with keys for each season and values of the number of sequences
        # to attempt.
        season_vals = {
            0: 10,
            1: ss,
            2: ss,
            3: ss,
            4: ss,
            5: ss,
            6: ss,
            7: ss,
            8: ss,
            9: ss,
            10: 10,
        }

        round_season = np.floor(night_season)

        cumulative_desired = np.zeros(raw_obs.size, dtype=float)
        for season in np.unique(round_season):
            in_season = np.where(round_season == season)
            cumulative = np.cumsum(raw_obs[in_season])
            cumulative = cumulative / cumulative.max() * season_vals[season]
            cumulative_desired[in_season] = cumulative + np.max(cumulative_desired)

    if ddf_name == "EDFS_a":
        season_vals = {
            0: 10,
            1: ss,
            2: ss,
            3: ss,
            4: ss,
            5: ss,
            6: ss,
            7: ss,
            8: ss,
            9: ss,
            10: 10,
        }

        round_season = np.floor(night_season)

        cumulative_desired = np.zeros(raw_obs.size, dtype=float)
        for season in np.unique(round_season):
            in_season = np.where(round_season == season)
            cumulative = np.cumsum(raw_obs[in_season])
            cumulative = cumulative / cumulative.max() * season_vals[season]
            cumulative_desired[in_season] = cumulative + np.max(cumulative_desired)

    if ddf_name == "COSMOS":
        # looks like COSMOS has no in-season time for 10 at the current start mjd.
        season_vals = {
            0: 10,
            1: ss * 5,
            2: ss * 5,
            3: ss * 2,
            4: ss,
            5: ss,
            6: ss,
            7: ss,
            8: ss,
            9: ss,
            10: 10,
        }

        round_season = np.floor(night_season)

        cumulative_desired = np.zeros(raw_obs.size, dtype=float)
        for season in np.unique(round_season):
            in_season = np.where(round_season == season)[0]
            cumulative = np.cumsum(raw_obs[in_season])
            if cumulative.max() > 0:
                cumulative = cumulative / cumulative.max() * season_vals[season]
                cumulative_desired[in_season] = cumulative + np.max(cumulative_desired)

    return cumulative_desired


def match_cumulative(cumulative_desired, mask=None, no_duplicate=True):
    """Generate a schedule that tries to match the desired cumulative distribution given a mask

    Parameters
    ----------
    cumulative_desired : `np.array`, float
        An array with the cumulative number of desired observations. Elements
        are assumed to be evenly spaced.
    mask : `np.array`, bool or int (None)
        Set to zero for indices that cannot be scheduled
    no_duplicate : bool (True)
        If True, only 1 event can be scheduled per element

    Returns
    -------
    schedule : `np.array`
        The resulting schedule, with values marking number of events in that cell.
    """

    rounded_desired = np.round(cumulative_desired)
    sched = cumulative_desired * 0
    if mask is None:
        mask = np.ones(sched.size)

    valid = np.where(mask > 0)[0].tolist()
    x = np.arange(sched.size)

    drd = np.diff(rounded_desired)
    step_points = np.where(drd > 0)[0] + 1

    # would be nice to eliminate this loop, but it's not too bad.
    # can't just use searchsorted on the whole array, because then there
    # can be duplicate values, and array[[n,n]] = 1 means that extra match gets lost.
    for indx in step_points:
        left = np.searchsorted(x[valid], indx)
        right = np.searchsorted(x[valid], indx, side="right")
        d1 = indx - left
        d2 = right - indx
        if d1 < d2:
            sched_at = left
        else:
            sched_at = right

        # If we are off the end
        if sched_at >= len(valid):
            sched_at -= 1

        sched[valid[sched_at]] += 1
        if no_duplicate:
            valid.pop(sched_at)

    return sched


def optimize_ddf_times(
    ddf_name,
    ddf_RA,
    ddf_grid,
    sun_limit=-18,
    airmass_limit=2.5,
    sky_limit=None,
    g_depth_limit=23.5,
    season_unobs_frac=0.1,
):
    """

    Parameters
    ----------
    ddf : `str`
        The name of the DDF
    ddf_grid : `np.array`
        An array with info for the DDFs. Generated by the rubin_sim/scheduler/surveys/generate_ddf_grid.py` script
    season_unobs_frac : `float`
        7.2 month observing season if season_unobs_frac = 0.2 (shaves 20% off each end of the full year)
    """
    sun_limit = np.radians(sun_limit)

    # XXX-- double check that I got this right
    ack = ddf_grid["sun_alt"][0:-1] * ddf_grid["sun_alt"][1:]
    night = np.zeros(ddf_grid.size, dtype=int)
    night[np.where((ddf_grid["sun_alt"][1:] >= 0) & (ack < 0))] += 1
    night = np.cumsum(night)
    ngrid = ddf_grid["mjd"].size

    # set a sun, airmass, sky masks
    sun_mask = np.ones(ngrid, dtype=int)
    sun_mask[np.where(ddf_grid["sun_alt"] >= sun_limit)] = 0

    airmass_mask = np.ones(ngrid, dtype=int)
    airmass_mask[np.where(ddf_grid["%s_airmass" % ddf_name] >= airmass_limit)] = 0

    sky_mask = np.ones(ngrid, dtype=int)
    if sky_limit is not None:
        sky_mask[np.where(ddf_grid["%s_sky_g" % ddf_name] <= sky_limit)] = 0
        sky_mask[np.where(np.isnan(ddf_grid["%s_sky_g" % ddf_name]) == True)] = 0

    m5_mask = np.zeros(ngrid, dtype=bool)
    m5_mask[np.isfinite(ddf_grid["%s_m5_g" % ddf_name])] = 1

    if g_depth_limit is not None:
        m5_mask[np.where(ddf_grid["%s_m5_g" % ddf_name] < g_depth_limit)] = 0

    big_mask = sun_mask * airmass_mask * sky_mask * m5_mask

    potential_nights = np.unique(night[np.where(big_mask > 0)])

    # prevent a repeat sequence in a night
    unights, indx = np.unique(night, return_index=True)
    night_mjd = ddf_grid["mjd"][indx]
    # The season of each night
    night_season = calc_season(ddf_RA, night_mjd)

    raw_obs = np.ones(unights.size)
    # take out the ones that are out of season
    season_mod = night_season % 1

    out_season = np.where(
        (season_mod < season_unobs_frac) | (season_mod > (1.0 - season_unobs_frac))
    )
    raw_obs[out_season] = 0

    cumulative_desired = ddf_slopes(ddf_name, raw_obs, night_season)

    night_mask = unights * 0
    night_mask[potential_nights] = 1

    unight_sched = match_cumulative(cumulative_desired, mask=night_mask)
    cumulative_sched = np.cumsum(unight_sched)

    nights_to_use = unights[np.where(unight_sched == 1)]

    # For each night, find the best time in the night.
    # XXX--probably need to expand this part to resolve the times when multiple things get scheduled
    mjds = []
    for night_check in nights_to_use:
        in_night = np.where(
            (night == night_check) & (np.isfinite(ddf_grid["%s_m5_g" % ddf_name]))
        )[0]
        m5s = ddf_grid["%s_m5_g" % ddf_name][in_night]
        # we could intorpolate this to get even better than 15 min resolution on when to observe
        max_indx = np.where(m5s == m5s.max())[0].min()
        mjds.append(ddf_grid["mjd"][in_night[max_indx]])

    return mjds, night_mjd, cumulative_desired, cumulative_sched


def split_times(mjds, ddf_name, ddf_grid, hours_lag=3.,
                g_depth_limit=23.5, sun_limit=-18, airmass_limit=2.5):
    """Take some mjds and split them if we can
    """

    lag_time = hours_lag/24.  # to days

    out_mjds = []
    split = []

    # XXX-Should probably eliminate loop
    for mjd in mjds:
        check_mjds = np.array([mjd-lag_time, mjd+lag_time])
        sun_alts = np.interp(check_mjds, ddf_grid['mjd'], ddf_grid['sun_alt'])
        g_depths = np.interp(check_mjds, ddf_grid['mjd'], ddf_grid["%s_m5_g" % ddf_name])
        airmasses = np.interp(check_mjds, ddf_grid['mjd'], ddf_grid["%s_airmass" % ddf_name])

        possible = np.where((sun_alts > sun_limit) & (g_depths > g_depth_limit) & (airmasses < airmass_limit))[0]
        if np.size(possible) == 0:
            out_mjds.append(mjd)
            split.append(False)
        elif np.size(possible) == 2:
            best = np.where(g_depths == np.max(g_depths))[0]
            possible = possible[best]
        if np.size(possible) == 1:
            out_mjds.append(mjd)
            out_mjds.append(np.max(check_mjds[possible]))
            split.extend([True, True])

    order = np.argsort(out_mjds)
    out_mjds = np.array(out_mjds)[order]
    split = np.array(split)[order]

    return out_mjds, split


def generate_ddf_scheduled_obs(
    data_file=None,
    flush_length=2,
    mjd_tol=5,
    expt=30.0,
    alt_min=25,
    alt_max=85,
    HA_min=21.0,
    HA_max=3.0,
    sun_alt_max=-18,
    dist_tol=3.0,
    season_unobs_frac=0.1,
    nvis_master=[8, 10, 20, 20, 24, 18],
    filters="ugrizy",
    nsnaps=[1, 2, 2, 2, 2, 2],
    hours_lag=3.
):
    """

    Parameters
    ----------
    data_file : path (None)
        The data file to use for DDF airmass, m5, etc. Defaults to using whatever is in
        rubin_sim_data/scheduler directory.
    flush_length : float (2)
        How long to keep a scheduled observation around before it is considered failed
        and flushed (days).
    mjd_tol : float (15)
        How close an observation must be in time to be considered matching a scheduled
        observation (minutes).
    expt : float (30)
        Total exposure time per visit (seconds).
    alt_min/max : float (25, 85)
        The minimum and maximum altitudes to permit observations to happen (degrees).
    HA_min/max : float (21, 3)
        The hour angle limits to permit observations to happen (hours).
    dist_tol : float (3)
        The distance tolerance for a visit to be considered matching a scheduled observation
        (degrees).
    season_unobs_frac : float (0.1)
        What fraction of the season should the DDF be considered unobservable. Taken off both the
        start and end of the year, so a season frac of 0.1 means 20% of the time the DDF is considered
        unobservable, so it will be in-season for 9.6 months.
    nvis_master : list of ints ([8, 10, 20, 20, 24, 18])
        The number of visits to make per filter
    filters : str (ugrizy)
        The filter names.
    nsnaps : list of ints ([1, 2, 2, 2, 2, 2])
        The number of snaps to use per filter
    hours_lag : float (2.)
        How much time to try and have between DDF sequences on the same night
    """
    if data_file is None:
        data_file = os.path.join(get_data_dir(), "scheduler", "ddf_grid.npz")

    flush_length = flush_length  # days
    mjd_tol = mjd_tol / 60 / 24.0  # minutes to days
    expt = expt
    alt_min = np.radians(alt_min)
    alt_max = np.radians(alt_max)
    dist_tol = np.radians(dist_tol)
    sun_alt_max = np.radians(sun_alt_max)

    ddfs = ddf_locations()
    ddf_data = np.load(data_file)
    ddf_grid = ddf_data["ddf_grid"].copy()

    all_scheduled_obs = []
    for ddf_name in ["ELAISS1", "XMM_LSS", "ECDFS", "COSMOS", "EDFS_a"]:
        print("Optimizing %s" % ddf_name)

        # 'ID', 'RA', 'dec', 'mjd', 'flush_by_mjd', 'exptime', 'filter', 'rotSkyPos', 'nexp',
        #         'note'
        # 'mjd_tol', 'dist_tol', 'alt_min', 'alt_max', 'HA_max', 'HA_min', 'observed'
        mjds = optimize_ddf_times(
            ddf_name,
            ddfs[ddf_name][0],
            ddf_grid,
            season_unobs_frac=season_unobs_frac,
        )[0]

        # now we have observations at the best times per night.
        # for each of those, try to add an mjd N hours before or after
        mjds, splits = split_times(mjds, ddf_name, ddf_grid, hours_lag=hours_lag)

        for mjd, split in zip(mjds, splits):
            for filtername, nvis, nexp in zip(filters, nvis_master, nsnaps):
                if "EDFS" in ddf_name:
                    n_final = nvis / 2 # Double field
                    if split:
                        n_final = n_final/2.

                    obs = scheduled_observation(n=int(n_final))
                    obs["RA"] = np.radians(ddfs[ddf_name][0])
                    obs["dec"] = np.radians(ddfs[ddf_name][1])
                    obs["mjd"] = mjd
                    obs["flush_by_mjd"] = mjd + flush_length
                    obs["exptime"] = expt
                    obs["filter"] = filtername
                    obs["nexp"] = nexp
                    obs["note"] = "DD:%s" % ddf_name

                    obs["mjd_tol"] = mjd_tol
                    obs["dist_tol"] = dist_tol
                    # Need to set something for HA limits
                    obs["HA_min"] = HA_min
                    obs["HA_max"] = HA_max
                    obs["alt_min"] = alt_min
                    obs["alt_max"] = alt_max
                    obs["sun_alt_max"] = sun_alt_max
                    all_scheduled_obs.append(obs)

                    obs = scheduled_observation(n=int(n_final))
                    obs["RA"] = np.radians(ddfs[ddf_name.replace("_a", "_b")][0])
                    obs["dec"] = np.radians(ddfs[ddf_name.replace("_a", "_b")][1])
                    obs["mjd"] = mjd
                    obs["flush_by_mjd"] = mjd + flush_length
                    obs["exptime"] = expt
                    obs["filter"] = filtername
                    obs["nexp"] = nexp
                    obs["note"] = "DD:%s" % ddf_name.replace("_a", "_b")

                    obs["mjd_tol"] = mjd_tol
                    obs["dist_tol"] = dist_tol
                    # Need to set something for HA limits
                    obs["HA_min"] = HA_min
                    obs["HA_max"] = HA_max
                    obs["alt_min"] = alt_min
                    obs["alt_max"] = alt_max
                    obs["sun_alt_max"] = sun_alt_max
                    all_scheduled_obs.append(obs)

                else:
                    n_final = nvis
                    if split:
                        n_final = n_final/2.
                    obs = scheduled_observation(n=int(n_final))
                    obs["RA"] = np.radians(ddfs[ddf_name][0])
                    obs["dec"] = np.radians(ddfs[ddf_name][1])
                    obs["mjd"] = mjd
                    obs["flush_by_mjd"] = mjd + flush_length
                    obs["exptime"] = expt
                    obs["filter"] = filtername
                    obs["nexp"] = nexp
                    obs["note"] = "DD:%s" % ddf_name

                    obs["mjd_tol"] = mjd_tol
                    obs["dist_tol"] = dist_tol
                    # Need to set something for HA limits
                    obs["HA_min"] = HA_min
                    obs["HA_max"] = HA_max
                    obs["alt_min"] = alt_min
                    obs["alt_max"] = alt_max
                    obs["sun_alt_max"] = sun_alt_max
                    all_scheduled_obs.append(obs)

    result = np.concatenate(all_scheduled_obs)
    # Put in the scripted ID so it's easier to track which ones fail.
    result["scripted_id"] = np.arange(result.size)
    return result


def blob_for_long(
    nside,
    nexp=2,
    exptime=30.0,
    filter1s=["g"],
    filter2s=["i"],
    pair_time=33.0,
    camera_rot_limits=[-80.0, 80.0],
    n_obs_template=None,
    season=300.0,
    season_start_hour=-4.0,
    season_end_hour=2.0,
    shadow_minutes=60.0,
    max_alt=76.0,
    moon_distance=30.0,
    ignore_obs=["DD", "twilight_neo"],
    m5_weight=6.0,
    footprint_weight=1.5,
    slewtime_weight=3.0,
    stayfilter_weight=3.0,
    template_weight=12.0,
    u_template_weight=24.0,
    footprints=None,
    u_nexp1=True,
    night_pattern=[True, True],
    time_after_twi=30.0,
    HA_min=12,
    HA_max=24 - 3.5,
    nights_delayed=-1,
):
    """
    Generate surveys that take observations in blobs.

    Parameters
    ----------
    nside : int (32)
        The HEALpix nside to use
    nexp : int (1)
        The number of exposures to use in a visit.
    exptime : float (30.)
        The exposure time to use per visit (seconds)
    filter1s : list of str
        The filternames for the first set
    filter2s : list of str
        The filter names for the second in the pair (None if unpaired)
    pair_time : float (33)
        The ideal time between pairs (minutes)
    camera_rot_limits : list of float ([-80., 80.])
        The limits to impose when rotationally dithering the camera (degrees).
    n_obs_template : int (5)
        The number of observations to take every season in each filter
    season : float (300)
        The length of season (i.e., how long before templates expire) (days)
    season_start_hour : float (-4.)
        For weighting how strongly a template image needs to be observed (hours)
    sesason_end_hour : float (2.)
        For weighting how strongly a template image needs to be observed (hours)
    shadow_minutes : float (60.)
        Used to mask regions around zenith (minutes)
    max_alt : float (76.
        The maximium altitude to use when masking zenith (degrees)
    moon_distance : float (30.)
        The mask radius to apply around the moon (degrees)
    ignore_obs : str or list of str ('DD')
        Ignore observations by surveys that include the given substring(s).
    m5_weight : float (3.)
        The weight for the 5-sigma depth difference basis function
    footprint_weight : float (0.3)
        The weight on the survey footprint basis function.
    slewtime_weight : float (3.)
        The weight on the slewtime basis function
    stayfilter_weight : float (3.)
        The weight on basis function that tries to stay avoid filter changes.
    template_weight : float (12.)
        The weight to place on getting image templates every season
    u_template_weight : float (24.)
        The weight to place on getting image templates in u-band. Since there
        are so few u-visits, it can be helpful to turn this up a little higher than
        the standard template_weight kwarg.
    u_nexp1 : bool (True)
        Add a detailer to make sure the number of expossures in a visit is always 1 for u observations.
    """

    BlobSurvey_params = {
        "slew_approx": 7.5,
        "filter_change_approx": 140.0,
        "read_approx": 2.0,
        "min_pair_time": 15.0,
        "search_radius": 30.0,
        "alt_max": 85.0,
        "az_range": 90.0,
        "flush_time": 30.0,
        "smoothing_kernel": None,
        "nside": nside,
        "seed": 42,
        "dither": True,
        "twilight_scale": True,
    }

    surveys = []
    if n_obs_template is None:
        n_obs_template = {"u": 3, "g": 3, "r": 3, "i": 3, "z": 3, "y": 3}

    times_needed = [pair_time, pair_time * 2]
    for filtername, filtername2 in zip(filter1s, filter2s):
        detailer_list = []
        detailer_list.append(
            detailers.CameraRotDetailer(
                min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
            )
        )
        detailer_list.append(detailers.CloseAltDetailer())
        # List to hold tuples of (basis_function_object, weight)
        bfs = []

        if filtername2 is not None:
            bfs.append(
                (
                    bf.M5DiffBasisFunction(filtername=filtername, nside=nside),
                    m5_weight / 2.0,
                )
            )
            bfs.append(
                (
                    bf.M5DiffBasisFunction(filtername=filtername2, nside=nside),
                    m5_weight / 2.0,
                )
            )

        else:
            bfs.append(
                (bf.M5DiffBasisFunction(filtername=filtername, nside=nside), m5_weight)
            )

        if filtername2 is not None:
            bfs.append(
                (
                    bf.FootprintBasisFunction(
                        filtername=filtername,
                        footprint=footprints,
                        out_of_bounds_val=np.nan,
                        nside=nside,
                    ),
                    footprint_weight / 2.0,
                )
            )
            bfs.append(
                (
                    bf.FootprintBasisFunction(
                        filtername=filtername2,
                        footprint=footprints,
                        out_of_bounds_val=np.nan,
                        nside=nside,
                    ),
                    footprint_weight / 2.0,
                )
            )
        else:
            bfs.append(
                (
                    bf.FootprintBasisFunction(
                        filtername=filtername,
                        footprint=footprints,
                        out_of_bounds_val=np.nan,
                        nside=nside,
                    ),
                    footprint_weight,
                )
            )

        bfs.append(
            (
                bf.SlewtimeBasisFunction(filtername=filtername, nside=nside),
                slewtime_weight,
            )
        )
        bfs.append(
            (bf.StrictFilterBasisFunction(filtername=filtername), stayfilter_weight)
        )

        if filtername2 is not None:
            bfs.append(
                (
                    bf.NObsPerYearBasisFunction(
                        filtername=filtername,
                        nside=nside,
                        footprint=footprints.get_footprint(filtername),
                        n_obs=n_obs_template[filtername],
                        season=season,
                        season_start_hour=season_start_hour,
                        season_end_hour=season_end_hour,
                    ),
                    template_weight / 2.0,
                )
            )
            bfs.append(
                (
                    bf.NObsPerYearBasisFunction(
                        filtername=filtername2,
                        nside=nside,
                        footprint=footprints.get_footprint(filtername2),
                        n_obs=n_obs_template[filtername2],
                        season=season,
                        season_start_hour=season_start_hour,
                        season_end_hour=season_end_hour,
                    ),
                    template_weight / 2.0,
                )
            )
        else:
            bfs.append(
                (
                    bf.NObsPerYearBasisFunction(
                        filtername=filtername,
                        nside=nside,
                        footprint=footprints.get_footprint(filtername),
                        n_obs=n_obs_template[filtername],
                        season=season,
                        season_start_hour=season_start_hour,
                        season_end_hour=season_end_hour,
                    ),
                    template_weight,
                )
            )
        # Masks, give these 0 weight
        bfs.append(
            (
                bf.ZenithShadowMaskBasisFunction(
                    nside=nside,
                    shadow_minutes=shadow_minutes,
                    max_alt=max_alt,
                    penalty=np.nan,
                    site="LSST",
                ),
                0.0,
            )
        )
        bfs.append(
            (
                bf.MoonAvoidanceBasisFunction(nside=nside, moon_distance=moon_distance),
                0.0,
            )
        )
        filternames = [fn for fn in [filtername, filtername2] if fn is not None]
        bfs.append((bf.FilterLoadedBasisFunction(filternames=filternames), 0))
        if filtername2 is None:
            time_needed = times_needed[0]
        else:
            time_needed = times_needed[1]
        bfs.append((bf.TimeToTwilightBasisFunction(time_needed=time_needed), 0.0))
        bfs.append((bf.NotTwilightBasisFunction(), 0.0))
        bfs.append((bf.PlanetMaskBasisFunction(nside=nside), 0.0))
        bfs.append((bf.AfterEveningTwiBasisFunction(time_after=time_after_twi), 0.0))
        # XXX--move kwargs up
        bfs.append((bf.HaMaskBasisFunction(ha_min=HA_min, ha_max=HA_max), 0.0))
        # don't execute every night
        bfs.append((bf.NightModuloBasisFunction(night_pattern), 0.0))
        # possibly force things to delay
        bfs.append((bf.DelayStartBasisFunction(nights_delay=nights_delayed), 0.0))

        # unpack the basis functions and weights
        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]
        if filtername2 is None:
            survey_name = "blob_long, %s" % filtername
        else:
            survey_name = "blob_long, %s%s" % (filtername, filtername2)
        if filtername2 is not None:
            detailer_list.append(detailers.TakeAsPairsDetailer(filtername=filtername2))

        if u_nexp1:
            detailer_list.append(detailers.FilterNexp(filtername="u", nexp=1))
        surveys.append(
            BlobSurvey(
                basis_functions,
                weights,
                filtername1=filtername,
                filtername2=filtername2,
                exptime=exptime,
                ideal_pair_time=pair_time,
                survey_note=survey_name,
                ignore_obs=ignore_obs,
                nexp=nexp,
                detailers=detailer_list,
                **BlobSurvey_params
            )
        )

    return surveys


def gen_long_gaps_survey(
    footprints,
    nside=32,
    night_pattern=[True, True],
    gap_range=[2, 7],
    HA_min=12,
    HA_max=24 - 3.5,
    time_after_twi=120,
    nights_delayed=-1,
):
    """
    Paramterers
    -----------
    HA_min(_max) : float
        The hour angle limits passed to the initial blob scheduler.
    """

    surveys = []
    f1 = ["g", "r", "i"]
    f2 = ["r", "i", "z"]
    # Maybe force scripted to not go in twilight?

    for filtername1, filtername2 in zip(f1, f2):
        blob = blob_for_long(
            footprints=footprints,
            nside=nside,
            filter1s=[filtername1],
            filter2s=[filtername2],
            night_pattern=night_pattern,
            time_after_twi=time_after_twi,
            HA_min=HA_min,
            HA_max=HA_max,
            nights_delayed=nights_delayed,
        )
        scripted = ScriptedSurvey([], nside=nside, ignore_obs=["blob", "DDF", "twi"])
        surveys.append(
            LongGapSurvey(blob[0], scripted, gap_range=gap_range, avoid_zenith=True)
        )

    return surveys


def gen_GreedySurveys(
    nside=32,
    nexp=2,
    exptime=30.0,
    filters=["r", "i", "z", "y"],
    camera_rot_limits=[-80.0, 80.0],
    shadow_minutes=60.0,
    max_alt=76.0,
    moon_distance=30.0,
    ignore_obs=["DD", "twilight_neo"],
    m5_weight=3.0,
    footprint_weight=0.75,
    slewtime_weight=3.0,
    stayfilter_weight=3.0,
    repeat_weight=-1.0,
    footprints=None,
):
    """
    Make a quick set of greedy surveys

    This is a convienence function to generate a list of survey objects that can be used with
    rubin_sim.scheduler.schedulers.Core_scheduler.
    To ensure we are robust against changes in the sims_featureScheduler codebase, all kwargs are
    explicitly set.

    Parameters
    ----------
    nside : int (32)
        The HEALpix nside to use
    nexp : int (1)
        The number of exposures to use in a visit.
    exptime : float (30.)
        The exposure time to use per visit (seconds)
    filters : list of str (['r', 'i', 'z', 'y'])
        Which filters to generate surveys for.
    camera_rot_limits : list of float ([-80., 80.])
        The limits to impose when rotationally dithering the camera (degrees).
    shadow_minutes : float (60.)
        Used to mask regions around zenith (minutes)
    max_alt : float (76.
        The maximium altitude to use when masking zenith (degrees)
    moon_distance : float (30.)
        The mask radius to apply around the moon (degrees)
    ignore_obs : str or list of str ('DD')
        Ignore observations by surveys that include the given substring(s).
    m5_weight : float (3.)
        The weight for the 5-sigma depth difference basis function
    footprint_weight : float (0.3)
        The weight on the survey footprint basis function.
    slewtime_weight : float (3.)
        The weight on the slewtime basis function
    stayfilter_weight : float (3.)
        The weight on basis function that tries to stay avoid filter changes.
    """
    # Define the extra parameters that are used in the greedy survey. I
    # think these are fairly set, so no need to promote to utility func kwargs
    greed_survey_params = {
        "block_size": 1,
        "smoothing_kernel": None,
        "seed": 42,
        "camera": "LSST",
        "dither": True,
        "survey_name": "greedy",
    }

    surveys = []
    detailer_list = [
        detailers.CameraRotDetailer(
            min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
        )
    ]
    detailer_list.append(detailers.Rottep2RotspDesiredDetailer())

    for filtername in filters:
        bfs = []
        bfs.append(
            (bf.M5DiffBasisFunction(filtername=filtername, nside=nside), m5_weight)
        )
        bfs.append(
            (
                bf.FootprintBasisFunction(
                    filtername=filtername,
                    footprint=footprints,
                    out_of_bounds_val=np.nan,
                    nside=nside,
                ),
                footprint_weight,
            )
        )
        bfs.append(
            (
                bf.SlewtimeBasisFunction(filtername=filtername, nside=nside),
                slewtime_weight,
            )
        )
        bfs.append(
            (bf.StrictFilterBasisFunction(filtername=filtername), stayfilter_weight)
        )
        bfs.append(
            (
                bf.VisitRepeatBasisFunction(
                    gap_min=0, gap_max=2 * 60.0, filtername=None, nside=nside, npairs=20
                ),
                repeat_weight,
            )
        )
        # Masks, give these 0 weight
        bfs.append(
            (
                bf.ZenithShadowMaskBasisFunction(
                    nside=nside, shadow_minutes=shadow_minutes, max_alt=max_alt
                ),
                0,
            )
        )
        bfs.append(
            (bf.MoonAvoidanceBasisFunction(nside=nside, moon_distance=moon_distance), 0)
        )

        bfs.append((bf.FilterLoadedBasisFunction(filternames=filtername), 0))
        bfs.append((bf.PlanetMaskBasisFunction(nside=nside), 0))

        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]
        surveys.append(
            GreedySurvey(
                basis_functions,
                weights,
                exptime=exptime,
                filtername=filtername,
                nside=nside,
                ignore_obs=ignore_obs,
                nexp=nexp,
                detailers=detailer_list,
                **greed_survey_params
            )
        )

    return surveys


def generate_blobs(
    nside,
    nexp=2,
    exptime=30.0,
    filter1s=["u", "u", "g", "r", "i", "z", "y"],
    filter2s=["g", "r", "r", "i", "z", "y", "y"],
    pair_time=33.0,
    camera_rot_limits=[-80.0, 80.0],
    n_obs_template=None,
    season=300.0,
    season_start_hour=-4.0,
    season_end_hour=2.0,
    shadow_minutes=60.0,
    max_alt=76.0,
    moon_distance=30.0,
    ignore_obs=["DD", "twilight_neo"],
    m5_weight=6.0,
    footprint_weight=1.5,
    slewtime_weight=3.0,
    stayfilter_weight=3.0,
    template_weight=12.0,
    u_template_weight=24.0,
    footprints=None,
    u_nexp1=True,
    scheduled_respect=45.0,
    good_seeing={"g": 3, "r": 3, "i": 3},
    good_seeing_weight=3.0,
    mjd_start=1,
    repeat_weight=-20,
):
    """
    Generate surveys that take observations in blobs.

    Parameters
    ----------
    nside : int (32)
        The HEALpix nside to use
    nexp : int (1)
        The number of exposures to use in a visit.
    exptime : float (30.)
        The exposure time to use per visit (seconds)
    filter1s : list of str
        The filternames for the first set
    filter2s : list of str
        The filter names for the second in the pair (None if unpaired)
    pair_time : float (33)
        The ideal time between pairs (minutes)
    camera_rot_limits : list of float ([-80., 80.])
        The limits to impose when rotationally dithering the camera (degrees).
    n_obs_template : int (3)
        The number of observations to take every season in each filter
    season : float (300)
        The length of season (i.e., how long before templates expire) (days)
    season_start_hour : float (-4.)
        For weighting how strongly a template image needs to be observed (hours)
    sesason_end_hour : float (2.)
        For weighting how strongly a template image needs to be observed (hours)
    shadow_minutes : float (60.)
        Used to mask regions around zenith (minutes)
    max_alt : float (76.
        The maximium altitude to use when masking zenith (degrees)
    moon_distance : float (30.)
        The mask radius to apply around the moon (degrees)
    ignore_obs : str or list of str ('DD')
        Ignore observations by surveys that include the given substring(s).
    m5_weight : float (3.)
        The weight for the 5-sigma depth difference basis function
    footprint_weight : float (0.3)
        The weight on the survey footprint basis function.
    slewtime_weight : float (3.)
        The weight on the slewtime basis function
    stayfilter_weight : float (3.)
        The weight on basis function that tries to stay avoid filter changes.
    template_weight : float (12.)
        The weight to place on getting image templates every season
    u_template_weight : float (24.)
        The weight to place on getting image templates in u-band. Since there
        are so few u-visits, it can be helpful to turn this up a little higher than
        the standard template_weight kwarg.
    u_nexp1 : bool (True)
        Add a detailer to make sure the number of expossures in a visit is always 1 for u observations.
    scheduled_respect : float (45)
        How much time to require there be before a pre-scheduled observation (minutes)
    """

    template_weights = {
        "u": u_template_weight,
        "g": template_weight,
        "r": template_weight,
        "i": template_weight,
        "z": template_weight,
        "y": template_weight,
    }

    BlobSurvey_params = {
        "slew_approx": 7.5,
        "filter_change_approx": 140.0,
        "read_approx": 2.0,
        "min_pair_time": 15.0,
        "search_radius": 30.0,
        "alt_max": 85.0,
        "az_range": 90.0,
        "flush_time": 30.0,
        "smoothing_kernel": None,
        "nside": nside,
        "seed": 42,
        "dither": True,
        "twilight_scale": False,
    }

    if n_obs_template is None:
        n_obs_template = {"u": 3, "g": 3, "r": 3, "i": 3, "z": 3, "y": 3}

    surveys = []

    times_needed = [pair_time, pair_time * 2]
    for filtername, filtername2 in zip(filter1s, filter2s):
        detailer_list = []
        detailer_list.append(
            detailers.CameraRotDetailer(
                min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
            )
        )
        detailer_list.append(detailers.Rottep2RotspDesiredDetailer())
        detailer_list.append(detailers.CloseAltDetailer())
        detailer_list.append(detailers.FlushForSchedDetailer())
        # List to hold tuples of (basis_function_object, weight)
        bfs = []

        if filtername2 is not None:
            bfs.append(
                (
                    bf.M5DiffBasisFunction(filtername=filtername, nside=nside),
                    m5_weight / 2.0,
                )
            )
            bfs.append(
                (
                    bf.M5DiffBasisFunction(filtername=filtername2, nside=nside),
                    m5_weight / 2.0,
                )
            )

        else:
            bfs.append(
                (bf.M5DiffBasisFunction(filtername=filtername, nside=nside), m5_weight)
            )

        if filtername2 is not None:
            bfs.append(
                (
                    bf.FootprintBasisFunction(
                        filtername=filtername,
                        footprint=footprints,
                        out_of_bounds_val=np.nan,
                        nside=nside,
                    ),
                    footprint_weight / 2.0,
                )
            )
            bfs.append(
                (
                    bf.FootprintBasisFunction(
                        filtername=filtername2,
                        footprint=footprints,
                        out_of_bounds_val=np.nan,
                        nside=nside,
                    ),
                    footprint_weight / 2.0,
                )
            )
        else:
            bfs.append(
                (
                    bf.FootprintBasisFunction(
                        filtername=filtername,
                        footprint=footprints,
                        out_of_bounds_val=np.nan,
                        nside=nside,
                    ),
                    footprint_weight,
                )
            )

        bfs.append(
            (
                bf.SlewtimeBasisFunction(filtername=filtername, nside=nside),
                slewtime_weight,
            )
        )
        bfs.append(
            (bf.StrictFilterBasisFunction(filtername=filtername), stayfilter_weight)
        )
        bfs.append(
            (
                bf.VisitRepeatBasisFunction(
                    gap_min=0, gap_max=3 * 60.0, filtername=None, nside=nside, npairs=20
                ),
                repeat_weight,
            )
        )

        if filtername2 is not None:
            bfs.append(
                (
                    bf.NObsPerYearBasisFunction(
                        filtername=filtername,
                        nside=nside,
                        footprint=footprints.get_footprint(filtername),
                        n_obs=n_obs_template[filtername],
                        season=season,
                        season_start_hour=season_start_hour,
                        season_end_hour=season_end_hour,
                    ),
                    template_weights[filtername] / 2.0,
                )
            )
            bfs.append(
                (
                    bf.NObsPerYearBasisFunction(
                        filtername=filtername2,
                        nside=nside,
                        footprint=footprints.get_footprint(filtername2),
                        n_obs=n_obs_template[filtername2],
                        season=season,
                        season_start_hour=season_start_hour,
                        season_end_hour=season_end_hour,
                    ),
                    template_weights[filtername2] / 2.0,
                )
            )
        else:
            bfs.append(
                (
                    bf.NObsPerYearBasisFunction(
                        filtername=filtername,
                        nside=nside,
                        footprint=footprints.get_footprint(filtername),
                        n_obs=n_obs_template[filtername],
                        season=season,
                        season_start_hour=season_start_hour,
                        season_end_hour=season_end_hour,
                    ),
                    template_weight,
                )
            )

        # Insert things for getting good seeing templates
        if filtername2 is not None:
            if filtername in list(good_seeing.keys()):
                bfs.append(
                    (
                        bf.NGoodSeeingBasisFunction(
                            filtername=filtername,
                            nside=nside,
                            mjd_start=mjd_start,
                            footprint=footprints.get_footprint(filtername),
                            n_obs_desired=good_seeing[filtername],
                        ),
                        good_seeing_weight,
                    )
                )
            if filtername2 in list(good_seeing.keys()):
                bfs.append(
                    (
                        bf.NGoodSeeingBasisFunction(
                            filtername=filtername2,
                            nside=nside,
                            mjd_start=mjd_start,
                            footprint=footprints.get_footprint(filtername2),
                            n_obs_desired=good_seeing[filtername2],
                        ),
                        good_seeing_weight,
                    )
                )
        else:
            if filtername in list(good_seeing.keys()):
                bfs.append(
                    (
                        bf.NGoodSeeingBasisFunction(
                            filtername=filtername,
                            nside=nside,
                            mjd_start=mjd_start,
                            footprint=footprints.get_footprint(filtername),
                            n_obs_desired=good_seeing[filtername],
                        ),
                        good_seeing_weight,
                    )
                )
        # Make sure we respect scheduled observations
        bfs.append((bf.TimeToScheduledBasisFunction(time_needed=scheduled_respect), 0))
        # Masks, give these 0 weight
        bfs.append(
            (
                bf.ZenithShadowMaskBasisFunction(
                    nside=nside,
                    shadow_minutes=shadow_minutes,
                    max_alt=max_alt,
                    penalty=np.nan,
                    site="LSST",
                ),
                0.0,
            )
        )
        bfs.append(
            (
                bf.MoonAvoidanceBasisFunction(nside=nside, moon_distance=moon_distance),
                0.0,
            )
        )
        filternames = [fn for fn in [filtername, filtername2] if fn is not None]
        bfs.append((bf.FilterLoadedBasisFunction(filternames=filternames), 0))
        if filtername2 is None:
            time_needed = times_needed[0]
        else:
            time_needed = times_needed[1]
        bfs.append((bf.TimeToTwilightBasisFunction(time_needed=time_needed), 0.0))
        bfs.append((bf.NotTwilightBasisFunction(), 0.0))
        bfs.append((bf.PlanetMaskBasisFunction(nside=nside), 0.0))

        # unpack the basis functions and weights
        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]
        if filtername2 is None:
            survey_name = "blob, %s" % filtername
        else:
            survey_name = "blob, %s%s" % (filtername, filtername2)
        if filtername2 is not None:
            detailer_list.append(detailers.TakeAsPairsDetailer(filtername=filtername2))

        if u_nexp1:
            detailer_list.append(detailers.FilterNexp(filtername="u", nexp=1))
        surveys.append(
            BlobSurvey(
                basis_functions,
                weights,
                filtername1=filtername,
                filtername2=filtername2,
                exptime=exptime,
                ideal_pair_time=pair_time,
                survey_note=survey_name,
                ignore_obs=ignore_obs,
                nexp=nexp,
                detailers=detailer_list,
                **BlobSurvey_params
            )
        )

    return surveys


def generate_twi_blobs(
    nside,
    nexp=2,
    exptime=30.0,
    filter1s=["r", "i", "z", "y"],
    filter2s=["i", "z", "y", "y"],
    pair_time=15.0,
    camera_rot_limits=[-80.0, 80.0],
    n_obs_template=None,
    season=300.0,
    season_start_hour=-4.0,
    season_end_hour=2.0,
    shadow_minutes=60.0,
    max_alt=76.0,
    moon_distance=30.0,
    ignore_obs=["DD", "twilight_neo"],
    m5_weight=6.0,
    footprint_weight=1.5,
    slewtime_weight=3.0,
    stayfilter_weight=3.0,
    template_weight=12.0,
    footprints=None,
    repeat_night_weight=None,
    wfd_footprint=None,
    scheduled_respect=15.0,
    repeat_weight=-1.0,
    night_pattern=None,
):
    """
    Generate surveys that take observations in blobs.

    Parameters
    ----------
    nside : int (32)
        The HEALpix nside to use
    nexp : int (1)
        The number of exposures to use in a visit.
    exptime : float (30.)
        The exposure time to use per visit (seconds)
    filter1s : list of str
        The filternames for the first set
    filter2s : list of str
        The filter names for the second in the pair (None if unpaired)
    pair_time : float (22)
        The ideal time between pairs (minutes)
    camera_rot_limits : list of float ([-80., 80.])
        The limits to impose when rotationally dithering the camera (degrees).
    n_obs_template : int (3)
        The number of observations to take every season in each filter
    season : float (300)
        The length of season (i.e., how long before templates expire) (days)
    season_start_hour : float (-4.)
        For weighting how strongly a template image needs to be observed (hours)
    sesason_end_hour : float (2.)
        For weighting how strongly a template image needs to be observed (hours)
    shadow_minutes : float (60.)
        Used to mask regions around zenith (minutes)
    max_alt : float (76.
        The maximium altitude to use when masking zenith (degrees)
    moon_distance : float (30.)
        The mask radius to apply around the moon (degrees)
    ignore_obs : str or list of str ('DD')
        Ignore observations by surveys that include the given substring(s).
    m5_weight : float (3.)
        The weight for the 5-sigma depth difference basis function
    footprint_weight : float (0.3)
        The weight on the survey footprint basis function.
    slewtime_weight : float (3.)
        The weight on the slewtime basis function
    stayfilter_weight : float (3.)
        The weight on basis function that tries to stay avoid filter changes.
    template_weight : float (12.)
        The weight to place on getting image templates every season
    u_template_weight : float (24.)
        The weight to place on getting image templates in u-band. Since there
        are so few u-visits, it can be helpful to turn this up a little higher than
        the standard template_weight kwarg.
    """

    BlobSurvey_params = {
        "slew_approx": 7.5,
        "filter_change_approx": 140.0,
        "read_approx": 2.0,
        "min_pair_time": 10.0,
        "search_radius": 30.0,
        "alt_max": 85.0,
        "az_range": 90.0,
        "flush_time": 30.0,
        "smoothing_kernel": None,
        "nside": nside,
        "seed": 42,
        "dither": True,
        "twilight_scale": False,
        "in_twilight": True,
    }

    surveys = []

    if n_obs_template is None:
        n_obs_template = {"u": 3, "g": 3, "r": 3, "i": 3, "z": 3, "y": 3}

    times_needed = [pair_time, pair_time * 2]
    for filtername, filtername2 in zip(filter1s, filter2s):
        detailer_list = []
        detailer_list.append(
            detailers.CameraRotDetailer(
                min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
            )
        )
        detailer_list.append(detailers.Rottep2RotspDesiredDetailer())
        detailer_list.append(detailers.CloseAltDetailer())
        detailer_list.append(detailers.FlushForSchedDetailer())
        # List to hold tuples of (basis_function_object, weight)
        bfs = []

        if filtername2 is not None:
            bfs.append(
                (
                    bf.M5DiffBasisFunction(filtername=filtername, nside=nside),
                    m5_weight / 2.0,
                )
            )
            bfs.append(
                (
                    bf.M5DiffBasisFunction(filtername=filtername2, nside=nside),
                    m5_weight / 2.0,
                )
            )

        else:
            bfs.append(
                (bf.M5DiffBasisFunction(filtername=filtername, nside=nside), m5_weight)
            )

        if filtername2 is not None:
            bfs.append(
                (
                    bf.FootprintBasisFunction(
                        filtername=filtername,
                        footprint=footprints,
                        out_of_bounds_val=np.nan,
                        nside=nside,
                    ),
                    footprint_weight / 2.0,
                )
            )
            bfs.append(
                (
                    bf.FootprintBasisFunction(
                        filtername=filtername2,
                        footprint=footprints,
                        out_of_bounds_val=np.nan,
                        nside=nside,
                    ),
                    footprint_weight / 2.0,
                )
            )
        else:
            bfs.append(
                (
                    bf.FootprintBasisFunction(
                        filtername=filtername,
                        footprint=footprints,
                        out_of_bounds_val=np.nan,
                        nside=nside,
                    ),
                    footprint_weight,
                )
            )

        bfs.append(
            (
                bf.SlewtimeBasisFunction(filtername=filtername, nside=nside),
                slewtime_weight,
            )
        )
        bfs.append(
            (bf.StrictFilterBasisFunction(filtername=filtername), stayfilter_weight)
        )
        bfs.append(
            (
                bf.VisitRepeatBasisFunction(
                    gap_min=0, gap_max=2 * 60.0, filtername=None, nside=nside, npairs=20
                ),
                repeat_weight,
            )
        )

        if filtername2 is not None:
            bfs.append(
                (
                    bf.NObsPerYearBasisFunction(
                        filtername=filtername,
                        nside=nside,
                        footprint=footprints.get_footprint(filtername),
                        n_obs=n_obs_template[filtername],
                        season=season,
                        season_start_hour=season_start_hour,
                        season_end_hour=season_end_hour,
                    ),
                    template_weight / 2.0,
                )
            )
            bfs.append(
                (
                    bf.NObsPerYearBasisFunction(
                        filtername=filtername2,
                        nside=nside,
                        footprint=footprints.get_footprint(filtername2),
                        n_obs=n_obs_template[filtername2],
                        season=season,
                        season_start_hour=season_start_hour,
                        season_end_hour=season_end_hour,
                    ),
                    template_weight / 2.0,
                )
            )
        else:
            bfs.append(
                (
                    bf.NObsPerYearBasisFunction(
                        filtername=filtername,
                        nside=nside,
                        footprint=footprints.get_footprint(filtername),
                        n_obs=n_obs_template[filtername],
                        season=season,
                        season_start_hour=season_start_hour,
                        season_end_hour=season_end_hour,
                    ),
                    template_weight,
                )
            )
        if repeat_night_weight is not None:
            bfs.append(
                (
                    bf.AvoidLongGapsBasisFunction(
                        nside=nside,
                        filtername=None,
                        min_gap=0.0,
                        max_gap=10.0 / 24.0,
                        ha_limit=3.5,
                        footprint=wfd_footprint,
                    ),
                    repeat_night_weight,
                )
            )
        # Make sure we respect scheduled observations
        bfs.append((bf.TimeToScheduledBasisFunction(time_needed=scheduled_respect), 0))
        # Masks, give these 0 weight
        bfs.append(
            (
                bf.ZenithShadowMaskBasisFunction(
                    nside=nside,
                    shadow_minutes=shadow_minutes,
                    max_alt=max_alt,
                    penalty=np.nan,
                    site="LSST",
                ),
                0.0,
            )
        )
        bfs.append(
            (
                bf.MoonAvoidanceBasisFunction(nside=nside, moon_distance=moon_distance),
                0.0,
            )
        )
        filternames = [fn for fn in [filtername, filtername2] if fn is not None]
        bfs.append((bf.FilterLoadedBasisFunction(filternames=filternames), 0))
        if filtername2 is None:
            time_needed = times_needed[0]
        else:
            time_needed = times_needed[1]
        bfs.append(
            (bf.TimeToTwilightBasisFunction(time_needed=time_needed, alt_limit=12), 0.0)
        )
        bfs.append((bf.PlanetMaskBasisFunction(nside=nside), 0.0))

        # Let's turn off twilight blobs on nights where we are
        # doing NEO hunts
        bfs.append((bf.NightModuloBasisFunction(pattern=night_pattern), 0))

        # unpack the basis functions and weights
        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]
        if filtername2 is None:
            survey_name = "blob_twi, %s" % filtername
        else:
            survey_name = "blob_twi, %s%s" % (filtername, filtername2)
        if filtername2 is not None:
            detailer_list.append(detailers.TakeAsPairsDetailer(filtername=filtername2))
        surveys.append(
            BlobSurvey(
                basis_functions,
                weights,
                filtername1=filtername,
                filtername2=filtername2,
                exptime=exptime,
                ideal_pair_time=pair_time,
                survey_note=survey_name,
                ignore_obs=ignore_obs,
                nexp=nexp,
                detailers=detailer_list,
                **BlobSurvey_params
            )
        )

    return surveys


def ddf_surveys(detailers=None, season_unobs_frac=0.2, euclid_detailers=None):
    obs_array = generate_ddf_scheduled_obs(season_unobs_frac=season_unobs_frac)

    euclid_obs = np.where(
        (obs_array["note"] == "DD:EDFS_b") | (obs_array["note"] == "DD:EDFS_a")
    )[0]
    all_other = np.where(
        (obs_array["note"] != "DD:EDFS_b") & (obs_array["note"] != "DD:EDFS_a")
    )[0]

    survey1 = ScriptedSurvey([], detailers=detailers)
    survey1.set_script(obs_array[all_other])

    survey2 = ScriptedSurvey([], detailers=euclid_detailers)
    survey2.set_script(obs_array[euclid_obs])

    return [survey1, survey2]


def ecliptic_target(nside=32, dist_to_eclip=40.0, dec_max=30.0, mask=None):
    """Generate a target_map for the area around the ecliptic"""

    ra, dec = _hpid2_ra_dec(nside, np.arange(hp.nside2npix(nside)))
    result = np.zeros(ra.size)
    coord = SkyCoord(ra=ra * u.rad, dec=dec * u.rad)
    eclip_lat = coord.barycentrictrueecliptic.lat.radian
    good = np.where(
        (np.abs(eclip_lat) < np.radians(dist_to_eclip)) & (dec < np.radians(dec_max))
    )
    result[good] += 1

    if mask is not None:
        result *= mask

    return result


def generate_twilight_neo(
    nside,
    night_pattern=None,
    nexp=1,
    exptime=15,
    ideal_pair_time=5.0,
    max_airmass=2.0,
    camera_rot_limits=[-80.0, 80.0],
    time_needed=10,
    footprint_mask=None,
    footprint_weight=0.1,
    slewtime_weight=3.0,
    stayfilter_weight=3.0,
    area_required=None,
    filters="riz",
    n_repeat=4,
    sun_alt_limit=-14.8,
    slew_estimate=4.5,
):
    """Generate a survey for observing NEO objects in twilight

    Parameters
    ----------
    night_pattern : list of bool (None)
        A list of bools that set when the survey will be active. e.g., [True, False]
        for every-other night, [True, False, False] for every third night.
    nexp : int (1)
        Number of snaps in a visit
    exptime : float (15)
        Exposure time of visits
    ideal_pair_time : float (5)
        Ideal time between repeat visits (minutes).
    max_airmass : float (2)
        Maximum airmass to attempt (unitless).
    camera_rot_limits : list of float ([-80, 80])
        The camera rotation limits to use (degrees).
    time_needed : float (10)
        How much time should be available (e.g., before twilight ends) (minutes).
    footprint_mask : np.array (None)
        Mask to apply to the constructed ecliptic target mask (None).
    footprint_weight : float (0.1)
        Weight for footprint basis function
    slewtime_weight : float (3.)
        Weight for slewtime basis function
    stayfilter_weight : float (3.)
        Weight for staying in the same filter basis function
    area_required : float (None)
        The area that needs to be available before the survey will return observations (sq degrees?)
    filters : str ('riz')
        The filters to use, default 'riz'
    n_repeat : int (4)
        The number of times a blob should be repeated, default 4.
    sun_alt_limit : float (-14.8)
        Do not start unless sun is higher than this limit (degrees)
    slew_estimate : float (4.5)
        An estimate of how long it takes to slew between neighboring fields (seconds).
    """
    # XXX finish eliminating magic numbers and document this one
    moon_distance = 30.0
    shadow_minutes = 60.0
    max_alt = 76.0
    max_elong = 60.0
    az_range = 180.0

    survey_name = "twilight_neo"
    footprint = ecliptic_target(nside=nside, mask=footprint_mask)
    constant_fp = ConstantFootprint()
    for filtername in filters:
        constant_fp.set_footprint(filtername, footprint)

    surveys = []
    for filtername in filters:
        detailer_list = []
        detailer_list.append(
            detailers.CameraRotDetailer(
                min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
            )
        )
        detailer_list.append(detailers.CloseAltDetailer())
        # Should put in a detailer so things start at lowest altitude
        detailer_list.append(
            detailers.TwilightTripleDetailer(
                slew_estimate=slew_estimate, n_repeat=n_repeat
            )
        )
        bfs = []

        bfs.append(
            (
                bf.FootprintBasisFunction(
                    filtername=filtername,
                    footprint=constant_fp,
                    out_of_bounds_val=np.nan,
                    nside=nside,
                ),
                footprint_weight,
            )
        )

        bfs.append(
            (
                bf.SlewtimeBasisFunction(filtername=filtername, nside=nside),
                slewtime_weight,
            )
        )
        bfs.append(
            (bf.StrictFilterBasisFunction(filtername=filtername), stayfilter_weight)
        )
        # Need a toward the sun, reward high airmass, with an airmass cutoff basis function.
        bfs.append(
            (bf.NearSunTwilightBasisFunction(nside=nside, max_airmass=max_airmass), 0)
        )
        bfs.append(
            (
                bf.ZenithShadowMaskBasisFunction(
                    nside=nside, shadow_minutes=shadow_minutes, max_alt=max_alt
                ),
                0,
            )
        )
        bfs.append(
            (bf.MoonAvoidanceBasisFunction(nside=nside, moon_distance=moon_distance), 0)
        )
        bfs.append((bf.FilterLoadedBasisFunction(filternames=filtername), 0))
        bfs.append((bf.PlanetMaskBasisFunction(nside=nside), 0))
        bfs.append(
            (
                bf.SolarElongationMaskBasisFunction(
                    min_elong=0.0, max_elong=max_elong, nside=nside
                ),
                0,
            )
        )

        bfs.append((bf.NightModuloBasisFunction(pattern=night_pattern), 0))
        # Do not attempt unless the sun is getting high
        bfs.append(((bf.SunAltHighLimitBasisFunction(alt_limit=sun_alt_limit)), 0))

        # unpack the basis functions and weights
        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]

        # Set huge ideal pair time and use the detailer to cut down the list of observations to fit twilight?
        surveys.append(
            BlobSurvey(
                basis_functions,
                weights,
                filtername1=filtername,
                filtername2=None,
                ideal_pair_time=ideal_pair_time,
                nside=nside,
                exptime=exptime,
                survey_note=survey_name,
                ignore_obs=["DD", "greedy", "blob"],
                dither=True,
                nexp=nexp,
                detailers=detailer_list,
                az_range=az_range,
                twilight_scale=False,
                area_required=area_required,
            )
        )
    return surveys


def run_sched(
    surveys,
    survey_length=365.25,
    nside=32,
    fileroot="baseline_",
    verbose=False,
    extra_info=None,
    illum_limit=40.0,
    mjd_start=60676.0,
):
    years = np.round(survey_length / 365.25)
    scheduler = CoreScheduler(surveys, nside=nside)
    n_visit_limit = None
    fs = FilterSchedUzy(illum_limit=illum_limit)
    observatory = ModelObservatory(nside=nside, mjd_start=mjd_start)
    observatory, scheduler, observations = sim_runner(
        observatory,
        scheduler,
        survey_length=survey_length,
        filename=fileroot + "%iyrs.db" % years,
        delete_past=True,
        n_visit_limit=n_visit_limit,
        verbose=verbose,
        extra_info=extra_info,
        filter_scheduler=fs,
    )

    return observatory, scheduler, observations


def main(args):
    survey_length = args.survey_length  # Days
    outDir = args.outDir
    verbose = args.verbose
    max_dither = args.maxDither
    illum_limit = args.moon_illum_limit
    nexp = args.nexp
    nslice = args.rolling_nslice
    rolling_scale = args.rolling_strength
    dbroot = args.dbroot
    gsw = args.gsw
    nights_off = args.nights_off
    nights_delayed = args.nights_delayed
    neo_night_pattern = args.neo_night_pattern
    neo_filters = args.neo_filters
    neo_repeat = args.neo_repeat
    ddf_season_frac = args.ddf_season_frac

    # Be sure to also update and regenerate DDF grid save file if changing mjd_start
    mjd_start = 60676.0
    nside = 32
    per_night = True  # Dither DDF per night

    camera_ddf_rot_limit = 75.0  # degrees

    extra_info = {}
    exec_command = ""
    for arg in sys.argv:
        exec_command += " " + arg
    extra_info["exec command"] = exec_command
    try:
        extra_info["git hash"] = subprocess.check_output(["git", "rev-parse", "HEAD"])
    except subprocess.CalledProcessError:
        extra_info["git hash"] = "Not in git repo"

    extra_info["file executed"] = os.path.realpath(__file__)
    try:
        rs_path = rubin_sim.__path__[0]
        hash_file = os.path.join(rs_path, "../", ".git/refs/heads/main")
        extra_info["rubin_sim git hash"] = subprocess.check_output(["cat", hash_file])
    except subprocess.CalledProcessError:
        pass

    # Use the filename of the script to name the output database
    if dbroot is None:
        fileroot = os.path.basename(sys.argv[0]).replace(".py", "") + "_"
    else:
        fileroot = dbroot + "_"
    file_end = "v3.1_"

    pattern_dict = {
        1: [True],
        2: [True, False],
        3: [True, False, False],
        4: [True, False, False, False],
        # 4 on, 4 off
        5: [True, True, True, True, False, False, False, False],
        # 3 on 4 off
        6: [True, True, True, False, False, False, False],
        7: [True, True, False, False, False, False],
    }
    neo_night_pattern = pattern_dict[neo_night_pattern]
    reverse_neo_night_pattern = [not val for val in neo_night_pattern]

    # Modify the footprint
    sky = EuclidOverlapFootprint(nside=nside, smc_radius=4, lmc_radius=6)
    footprints_hp_array, labels = sky.return_maps()

    wfd_indx = np.where(
        (labels == "lowdust") | (labels == "LMC_SMC") | (labels == "virgo")
    )[0]
    wfd_footprint = footprints_hp_array["r"] * 0
    wfd_footprint[wfd_indx] = 1

    footprints_hp = {}
    for key in footprints_hp_array.dtype.names:
        footprints_hp[key] = footprints_hp_array[key]

    footprint_mask = footprints_hp["r"] * 0
    footprint_mask[np.where(footprints_hp["r"] > 0)] = 1

    repeat_night_weight = None

    observatory = ModelObservatory(nside=nside, mjd_start=mjd_start)
    conditions = observatory.return_conditions()

    footprints = make_rolling_footprints(
        fp_hp=footprints_hp,
        mjd_start=conditions.mjd_start,
        sun_ra_start=conditions.sun_ra_start,
        nslice=nslice,
        scale=rolling_scale,
        nside=nside,
        wfd_indx=wfd_indx,
        order_roll=1,
        n_cycles=4,
    )

    gaps_night_pattern = [True] + [False] * nights_off

    long_gaps = gen_long_gaps_survey(
        nside=nside,
        footprints=footprints,
        night_pattern=gaps_night_pattern,
        nights_delayed=nights_delayed,
    )

    # Set up the DDF surveys to dither
    u_detailer = detailers.FilterNexp(filtername="u", nexp=1)
    dither_detailer = detailers.DitherDetailer(
        per_night=per_night, max_dither=max_dither
    )
    details = [
        detailers.CameraRotDetailer(
            min_rot=-camera_ddf_rot_limit, max_rot=camera_ddf_rot_limit
        ),
        dither_detailer,
        u_detailer,
        detailers.Rottep2RotspDesiredDetailer(),
    ]
    euclid_detailers = [
        detailers.CameraRotDetailer(
            min_rot=-camera_ddf_rot_limit, max_rot=camera_ddf_rot_limit
        ),
        detailers.EuclidDitherDetailer(),
        u_detailer,
        detailers.Rottep2RotspDesiredDetailer(),
    ]
    ddfs = ddf_surveys(
        detailers=details,
        season_unobs_frac=ddf_season_frac,
        euclid_detailers=euclid_detailers,
    )

    greedy = gen_GreedySurveys(nside, nexp=nexp, footprints=footprints)
    neo = generate_twilight_neo(
        nside,
        night_pattern=neo_night_pattern,
        filters=neo_filters,
        n_repeat=neo_repeat,
        footprint_mask=footprint_mask,
    )
    blobs = generate_blobs(
        nside,
        nexp=nexp,
        footprints=footprints,
        mjd_start=conditions.mjd_start,
        good_seeing_weight=gsw,
    )
    twi_blobs = generate_twi_blobs(
        nside,
        nexp=nexp,
        footprints=footprints,
        wfd_footprint=wfd_footprint,
        repeat_night_weight=repeat_night_weight,
        night_pattern=reverse_neo_night_pattern,
    )
    surveys = [ddfs, long_gaps, blobs, twi_blobs, neo, greedy]
    observatory, scheduler, observations = run_sched(
        surveys,
        survey_length=survey_length,
        verbose=verbose,
        fileroot=os.path.join(outDir, fileroot + file_end),
        extra_info=extra_info,
        nside=nside,
        illum_limit=illum_limit,
        mjd_start=mjd_start,
    )

    return observatory, scheduler, observations


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", dest="verbose", action="store_true")
    parser.set_defaults(verbose=False)
    parser.add_argument("--survey_length", type=float, default=365.25 * 10)
    parser.add_argument("--outDir", type=str, default="")
    parser.add_argument(
        "--maxDither", type=float, default=0.7, help="Dither size for DDFs (deg)"
    )
    parser.add_argument(
        "--moon_illum_limit",
        type=float,
        default=40.0,
        help="illumination limit to remove u-band",
    )
    parser.add_argument("--nexp", type=int, default=2)
    parser.add_argument("--rolling_nslice", type=int, default=2)
    parser.add_argument("--rolling_strength", type=float, default=0.9)
    parser.add_argument("--dbroot", type=str)
    parser.add_argument("--gsw", type=float, default=3.0, help="good seeing weight")
    parser.add_argument("--ddf_season_frac", type=float, default=0.2)
    parser.add_argument("--nights_off", type=int, default=6)
    parser.add_argument("--nights_delayed", type=int, default=-1)
    parser.add_argument("--neo_night_pattern", type=int, default=4)
    parser.add_argument("--neo_filters", type=str, default="riz")
    parser.add_argument("--neo_repeat", type=int, default=4)

    args = parser.parse_args()
    main(args)
