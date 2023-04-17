

if __name__ == "__main__":
    runs = ['../starter/starter_v3.1_10yrs.db',
    '../noroll/noroll_v3.1_10yrs.db',
    "../selfcal_footprints/heavy_plane_v3.1_10yrs.db",
    "../selfcal_footprints/no_plane_v3.1_10yrs.db",
    "../selfcal_footprints/light_plane_v3.1_10yrs.db"]

    nights = [730.5, 365.25]

    filters = 'ugrizy'

    with open('selfcal.sh', 'w') as f:
        for night in nights:
            for filtername in filters:
                for run in runs:
                    print("./run_selfcal.py --pointings_db %s --filtername %s --night_max %i" % (run, filtername, night), file=f)
        for filtername in filters:
            print("./run_selfcal.py --pointings_db ../starter/starter_v3.1_10yrs.db --filtername %s --night_max 365 --n_patches 4" % filtername, file=f)