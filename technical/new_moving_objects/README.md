Trying out refactoring moving objects

for the original code:
time make_lsst_obs --simulation_db starter_v3.1_0yrs.db --orbit_file /Users/yoachim/rubin_sim_data/orbits/trojan_5k.txt

 360.78s user 7.69s system 99% cpu 6:08.52 total


After dumping the write:
354.65s user 6.53s system 99% cpu 6:01.19 total

So a modest improvement, but better than nothing. And this should make it possible for things to scale.

killing logging:
 354.21s user 4.77s system 99% cpu 5:59.00 total


 doing dmag with arrays:
 353.23s user 3.71s system 99% cpu 5:56.96 total



 ----
 run_moving_calc --obs_file starter_v3.1_0yrs__trojan_5k_obs.npz --simulation_db starter_v3.1_0yrs.db --orbit_file /Users/yoachim/rubin_sim_data/orbits/trojan_5k.txt --out_dir starter_v3.1_0yrs_ss --objtype Trojan --start_time 60218.0 

 run_moving_fractions --work_dir starter_v3.1_0yrs_ss --metadata Trojan --start_time 60218.0


 --

 make_lsst_obs --simulation_db starter_v3.1_0yrs.db --orbit_file /Users/yoachim/rubin_sim_data/orbits/mba_5k.txt --positions_file=/Users/yoachim/rubin_sim_data/orbits_precompute/mba_5k.npz
