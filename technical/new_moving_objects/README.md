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


time make_lsst_obs --simulation_db starter_v3.1_10yrs.db --orbit_file /Users/yoachim/rubin_sim_data/orbits/mba_5k.txt --positions_file=/Users/yoachim/rubin_sim_data/orbits_precompute/mba_5k.npz

2415.72s user 126.18s system 99% cpu 42:23.40 total
1650.16s user 140.29s system 99% cpu 29:50.55 total
1953.97s user 158.60s system 99% cpu 35:12.65 total -- this one gets all the obejcts!



now the original:
5572.44s user 204.47s system 99% cpu 1:36:17.70 total

So more than a factor of 2 speedup!


Memory footprint for MBAs is 1 Gb. Grows pretty slowly, so not worried about it getting over 2 Gb.



----

checking for matches between txt and npz files:
`ls *.txt | xargs -I'{}' egrep -cv '#|^$' '{}'`
results:
1606944  43387563 730489834 ender_a1_v3.1_10yrs__granvik_5k_obs.txt
  1599161  43177422 727096237 ender_a1_v3.1_10yrs__granvik_pha_5k_obs.txt
   2884780   77889135 1336153718 ender_a1_v3.1_10yrs__l7_5k_obs.txt
   2298870   62069565 1049016992 ender_a1_v3.1_10yrs__mba_5k_obs.txt
   2368150   63940125 1090281653 ender_a1_v3.1_10yrs__occ_rmax20_5k_obs.txt
   2292906   61908537 1052082178 ender_a1_v3.1_10yrs__occ_rmax5_5k_obs.txt
   2556825   69034350 1168500640 ender_a1_v3.1_10yrs__trojan_5k_obs.txt
   228866   6179457 103179556 ender_a1_v3.1_10yrs__vatiras_granvik_10k_obs.txt
(Note, there is one line in the text files for column names that is not commented)

1606987   ender_a1_v3.1_10yrs__granvik_5k_obs.npz
1599311   ender_a1_v3.1_10yrs__granvik_pha_5k_obs.npz
2884779   ender_a1_v3.1_10yrs__l7_5k_obs.npz
2298869   ender_a1_v3.1_10yrs__mba_5k_obs.npz
2368149   ender_a1_v3.1_10yrs__occ_rmax20_5k_obs.npz
2292913   ender_a1_v3.1_10yrs__occ_rmax5_5k_obs.npz
2556824   ender_a1_v3.1_10yrs__trojan_5k_obs.npz
228865    ender_a1_v3.1_10yrs__vatiras_granvik_10k_obs.npz


