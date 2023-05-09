Things changed since baseline 3.0

* Minor update to the footprint to include small Euclid swath
* Update survey start time to MJD=60676, Jan 1 2025
* uzy filter loader rather than just uz
* Updated to newer cloud database
* Telescope now parks with r-filter in. Should be no filter change time added when coming out of park.
* Updated long gaps survey to only attempt log gaps if early pair was sucessful. Also only attempt one starting blob in a night.
* Changing the blobs to not have az limits.
* Fixed bug in projection for solving blob path. Might be a slight slewtime improvement from it.
* Fix bug that was letting DDFs execute in twilight time

This directory has starter.py, which is the initial simulation for v3.1 that other experiements can be compared to.

The nominal idea is that after running a serries of different simulations, we will make decisions and merge into a final nominal baseline simulation at the end.

