Things changed since baseline 3.0

* minor update to the footprint to include small Euclid swath
* update survey start time
* uzy filter loader rather than just uz
* updated to newer cloud database
* telescope now parks with r-filter in. Should be no filter change time added when coming out of park.
* updated long gaps survey to only attempt log gaps if early pair was sucessful
* potentially changing the blobs to not have azand radius limits.

This directory has starter.py, which is the initial simulation for v3.1 that other experiements can be compared to.

The nominal idea is that after running a serries of different simulations, we will make decisions and merge into a final nominal simulation at the end.

