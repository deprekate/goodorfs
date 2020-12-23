To apply patches navigate into the directory of the repo/tarball and run the command `git apply xxx.patch`

To patch Glimmer to use GOODORFS rather than LONGORFS, run the commands:
```
cd glimmer3.02
git apply /path/to/g3-from-scratch.patch
```

You can then run the Glimmer scripts as normal:
```
$ g3-from-scratch.csh NC_001416.fna NC_001416
Step 1 of 4:  Finding goodorfs for training
Step 2 of 4:  Extracting training sequences
Step 3 of 4:  Building ICM
Step 4 of 4:  Running Glimmer3
...
```
