To apply patches navigate into the directory of the repo/tarball and run the command `git apply xxx.patch`

An example for patching PHANOTATE to print the stop codons of the training gene set would be:
```
git clone https://github.com/deprekate/PHANOTATE.git
cd PHANOTATE
git apply /path/to/phanotate.patch
```

