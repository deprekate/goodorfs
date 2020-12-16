Introduction
------------

GOODORFS is a tool to classify open reading frames into coding and noncoding

It takes as input a fasta file representing the entire genome.  It then finds all potential open reading frames,  nd then for each, calculates the Energy Density Profile from the amino acid frequency.

To run `GOODORFS` simply provide the path to a fasta file.  The default output is the same format as Glimmers LONGORFS program, in order to serve as a drop in replacement.
The columns are: orf_id, start_location, stop_location, frame, a bunch of zeros as filler
```
$ python3 goodorfs.py tests/NC_001416.fna
Sequence file = tests/NC_001416.fna
Number of orfs = 502
00001     191     736  +2   0.000
00002     711    2636  +3   0.000
00003    2633    2839  +2   0.000
00004    3270    2830  -1   0.000
00005    2836    4437  +1   0.000
00006    5095    4604  -2   0.000
00007    4283    5737  +2   0.000
...
```


