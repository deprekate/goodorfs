Introduction
------------

GOODORFS is a tool to classify open reading frames into coding and noncoding

It takes as input a fasta file representing the entire genome.  It then finds all potential open reading frames, and then for each, calculates the Energy Density Profile from the amino acid frequency.

To install `GOODORFS`,
```
pip3 install goodorfs
```
or
```
git clone https://github.com/deprekate/goodorfs.git
cd goodorfs
python3 setup.py install

```

To run `GOODORFS` simply provide the path to a fasta file.  The default output is the same format as Glimmers LONGORFS program, in order to serve as a drop in replacement.
The columns are: orf_id, start_location, stop_location, frame, a bunch of zeros as filler
```
$ goodorfs.py tests/NC_001416.fna
00001     191     736  +2   0.000
00002     711    2636  +3   0.000
00003    2633    2839  +2   0.000
00004    3270    2830  -1   0.000
00005    2836    4437  +1   0.000
00006    5095    4604  -2   0.000
00007    4283    5737  +2   0.000
...
```

Additionally `GOODORFS` can also output the nucleotide sequences in fasta format for use in other applications:
```
$ good-orfs.py -Y fna tests/NC_001416.fna | head
>NC_001416_orf1 [START=191] [STOP=736]
ATGGAAGTCAACAAAAAGCAGCTGGCTGACATTTTCGGTGCGAGTATCCGTACCATTCA...
>NC_001416_orf2 [START=711] [STOP=2636]
GTGAATATATCGAACAGTCAGGTTAACAGGCTGCGGCATTTTGTCCGCGCCGGGCTTCG...
>NC_001416_orf3 [START=2633] [STOP=2839]
ATGACGCGACAGGAAGAACTTGCCGCTGCCCGTGCGGCACTGCATGACCTGATGACAGG...
>NC_001416_orf4 [START=3270] [STOP=2830]
GTGCATGGCCACACCTTCCCGAATCATCATGGTAAACGTGCGTTTTCGCTCAACGTCAA...
...
```

We have started testing GOODORFS to run on metagenomes. All that is needed is to bin reads according to their GC content and then run the bins through GOODROFS in batches in order to predict gene fragments within the reads. 

We have added a script to group the reads according to gc content. It prints out batches of 500 reads separated by the null terminator character, which allows commands to be chained to xargs.  To run GOODORFS on the supplied sample metagenome (which is in FASTA file format), run the command:
```
python3 scripts/bin_reads.py tests/ERR5004783_part.fasta | xargs -0 -I {} sh -c "echo '{}' | ./goodorfs.py -Y fna"
```

