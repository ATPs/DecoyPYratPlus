# DecoyPYratPlus

Download the the GitHub package directly

You need to install numpy.

## usage
python $PATH_of_decoyPYratPlus/decoyPYratPlus.py --help

`python /data/p/xiaolong/DecoyPYrat/decoypyrat/decoyPYratPlus.py  -h`
```ini
usage: decoyPYratPlus.py [-h] [--cleavage_sites CSITES] [--anti_cleavage_sites NOC] [--cleavage_position {c,n}] [--min_peptide_length MINLEN] [--max_peptide_length MAXLEN] [--max_iterations MAXIT]
                         [--miss_cleavage MISS_CLEAVAGE] [--do_not_shuffle] [--all_shuffle_mimic] [--do_not_switch] [--decoy_prefix DPREFIX] [--output_fasta DOUT] [--temp_file TOUT] [--no_isobaric]
                         [--memory_save] [--keep_names] [--target TARGET_FILE] [--checkSimilar]
                         *.fasta|*.fa|*.fasta.gz|*.fa.gz|*.txt|*.txt.gz [*.fasta|*.fa|*.fasta.gz|*.fa.gz|*.txt|*.txt.gz ...]

Create decoy protein sequences. Each protein is reversed and the cleavage sites switched with preceding amino acid. Peptides are checked for existence in target sequences if found the tool will attempt to
shuffle them. James.Wright@sanger.ac.uk 2015

positional arguments:
  *.fasta|*.fa|*.fasta.gz|*.fa.gz|*.txt|*.txt.gz
                        FASTA file of target proteins sequences for which to create decoys

optional arguments:
  -h, --help            show this help message and exit
  --cleavage_sites CSITES, -c CSITES
                        A list of amino acids at which to cleave during digestion. Default = KR
  --anti_cleavage_sites NOC, -a NOC
                        A list of amino acids at which not to cleave if following cleavage site ie. Proline. Default = none
  --cleavage_position {c,n}, -p {c,n}
                        Set cleavage to be c or n terminal of specified cleavage sites. Default = c
  --min_peptide_length MINLEN, -l MINLEN
                        Set minimum length of peptides to compare between target and decoy. Default = 6
  --max_peptide_length MAXLEN, -M MAXLEN
                        Set max length of peptides to compare between target and decoy. Default = 40
  --max_iterations MAXIT, -n MAXIT
                        Set maximum number of times to shuffle a peptide to make it non-target before failing. Default=100
  --miss_cleavage MISS_CLEAVAGE, -L MISS_CLEAVAGE
                        miss_cleavage when digesting protein. only work when checkSimilar is enabled. Default=2
  --do_not_shuffle, -x  Turn OFF shuffling of decoy peptides that are in the target database. Default=false
  --all_shuffle_mimic, -R
                        use random seq when generating decoy proteins. Similar method like mimic. Default=false
  --do_not_switch, -s   Turn OFF switching of cleavage site with preceding amino acid. Default=false
  --decoy_prefix DPREFIX, -d DPREFIX
                        Set accesion prefix for decoy proteins in output. Default=XXX
  --output_fasta DOUT, -o DOUT
                        Set file to write decoy proteins to. Default=decoy.fa
  --temp_file TOUT, -t TOUT
                        Set temporary file to write decoys prior to shuffling. Default=tmp.fa
  --no_isobaric, -i     Do not make decoy peptides isobaric. Default=false, I will be changed to L in decoy sequences
  --memory_save, -m     Slower but uses less memory (does not store decoy peptide list). Default=false
  --keep_names, -k      Keep sequence names in the decoy output. Default=false
  --target TARGET_FILE, -T TARGET_FILE
                        Combine and store the target file. I will be changed to L default. If no_isobaric, I will not be changed. Default="", do not save file
  --checkSimilar, -S    If set, ItoL is enabled automatically; output_fasta will include target sequences by changing I to L; allow overlapped digestion, and max_peptide_length will be used. In default
                        setting, the digested peptides do not overlap with each other. "Peptides are checked for existence in target sequences and if found the tool will attempt to shuffle them iterativly
                        until they are unique". Additionally, consider those amino acids were equal: N=D, Q=E, GG=N. If cannot solve after max_iterations of shuffling, introduce AA mutations, deletions or
                        insertions. AA not in cleavage_sites and anti_cleavage_sites. The number_of_changes (sum of mutations, deletions and insertions) <= 1 for each peptide. 3 * max_iterations times.
                        do_not_shuffle option will be ignored. Default=false

```

## updates compare to DecoyPYrat
Most decoy peptide were generated with the same method as DecoyPYrat
* gzip file supported, multiple input files
* if `checkSimilar` is set, DecoyPYratPlus further remove the possibility that a decoy peptide exists in the target database
  * consider those amino acids were equal: N=D, Q=E, GG=N. This will further reduce the chance that the decoy peptide is too similar to the target peptide.
  * allow miss-cleavage when get possible target peptides
    * DecoyPYrat do not include peptides from miss-cleaved sites when checking overlap of target and decoy peptides
  * DecoyPYratPlus introduces single mutation to the peptides which cannot be solved by shuffling. Mutation can be substitution, insertion or deletion. AA frequency in decoy database is almost unchanged. 
* an option of `--all_shuffle_mimic` is provided. The decoy sequences peptides will all be shuffled, nor just revert of the sequence.
  * **Note: It may be a bad idea, since shuffle every decoy peptides will make it much larger of the decoy searching space**

## how it works
* gzipped and multiple files were combined together and output to file defined by `--target` file. If not set, target sequencces won't be saved. Note: `I` will be changed to `L` in default.
* If `checkSimilar` is not set, it runs the same as `DecoyPYrat`
* If `checkSimilar` is enabled,
  *  if `all_shuffle_mimic` is not enabled
     *  it first run the same as `DecoyPYrat`
     *  target protein is digested with number of `miss_cleavage` allowed. Here, only cut at the C-terminal of cleavage sites. The target peptide filtered by `min_peptide_length` and `max_peptide_length`
     *  The decoy proteins were processed one by one. Get possible decoy peptides through digestion with miss-cleavage allowed. If after considering N=D, Q=E, GG=N, no decoy peptides were identified in target database, save the decoy protein. Otherwise, shuffle the decoy peptide only and do not change the other parts of the protein and get the new decoy protein. Keep the new decoy protein if has fewer peptides that can be identified in target database. Repeat at most 10 rounds until the decoy protein shares no common peptide with the target database. Repeat additionally `max_iterations * 2` rounds by shuffling decoy peptides with mutation allowed. Stop until the new decoy protein shared no peptides in target database.
  *  if `all_shuffle_mimic` is enabled, almost the same as above. the differences are:
     *  it first run the same as `DecoyPYrat`, but the peptides were shuffled instead of reversed when first generated.
     *  the max number of shuffling decoy peptides with mutation were changed to `max_iterations * 5`.

--------

# DecoyPYrat
DecoyPYrat - Fast Hybrid Decoy Sequence Database Creation for Proteomic Mass Spectrometery Analyses

Accurate statistical evaluation of sequence database peptide identifications from tandem mass spectra is essential in mass spectrometry based proteomics experiments. These statistics are dependent on accurately modelling random identifications.

The target-decoy approach has risen to become the de-facto approach to calculating false discovery rates (FDR) in proteomic datasets. The main principle of this approach is to search a set of decoy protein sequences that emulate the size and composition of the target protein sequences searched whilst not matching real proteins in the sample.

DecoyPYrat creates decoy protein sequences by following these steps: each protein is reversed and the cleavage sites switched with preceding amino acid. Peptides are checked for existence in target sequences and if found the tool will attempt to shuffle them iterativly until they are unique.

## Download and installation

### Bioconda

DecoyPYrat is available in the bioconda bioinformatics software repository. To access it, first install [Miniconda](https://docs.conda.io/en/latest/miniconda.html), and then run the command:

```shell
conda install -c bioconda decoypyrat
```

After this, you can invoke the software like this:

```shell
decoypyrat
```

### Direct script usage

You can clone this repository and invoke the software like this:

```shell
python decoypyrat/decoyPYrat.py
```

### Getting help

You can see the full usage instructions by specifying the "-h" argument:

If installing with Bioconda:
```shell
decoypyrat -h
```

If using the script directly:
```shell
python decoypyrat/decoyPYrat.py -h
```

## Citation:
[DecoyPyrat: Fast Non-redundant Hybrid Decoy Sequence Generation for Large Scale Proteomics.
J Proteomics Bioinform. 2016 Jun 27;9(6):176-180.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4941923/)
