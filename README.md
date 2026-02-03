# DecoyPYratPlus

## install
Download the the GitHub package directly

`git clone https://github.com/ATPs/DecoyPYratPlus.git`

You need to install numpy, tqdm.
```bash
pip install numpy
pip install tqdm
```
or use conda/mamba

```bash
mamba  install -c conda-forge tqdm numpy
```

## usage

`python PATH_to_FOLDER_OF_DecoyPYratPlus/decoypyrat/decoyPYratPlus.py  -h`

```ini
usage: decoyPYratPlus.py [-h] [--cleavage_sites CSITES] [--anti_cleavage_sites NOC] [--cleavage_position {c,n}]
                         [--min_peptide_length MINLEN] [--max_peptide_length MAXLEN] [--max_iterations MAXIT]
                         [--miss_cleavage MISS_CLEAVAGE] [--do_not_shuffle] [--all_shuffle_mimic] [--do_not_switch]
                         [--decoy_prefix DPREFIX] [--output_fasta DOUT] [--no_isobaric] [--target_I2L]
                         [--fast_digest] [--threads THREADS] [--keep_names] [--target TARGET_FILE]
                         [--checkSimilar CHECKSIMILAR] [--concat CONCAT] [--dedup]
                         *.fasta|*.fa|*.fasta.gz|*.fa.gz|*.txt|*.txt.gz
                         [*.fasta|*.fa|*.fasta.gz|*.fa.gz|*.txt|*.txt.gz ...]

Create decoy protein sequences. Each protein is reversed and the cleavage sites switched with preceding amino acid.
Peptides are checked for existence in target sequences if found the tool will attempt to shuffle them.

positional arguments:
  *.fasta|*.fa|*.fasta.gz|*.fa.gz|*.txt|*.txt.gz
                        FASTA file of target proteins sequences for which to create decoys

optional arguments:
  -h, --help            show this help message and exit
  --cleavage_sites CSITES, -c CSITES
                        A list of amino acids at which to cleave during digestion. Default = KR
  --anti_cleavage_sites NOC, -a NOC
                        A list of amino acids at which not to cleave if following cleavage site ie. Proline.
                        Default = "P"
  --cleavage_position {c,n}, -p {c,n}
                        Set cleavage to be c or n terminal of specified cleavage sites. Default = c
  --min_peptide_length MINLEN, -l MINLEN
                        Set minimum length of peptides to compare between target and decoy. Default = 6
  --max_peptide_length MAXLEN, -M MAXLEN
                        Set max length of peptides to compare between target and decoy. Default = 40
  --max_iterations MAXIT, -n MAXIT
                        Set maximum number of times to shuffle a peptide to make it non-target before failing.
                        Default=100
  --miss_cleavage MISS_CLEAVAGE, -L MISS_CLEAVAGE
                        miss_cleavage when digesting protein. only work when checkSimilar is enabled. Default=2
  --do_not_shuffle, -x  Turn OFF shuffling of decoy peptides that are in the target database. Default=False
  --all_shuffle_mimic, -R
                        use random seq when generating decoy proteins. Similar method like mimic. Default=False
  --do_not_switch, -s   If set, revert the whole sequence as decoy. Default=False. Use the pseudo-reverse method
                        used by Comet.
  --decoy_prefix DPREFIX, -d DPREFIX
                        Set accesion prefix for decoy proteins in output. Default=DECOY
  --output_fasta DOUT, -o DOUT
                        Set file to write decoy proteins to. Default=decoy.fa
  --no_isobaric, -i     Do not make decoy peptides isobaric. Default=False, I will be changed to L in decoy
                        sequences
  --target_I2L          Convert I to L in target output and dedup. Default=False
  --fast_digest         Use Cython-accelerated digestion if available. Default=False
  --threads THREADS, -N THREADS
                        number of threads to use. default 1. Note: currently only one thread is allowed
  --keep_names, -k      Keep sequence names in the decoy output. Default=False. if decoy_prefix = "DECOY", name
                        will be like DECOY_1, DECOY_2. The orginal name will be included in the header line,
                        separated by " ". If set, name will be like DECOY_{oringal_name}
  --target TARGET_FILE, -T TARGET_FILE
                        Combine and store the target file. Target sequences are unchanged by default; use
                        --target_I2L to convert I to L. Default="", do not save file
  --checkSimilar CHECKSIMILAR, -S CHECKSIMILAR
                        If set, I/L are treated as the same unless --no_isobaric is set; target output I to L is
                        controlled by --target_I2L. To enable this option, it is recommended to set
                        checkSimilar="GG=N,N=D,Q=E". Allow missed cleavage, and max_peptide_length will be used
                        when determining shared peptides in target and decoy sequences. In default setting, the
                        digested peptides do not overlap with each other. "Peptides are checked for existence in
                        target sequences and if found the tool will attempt to shuffle them iterativly until they
                        are unique". Additionally, consider those amino acids were equal: N=D, Q=E, GG=N. If cannot
                        solve after max_iterations of shuffling, introduce AA mutations, deletions or insertions.
                        AA not in cleavage_sites and anti_cleavage_sites. The number_of_changes (sum of mutations,
                        deletions and insertions) <= 1 for each peptide. 2 * max_iterations times. do_not_shuffle
                        option will be ignored. Default=False
  --concat CONCAT, -C CONCAT
                        default = "", do not append target sequences to decoy sequences. If concat is set, the
                        target sequence is appended to decoy sequences, with concat as the separator. use '*' for
                        comet or other engines which treat "*" as separator. the same as cleavage_sites ("R" in
                        most cases) for other search engines
  --dedup, -U           If set, try to remove duplicated sequences in the target database. The first sequence will
                        be kept. default = False.

```

## optional: fast digest (Cython)
`--fast_digest` will use a Cython-accelerated digestion module if available. To build it (Cython must already be installed in your environment):

```bash
cd /data/p/xiaolong/DecoyPYrat
python setup_fast_digest.py build_ext --inplace
```

If the module cannot be built or imported, the program falls back to the pure Python implementation.

### Example benchmark (user-provided)
Command:
```bash
python /data/p/xiaolong/DecoyPYrat/decoypyrat/decoyPYratPlus.py /data2/pub/proteome/web/protinsight/comet/proteins/20260206/protinsight_proteinseq.fasta --output_fasta /data2/pub/proteome/web/protinsight/comet/proteins/20260206/protinsight_proteinseq.decoy.fasta --keep_names --checkSimilar "GG=N,N=D,Q=E" --min_peptide_length 7 --max_peptide_length 40 --all_shuffle_mimic
```

Pure Python (no `--fast_digest`): total time ~214.82 seconds. Key stages:
- digest input proteins allowing miss cleavage: ~29.10 s
- get decoy proteins with peptides in target proteins: ~49.49 s
- digest decoy proteins allowing miss cleavage: ~22.77 s

Same command with `--fast_digest`:
```bash
python /data/p/xiaolong/DecoyPYrat/decoypyrat/decoyPYratPlus.py /data2/pub/proteome/web/protinsight/comet/proteins/20260206/protinsight_proteinseq.fasta --output_fasta /data2/pub/proteome/web/protinsight/comet/proteins/20260206/protinsight_proteinseq.decoy.fasta --keep_names --checkSimilar "GG=N,N=D,Q=E" --min_peptide_length 7 --max_peptide_length 40 --all_shuffle_mimic --fast_digest
```

With `--fast_digest`: total time ~158.10 seconds. Key stages:
- digest input proteins allowing miss cleavage: ~14.87 s
- get decoy proteins with peptides in target proteins: ~39.00 s
- digest decoy proteins allowing miss cleavage: ~13.19 s

## updates compare to DecoyPYrat
Most decoy peptide were generated with the same method as DecoyPYrat
* gzip file supported, multiple input files were supported.
* if `checkSimilar` is set, DecoyPYratPlus further remove the possibility that a decoy peptide exists in the target database
  * consider those amino acids were equal: I=L (unless `--no_isobaric` is set), N=D, Q=E, GG=N. `--target_I2L` only affects target output, not the overlap check.
  * allow miss-cleavage when get possible target peptides
    * DecoyPYrat do not include peptides from miss-cleaved sites when checking overlap of target and decoy peptides
  * DecoyPYratPlus introduces single mutation to the peptides which cannot be solved by shuffling. Mutation can be substitution, insertion or deletion. AA frequency in decoy database is almost unchanged. 
  * Human Proteome Project (HPP) Mass Spectrometry Data Interpretation GuidelinesVersion 3.0.0 – October 15, 2019
    > Even when very high confidence peptide identifications are demonstrated, consider alternate mappings of the peptide to proteins other than the claimed one. Consider isobaric sequence/mass modification variants, all known SAAVs, and unreported SAAVs. Even when a peptide identification is shown to be very highly confident, care should be taken when mapping it to a protein or novel coding element. Consider whether I=L, N[Deamidated]=D, Q[Deamidated]=E, GG=N, Q≈K, F≈M[Oxidation], or other isobaric or near isobaric substitutions could change the mapping of the peptide from an extraordinary result to a mapping to a commonly-observed protein. Consider if a known single amino-acid variation (SAAV) in neXtProt could turn an extraordinary result into an ordinary result. Consider if a single amino-acid change, not yet annotated in a well-known source, could turn an extraordinary result into a questionable result. Check more than one reference proteome (e.g., RefSeq may have entries that UniProt and Ensembl do not, and vice versa). A tool to assist with this analysis is available at neXtProt at https://www.nextprot.org/tools/peptide-uniqueness-checker (Unicity Checker), and another at PeptideAtlas at http://peptideatlas.org/map (ProteoMapper).
  * set `--checkSimilar "GG=N,N=D,Q=E"` is recommended based on the HPP standard. 
    * The program will perform the replacement sequentially. That means, it will first replace all "GG" with "N", then replace all "N" with "D" and then replace all "Q" with "E".
    * The setting should be like `--checkSimilar "GG=N,N=D,Q=E,K=E,F=M"` if considering `Q≈K, F≈M`. If the database is big, it is hard to find proper decoy sequences if many amino acids were treated as equal. It seems that although `Q≈K, F≈M`, they are different enough. We would not recommend to use this setting.
* an option of `--all_shuffle_mimic` is provided. The decoy sequences peptides will all be shuffled, not just revert of the sequence.
  * **Note: Shuffle every decoy peptides will make the decoy searching space much larger. But here we try to generate the same decoy peptide for each unique target peptide so the decoy database contained slightly more searchable peptides than in the target database.**
  * In [mimic](https://github.com/msaid-de/mimic/tree/master)
    > Mimic combines the advantages of both methods by shuffling peptides in a manner that conserves homolog peptides. The method first identifies all homolog peptides, shuffles the unique peptides and then redistributes them to each protein again again. In this way a set of peptides that are homologs before shuffling will be homologs after shuffling as well https://github.com/percolator/mimic/wiki.
  * Here, more shuffle is done. Since miss-cleavage is allowed, the number of peptide to check is larger. A dictionary is created. The key is the decoy peptide that need to be altered (they were identified in the target database). The value is a list of the alternative decoy peptides. Decoy proteins were checked one by one, each time the alternative decoy peptide were stored in the dictionary, and will be re-used in the following runs. A new alternative decoy peptide will be created if all existing alternatives cannot solve the problem.
  * This option is useful one the "averaging strategy" is used. 
    > Keich, U., Tamura, K. & Noble, W. S. Averaging Strategy To Reduce Variability in Target-Decoy Estimates of False Discovery Rate. J Proteome Res 18, 585–593 (2019).



* The method for revert the sequences was changed. Now used a method the same to Comet.
  * Comet generates decoys by reversing each target peptide sequence, keeping the N-terminal or C-terminal amino acid in place (depending on the "sense" value of the digestion enzyme specified by search_enzyme_number). For example, peptide DIGSESTK becomes decoy peptide TSESGIDK for a tryptic search and peptide DVINHKGGA becomes DAGGKHNIV for an Asp-N search. https://uwpr.github.io/Comet/parameters/parameters_202301/decoy_search.html

* `--dedup` option added. If enabled, only one one of the duplicated sequences were kept.
* `--concat` option added. If set `--concat * `, the decoy sequence will be joined with the target sequence by symbol defined in `--concat`. 
  * Use `*` if the search engine like Comet treats `*` as cutting site. 
  * Use `R` for Trypsion digestion.
  * This option can be used for multistage searches. 
    > Ivanov, M. V., Levitsky, L. I. & Gorshkov, M. V. Adaptation of Decoy Fusion Strategy for Existing Multi-Stage Search Workflows. J. Am. Soc. Mass Spectrom. 27, 1579–1582 (2016).


## how it works
* gzipped and multiple files were combined together and output to file defined by `--target` file. If not set, target sequencces won't be saved. Note: target `I` is unchanged by default; use `--target_I2L` to convert `I` to `L`.
* If `checkSimilar` is not set, it runs the same as `DecoyPYrat`
* If `checkSimilar` is enabled,
  *  if `all_shuffle_mimic` is not enabled
     *  it first run the same as `DecoyPYrat`
     *  target protein is digested with number of `miss_cleavage` allowed. Here, only cut at the C-terminal of cleavage sites. The target peptide filtered by `min_peptide_length` and `max_peptide_length`
     *  The decoy proteins were processed one by one. Get possible decoy peptides through digestion with miss-cleavage allowed. If after considering N=D, Q=E, GG=N, no decoy peptides were identified in target database, save the decoy protein. Otherwise, shuffle the decoy peptide only and do not change the other parts of the protein and get the new decoy protein. Keep the new decoy protein if has fewer peptides that can be identified in target database. Repeat at most 10 rounds until the decoy protein shares no common peptide with the target database. Repeat additionally `max_iterations * 2` rounds by shuffling decoy peptides with mutation allowed. Stop until the new decoy protein shared no peptides in target database.
  *  if `all_shuffle_mimic` is enabled, almost the same as above. the differences are:
     *  it first run the same as `DecoyPYrat`. Then the reversed peptides were shuffled, so that each reversed peptide will be changed to the same randomized peptide.
     *  the max number of shuffling decoy peptides with mutation were changed to `max_iterations * 5`.
  *  For each decoy protein, all its peptides (allowing miss-cleavage and considering equal peptides defined by `checkSimilar`, like `I=L` unless `--no_isobaric` is set, `GG=N,N=D,Q=E`) were calculated and the peptides overlapped with the target database were stored in `ls_decoy_tochange`.
     *  Randomly choose one peptide in `ls_decoy_tochange`. 
        *  If peptide already existed in the dictionary `dAlternative2`, use previous alternative peptides (each peptide may have multiple alternatives).
        *  If not, generate an alternative decoy peptide that do not exist in the target database.
     *  Replace the original peptide in the decoy sequence with the new alternative decoy peptide. Get the count of overlapped decoy peptides. If the number is smaller, save the peptide and alternative peptide in `dAlternative2`, and change the decoy protein. Repeat the process until the decoy protein does not contain peptides that exists in the target database.
        *  If cannot solve, Regenerate the alternative peptide and allowing introducing AA insertion, deletion or substitution.
     

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
