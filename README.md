# mkLTG
**Taxonomic assignment of metabarcoding sequences using variable % identity thresholds**

If you are using this tool, please, cite Meglécz, E. (2024). mkLTG: A command-line tool for taxonomic assignment of metabarcoding sequences using variable identity thresholds. [**Biologia Futura**.](https://link.springer.com/article/10.1007/s42977-024-00201-x) [Accepted manuscript](https://amu.hal.science/hal-04434060)

## Installation

**Perl interpreter** 
This is generally installed in all unix like systems. Otherwise see https://www.activestate.com/products/perl/

**ncbi-blast+**
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/


**mkLTG**
Clone the github repository
~~~~
git clone https://github.com/meglecz/mkLTG.git
~~~~

If you are not comfortable with git, just download the zipped package (**Code** at the top right of the github page https://github.com/meglecz/mkLTG, then **Download ZIP**), dezip the file and you are ready.



## Reference database and taxonomy files

### COI database
A ready to use COI database with its associated taxonomy file formated from [COInr](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13756) for mkLTG (and [VTAM](https://www.sciencedirect.com/science/article/pii/S200103702300034X?via%3Dihub) ) is available here:  https://osf.io/vrfwz/ . The COInr database is a comprehensive non-redundant database with COI sequences from NCBI-nt and BOLD databases.

Alternativelly, it is possible to download the COInr database from zenodo https://zenodo.org/record/6555985, and prepare a custom database using mkCOInr (https://github.com/meglecz/mkCOInr)

### Other databases
Any BLAST databases can be formated from a fasta file and a taxmap file using the following command:
~~~~
makeblastdb -dbtype nucl -in FASTA_FILE -parse_seqids -taxid_map TAXMAP_FILE
~~~~

Format of the **taxmap file** (Sequence identifiers and taxIDs separated by tabulation):

~~~~
HQ449147_1	1000596
6930156	1000596
HQ449148_1	1000597
...
~~~~

mkLTG also needs a **taxonomy** tsv with the following columns:

- tax_id: can be correct NCBI taxID, or arbitrary negative numbers for taxa not in NCBI
- parent_tax_id: taxiID of the closest parent of tax_id
- rank: taxonomic rank (e.g. species, genus, subgenus, no_rank)
- name_txt: Scientifc name of the taxon
- old_tax_id: taxIDs that have been merged to the tax_id by NCBI; if more than one for a given tax_id, make one line for each old_tax_id (can be empty)
- taxlevel index: 0 => root, 1=> superkingdom, 2=> kingdom, 3=> phylum, 4=> class, 5=> order, 6=> family, 7=> genus, 8=> species; Levels in between have 0.5 added to the next highest level (e.g. 5.5 for infraorder and for superfamily).
- synonyms: list of synonyms; optional

~~~
tax_id  parent_tax_id   rank    name_txt        old_tax_id      taxlevel        synonyms
1       1       no rank root            0
2       131567  superkingdom    Bacteria                1       Prokaryotae;Prokaryota;Procaryotae
6       335928  genus   Azorhizobium            7
7       6       species Azorhizobium caulinodans        395     8       Azotirhizobium caulinodans
-34968  2778801 subfamily       Callithamnioideae               6.5
-35035  -35034  species Fractonotus caelatus            8
-35036  -35030  family  Ramazzottidae           6
~~~

When all sequences in the database have NCBI taxIDs, the taxonomy file can be automatically produced by mkLTG form the NCBI tax dmp files (https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/)


## mkLTG

### Test files
For the sake of this tutorial I have the following file system:

~~~~
~/mkLTG
├── COInr_for_vtam_2022_05_06_dbV5
│   ├── COInr_for_vtam.ndb
│   ├── COInr_for_vtam.nhr
│   ├── COInr_for_vtam.nin
│   ├── COInr_for_vtam.nog
│   ├── COInr_for_vtam.nos
│   ├── COInr_for_vtam.not
│   ├── COInr_for_vtam.nsq
│   ├── COInr_for_vtam.ntf
│   ├── COInr_for_vtam.nto
│   ├── COInr_for_vtam_taxonomy.tsv
│   └── format_db.log
├── data
│   └── test.fas
│   ├── test.tsv
├── params
│   └── params_default.tsv
├── README.md
└── scripts
    ├── compair_params.R
    ├── make_param_files_pid_fix.pl
    ├── make_param_files_pid_variable.pl
    ├── make_test_files.pl
    ├── mkLTG.pl
    ├── mkLTG.pm
    └── run_ltg_with_multiple_param_setting.pl

~~~~

The files in COInr_for_vtam_2022_05_06_dbV5 were downloaded form https://osf.io/vrfwz/

All other files are provided with mkLTG.

### Run
**Unix like systems**

~~~
cd mkLTG

perl scripts/mkLTG.pl -in data/test.fas -taxonomy COInr_for_vtam_2022_05_06_dbV5/COInr_for_vtam_taxonomy.tsv -blast_db COInr_for_vtam_2022_05_06_dbV5/COInr_for_vtam -outdir out -out_name test -ltg_params params/params_default.tsv
~~~

The input file (-in) can also be a tsv file with headings. One of the columns must have 'sequence' as a heading.


**Windows**

~~~
cd mkLTG

perl scripts\mkLTG.pl -in data\test.fas -taxonomy COInr_for_vtam_2022_05_06_dbV5\COInr_for_vtam_taxonomy.tsv -ncbitax_dir '' -blast_db COInr_for_vtam_2022_05_06_dbV5\COInr_for_vtam -outdir out -out_name test -ltg_params params\params_default.tsv -windows 1
~~~
The input file (-in) can also be a tsv file with headings. One of the columns must have 'sequence' as a heading.

### Algorithm

Table 1 Default parameter setting

| pid  | pcov | phit | taxn | seqn | refres  | ltgres  |
| ---- | ---- | ---- | ---- | ---- | ------- | ------- |
| 100  | 70   | 70   | 1    | 1    | species | species |
| 97   | 70   | 70   | 1    | 1    | species | species |
| 95   | 70   | 70   | 2    | 2    | species | species |
| 90   | 70   | 70   | 3    | 3    | genus   | species |
| 85   | 70   | 70   | 4    | 4    | family  | genus   |
| 80   | 70   | 70   | 4    | 4    | family  | genus   |

Query sequences are BLASTed against the reference database. Hits are then validated if 
 - % of identity is >= **pid**
 - %of coverage is >=**pcov**
 - the resolution of the reference sequence is at least **refres**

The Lowest Taxonomic Group (**LTG**) is determined if among the validated hits, there are at least **taxn** different taxa and **seqn** sequences. LTG is the Lowest Taxonomic Group  that contains at least **phit** proportion of the validated hits.

mkLTG starts at the highest **pid** value of the parameter setting. If LTG cannot be established using this **pid** and the associated values of the other parameters (e.g. not enough valid hits) the procedure is reiterated using decreasing **pid** values till LTG is established, or no more **pid** is provided in the parameter file.

### Arguments and options

 - **in**: Name of the input fasta or tsv file containing the sequences to be assigned. The tsv file must have 'sequence' as a heading for one of the columns; Compulsory
 - **blast_db**: name of the blast database;  Compulsory
 - **taxonomy**: TSV file with the following tab separated columns: taxid, parent_taxid, taxlevel, taxname, merged_taxid, taxlevel_index;  Either taxonomy or ncbitax_dir is compulsory
 - **ncbitax_dir**: Directory of ncbi taxonomy dmp files (downloaded from https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/) Only necessary if no **taxonomy** file is provided. The taxonomy file is automatically created from the dmp files and can be used if all sequences have NCBI taxIDs;  Either taxonomy or ncbitax_dir is compulsory
 - **ltg_params**: tsv file wih the following tab separated columns: pid pcov phit taxn seqn refres ltgres ;  Compulsory
 - **outdir**: name of the otput directory;  Compulsory
 - **out_name**: alpha-numeric string for naming output files;  Optional
 - **blastout**: If empty, run blast; otherwise this file is used to make LTG ;  Optional
 - **delete_tmp**: 0/1; if 1 delete temporary files after the run ;  Optional, 1 by default
 - **windows**: 1/0 set to one if running on windows;  Optional, 0 by default

BLAST parameters
 - **task**: megablast/blastn ;  Optional, megablast by default
 - **blast_e**: E-value threshold;  Optional, 1e-20 by default
 - **dust**: yes/no; apply dust filters for blast;  Optional, yes by default
 - **max_target_seqs**: maximum number of blast hits per query ;  Optional, 500 by default
 - **num_threads**: number of threads when running blast ;  Optional, 8 by default
 - **qcov_hsp_perc**: minium query coverage in blast ;  Optional, 70 by default
 - **batch_size**: batch size for blast (to avoid too large output files) ;  Optional, 1000 by default

   ​



## Optimize parameters

### Make Test Files

Randomly select 1000 sequences (seqn) from the COInr.tsv (db_tsv) file (downloaded from https://zenodo.org/record/6555985). Select only sequences assigned to species (select_taxlevel).

The randomy selected sequences (random_seq.fas) will be used as queries for the evaluation of the parameter settings, and the remaining sequences are formatted for blast (reduced_DB.fas). The random_seq.fas is BLASTed against the reduced_DB.fas database.

~~~~
perl scripts/make_test_files.pl -db_tsv COInr/COInr.tsv -taxonomy COInr/taxonomy.tsv -seqn 1000 -select_taxlevel species -outdir optimize/series1
~~~~

### Make parameter files using variable pid

Make parameter files for a large number of parameter value combinations. The possibe parameter values can be set by editing lines 16-22  in the script. This script produces parameter files with multiple pid values. To reduce the number of parameter files, pcov and phits are fixed within each parameter settings (same values for different pids). seqn is the same as taxn, and ltgres is a resolution higher than refres.

~~~~
perl scripts/make_param_files_pid_variable.pl -outdir optimize/param_files_pid_variable
~~~~

###  Make parameter files using a singe pid per parameter setting

Make parameter files for a large number of parameter value combinations. The possibe parameter values can be set by editing lines 15-19  in the script. This script produces parameter files with a single pid value. To reduce the number of parameter files,  seqn is the same as taxn, and ltgres is a resolution higher than refres.

~~~~
perl scripts/make_param_files_pid_fix.pl -outdir optimize/param_files_pid_fix
~~~~

### Run mkLTG for all parameter files 

Run mkLTG for each of the above produced parameter settings using random_seq.fas as a query and the results of the BLAST against reduced_DB.fas. 

At each major taxonomic level, count the number of sequences correctly assigned (TP), incorrectly assigned (FP) and non-assigned (FN)

#### Variable pid 
~~~~
perl scripts/run_ltg_with_multiple_param_setting.pl -fas optimize/series1/random_seq.fas -taxonomy COInr/taxonomy.tsv -blastout optimize/series1/blastout_random_seq.tsv -outdir optimize/series1/LTG_pid_variable -ltg_params_dir optimize/param_files_pid_variable
~~~~

#### Fix pid
~~~~
perl scripts/run_ltg_with_multiple_param_setting.pl -fas optimize/series1/random_seq.fas -taxonomy COInr/taxonomy.tsv -blastout optimize/series1/blastout_random_seq.tsv -outdir optimize/series1/LTG_pid_fix -ltg_params_dir optimize/param_files_pid_fix
~~~~

### Calculate F1 score, get the parameter setting with the highest F1, make graphs

Adapt the **compare_params.R** script to your use



## Versions

**mkLTG-0.1.0** 

2023-02-24

**mkLTG-0.1.1**

2023-03-29

Works with more sequences than the batch size

**mkLTG-0.1.2**

2023-04-11

Input sequences can be provided in fasta or TSV format

