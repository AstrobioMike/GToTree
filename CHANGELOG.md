# Change Log

## v2.0.0 (DD-MMM-2026)

### MAJOR VERSION INCREASE!
GToTree v2 is a complete rewrite from the ground up. It used to be an ~5,000-line bash runner 
with a large collection of standalone bash and python scripts. Now it's a python package. 

Why did I do this? I don't know. It worked fine, and I've successfully maintained it for almost 7 years and
700 citations the way it was. But here we are, and I think it's much better now :) 

There will likely be some snags I missed that will need to be worked out. If you hit any problems 
or find yourself struggling to find something you used to use, please post an issue, and I will do 
my best to rapidly address anything!

### Added

#### Reference genomes by taxonomy (`-w`)!

- the new `-w`/`--wanted-ref-tax` parameter takes a taxon name at any rank (e.g. `-w Bacteria`,
  `-w Nitrospirota`, `-w "Escherichia coli"`) and adds that taxon's reference genomes to
  the run, so gathering an accession list up front is no longer necessary :)
  - may be given more than once to pool several taxa (e.g. `-w Bacteria -w Archaea`);
    each is resolved and dereplicated on its own, then merged
  - `--source` specifies which taxonomy the genomes are drawn from, `gtdb` (default) or `ncbi`
  - `--derep-rank` dereplicates the selection down to one genome per unique value of a
    rank (default `auto`, which uses two ranks finer than the target; `off` keeps them all)
    - this is particularly helpful when trying to make broad-level trees where you don't need 
      a bajillion of the same species
  - `--exclusion-list` allows you to specify NCBI accessions you definitely don't want included
    regardless of the `-w`/`--wanted-ref-tax` search (it does not affect accessions you explicitly
    specified in the file passed to the `-a` parameter)


#### More pre-packaged single-copy gene-sets!

- 46 pre-packaged SCG-sets now ship, newly built from GTDB r232
  - run `gtt hmms` to see them all
  - the original NCBI-taxonomy-based sets are still available via manual download from here


#### Make your own SCG-sets

- `gtt gen-scg-hmms` is available to produce your own new single-copy gene HMMs for any group of input genomes,
  that can then be passed to GToTree for treeing


#### SCG-set auto-selection

- when `-w | --wanted-ref-tax` is used and a single-copy gene-set isn't specified with `-H`, 
  an appropriate pre-packaged SCG-set is now selected automatically
  - the lowest-rank pre-built set that covers all the requested taxa will be used
    - e.g., if given multiple wanted reference taxa including some bacteria and archaea, 
      the prebuilt `Bacteria-and-Archaea` HMM will be set; if eukaryotes were passed to `-w`, 
      the `Universal-Hug-et-al` set will be used
  - the set chosen and the reason for it are reported in the run's opening banner
  - passing `-H` explicitly always overrides this


#### Resuming interrupted runs

- the new `-R`/`--resume` flag attempts to continue an interrupted run
  - this is available on `gtt gen-scg-hmms` and `gtt search-annotations` too


#### New helper subcommands

- All prior helper commands are now grouped under `gtt` as subcommands. Run `gtt` by itself to
  see an overview. Main ones include:
  - `gtt hmms` - list all available pre-packaged SCG-HMMs
  - `gtt gen-scg-hmms` — build a new SCG-HMM set from any set of input genomes
  - `gtt update-headers` - modify the headers/labels in a previously complete GToTree run to add taxonomy info or labels from a mapping file
    - takes the output directory of a completed GToTree run and writes out a 
      new tree, alignment, and `genomes-summary-info.tsv` carrying updated genome labels
    - can accept `-m`, `--add-gtdb-tax`, `--add-ncbi-tax`, and `-lineage-ranks` just like the main program
    - puts new outputs in a new directory, it does not overwrite the originals
    - creates a "label-map.tsv" tracking changed labels
    - regenerates any iToL files for any Pfam/KO searches the original run did
  - `gtt get-accs-from-gtdb` / `gtt get-accs-from-ncbi` - pull accessions based on taxonomy if you want them (though the main `gtotree` program does this for you now)
  - `gtt search-annotations` — search input genomes for target Pfams (`-p`) and/or KOs (`-K`)
  - `gtt itol` - various subcommands for generating iToL-compatible files that can be drag-and-dropped onto a tree in the iToL web interface
  - `gtt data get` / `gtt data locations` — download or update the databases GToTree uses,
    and check or set the data-location environment variables (this is all handled automatically anyway)
  - `gtt test` — a quick end-to-end test of the installed environment


#### New outputs

- auto-creates iToL-compatible files for different input sources
  - so you can, e.g., drag-and-drop an iToL file right onto the tree in iToL to highlight all your MAGs that came in as nucleotide fasta files alongside references


### Parameter changes or additions

#### Defaults

- `-j`/`--num-jobs` now defaults to 4 (was 1)
- default output directory is now `gtotree-output`
- `-L`/`--lineage-ranks` now defaults to `domain,phylum,class,genus,species` (strain is
  no longer included)


#### Filtering

- added a gene-representation cutoff (`-r`/`--gene-rep-cutoff`, default 0.1) parameter
  - this sets the minimum proportion of genomes that must have hits to a target gene 
    for that gene to be retained
  - it's calculated on the genomes remaining after initial genome processing but before 
    genome-level filtering
- the genome-hits cutoff (`-G`/`--genome-hits-cutoff`, default 0.5) is now calculated
  against the number of target genes remaining/that will actually contribute to the final tree,
  rather than the starting total
  - e.g., starting from 100 target genes with 20 filtered out for too few hits, the
    percentage is now taken against 80 rather than 100. Meaning, a genome needs to have 
    >= 40 genes in order to meet the 0.5 default, rather than 50


#### Output files

- per-source summary tables (accessions, fasta files, genbank files, amino-acid files)
  have been consolidated into one `genomes-summary-info.tsv` with a `source` column
- every genome dropped at any stage is now recorded in a single `removed-genomes.tsv`
  (`genome_id`, `input`, `source`, `stage_removed`, `reason_removed`) rather than in
  scattered per-stage files
- output files and directories that used underscores now use dashes (sorry, i've grown to hate underscores when they aren't necessary)
- estimated percent completion and redundancy are no longer reported
  - that was basically "free" information before, but tools like CheckM2 do this way better, so I'd rather not spit out naive numbers anymore even if they're free

#### Behavior

- taxids are no longer pulled from input GenBank files; taxonomy can only be added for those pulled in 
  via `-w | --wanted-ref-tax` or provided as input accessions (`-a`)
  - this was always somewhat hacky given the variability of GenBank files. If this was
    actually useful to anyone, open an issue and I'll look into adding it back in
- all reference data is now fetched over http rather than ftp, which helps users behind firewalls that might block ftp
- reference datasets are drastically slimmed down and hosted as GitHub release assets, greatly speeding up any needed downloads
- console output has more progress bars and fewer walls of text, and more steps are logged


---


## v1.8.19 (16-Jul-2026)

- removed high-genome notice that *required* user interaction to proceed (https://github.com/AstrobioMike/GToTree/issues/83#issuecomment-4989358508)
  - there is already a for the aligner used if there are a lot of genomes, and if they shut that off explicitly, they already saw my notice and suggestion
- added support for bailing on slow connections and restarting downloads for the gtdb/ncbi data tables (since they are hosted on github now and this happens from time-to-time)

--- 


## v1.8.18 (4-Jul-2026)

- robustness improvements to assembly downloading
- downloads pre-packaged, slimmed ncbi assembly summary and gtdb metadata tables now (much faster)
  - --use-ecogenomics flag dropped from gtt-get-accessions-from-GTDB as it's not longer relevant

--- 

## v1.8.17 (11-May-2026)

### Fixed
- minor updates to `gtt-gen-SCG-HMMs` for things that may have caused it to fail sometimes, re-uses the main GToTree-stored NCBI assembly data now, and also unlocked the PFam version again since the interpro folks added back in the data we needed (see https://github.com/AstrobioMike/GToTree/issues/104)

---

## v1.8.16 (09-Jul-2025)

### Fixed
- Fix to gtt-get-accessions-from-GTDB erroring out after initial use

---

## v1.8.15 (27-Jun-2025)

### Added
- added a flag to `gtt-get-accessions-from-GTDB` (`--use-ecogenomics`) to allow specifing to download from [data.gtdb.ecogenomic.org/releases](https://data.gtdb.ecogenomic.org/releases/) instead of the default [data.ace.uq.edu.au/public/gtdb/data/releases](https://data.ace.uq.edu.au/public/gtdb/data/releases/) thanks to the suggestion from @Stian-2rz (https://github.com/AstrobioMike/GToTree/issues/107)
- a partitions file in nexus format (`<outdir>/run_files/Partitions.nex`) is produced in addition to the regular text-formatted one (`<outdir>/run_files/Partitions.txt`), because for some iqtree model settings the txt format has leads to a bug – also thanks to @Stian-2rz! (https://github.com/AstrobioMike/GToTree/issues/108)

### Changed
- currently skipping certificate checking on GTDB downloads (prompted from https://github.com/AstrobioMike/GToTree/issues/107)
  - this impacts `gtt-get-accessions-from-GTDB` and the internally used `gtt-check-or-setup-GTDB-files`

---

## v1.8.14 (21-Apr-2025)

### Changed
- change to taxonkit call when adding NCBI tax info (now using `reformat2` and a pattern) in order to deal with NCBI tax-structure update

---

## v1.8.13 (18-Mar-2025)

### Changed
- changed `gtt-gen-SCG-HMMs` to only use Pfam 37.0 for now (as later versions don't have one of the required files currently; see https://github.com/AstrobioMike/GToTree/issues/104)

---

## v1.8.12 (11-Mar-2025)

### Changed
- changed GTDB download links from https://data.gtdb.ecogenomic.org/releases/ to https://data.ace.uq.edu.au/public/gtdb/data/releases/ due to the former becoming prohibitively slow recently

---

## v1.8.11 (10-Mar-2025)

### Added
- VeryFastTree is now an available treeing program (`-T`)

### Changed
- when using `gtt-get-accessions-from-GTDB`, if the requested taxon has spaces in it (e.g., `gtt-get-accessions-from-GTDB -t "Bacillus_A anthracis"`), the output files will have spaces replaced with dashes now
  - e.g., one of the outputs will now be "GTDB-Bacillus_A-anthracis-species-accs.txt" instead of "GTDB-Bacillus_A anthracis-species-accs.txt"

---

## v1.8.10 (3-Feb-2025)

### Added
- saving ncbi downloaded files is possible when debug flag (`-d`) is set as requested in https://github.com/AstrobioMike/GToTree/issues/95, implemented in https://github.com/AstrobioMike/GToTree/pull/102
  - with the debug flag set while running, it will keep specific files in `<output_dir>/<tmp_dir>/ncbi-downloads/`:
    - if amino-acid seqs are used, it will keep the downloaded amino-acid seqs
    - if there were no amino-acid seqs, and the genome had to be downloaded, it will keep the downloaded genome and the prodigal-called amino-acid seqs
    - if using nucleotide mode (`-z`), it will keep the downloaded genome and the prodigal-called nt cds and amino-acid seqs

---

## v1.8.9 (31-Jan-2025)

### Fixed
- added logic to catch, exit, and report when muscle doesn't successfully produce an alignment for a single-copy gene-set (thanks to https://github.com/AstrobioMike/GToTree/issues/101)

---

## v1.8.8 (7-Oct-2024)

### Changed
- updated the call to FastTree and FastTreeMP to be include -nt and -gtr when GToTree is run in nucleotide mode (-z)

### Fixed
- properly saving additional pfam target HMMs when that functionality is used

---

## v1.8.7 (29-Sep-2024)

### Added
- `gtt-gen-SCG-HMMs` now reports which version of PFAM was used (prints it out to the terminal and writes it to a file)

### Changed
- improvements to the "Universal" Hug et al. gene set thanks so much to @molly-kholodova digging in and reaching out!
    - PF00181 ("Ribosomal_L2") was changed to PF03947 ("Ribosomal_L2_C")
        - the C-terminal (which PF03947 covers) is better conserved
    - PF00827 ("Ribosomal_L15") was changed to PF00828 ("Ribosomal_L27A")
        - PF00827 was archaea/euk only, PF00828 holds the bac/arc L15 also
    - PF17135 ("Ribosomal_L18") was changed to PF00861 ("Ribosomal_L18p")
        - the PF00861 model is better distributed

---

## v1.8.6 (8-May-2024)

### Fixed
- fixed when taxonomy information wasn't being added to labels when running in nucleotide mode (`-z`; https://github.com/AstrobioMike/GToTree/issues/91)

---

## v1.8.5 (1-May-2024)

### Changed
- update to `gtt-gen-SCG-HMMs` to deal with ncbi assembly summary files having a column name of "#assembly_accession" instead of what was once "# assembly_accession"

---

## v1.8.4 (28-Nov-2023)

### Fixed
- fixed an issue that prevented moving forward when there were more than 12,500 input genomes (https://github.com/AstrobioMike/GToTree/issues/83)

---

## v1.8.3 (14-Oct-2023)

### Changed
- updated links to GTDB files as they switched from .tar.gz extensions to .tsv.gz extensions in latest release, thanks to note from @jmtsuji (https://github.com/AstrobioMike/GToTree/issues/81)

---

## v1.8.2 (26-Jul-2023)

### Added
- added http option to gtt-test.sh (`gtt-test.sh http`) thanks to https://github.com/AstrobioMike/GToTree/issues/78 (https://github.com/AstrobioMike/GToTree/commit/9eb248ad5a54563370978d3575727eb63ad93483)

### Fixed
- updated `gtt-get-ncbi-tax-data` to appropriately pull from http instead of ftp also thanks to https://github.com/AstrobioMike/GToTree/issues/78
- fix to check for ncbi assemblies "date-retrieved.txt" file, as also caught and fixed by @hyphaltip (https://github.com/AstrobioMike/GToTree/pull/80) 🙏 

---

Earlier version changes are tracked on the [releases page](https://github.com/AstrobioMike/GToTree/releases).

---
