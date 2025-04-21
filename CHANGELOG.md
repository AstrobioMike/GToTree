# Change Log

## v1.8.14 (21-Apr-2025)

### Changed
- change to taxonkit call when adding NCBI tax info (now using `reformat2` and a pattern) in order to deal with NCBI tax-structure update


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
