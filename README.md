<p align="center">
  <h1 align="center">
  bed2gtf
  </h1>

  <p align="center">
    <a href="https://img.shields.io/badge/version-0.1.0dev-green" target="_blank">
      <img alt="Version Badge" src="https://img.shields.io/badge/version-0.2.3-green">
    </a>
    <a href="https://crates.io/crates/gxf2bed" target="_blank">
      <img alt="Crates.io Version" src="https://img.shields.io/crates/v/bed2gtf?color=green">
    </a>
    <a href="https://github.com/alejandrogzi/gxf2bed" target="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/alejandrogzi/bed2gtf?color=blue">
    </a>
    <a href="https://crates.io/crates/gxf2bed" target="_blank">
      <img alt="Crates.io Total Downloads" src="https://img.shields.io/crates/d/bed2gtf">
    </a>
    <a href="https://anaconda.org/bioconda/gxf2bed" target="_blank">
      <img alt="Conda Platform" src="https://img.shields.io/conda/pn/bioconda/bed2gtf">
    </a>
  </p>


  <p align="center">
    A high-performance bed-to-gtf converter written in Rust.

    // translates

    chr27 17266469 17281218 ENST00000541931.8 1000 + 17266469 17281218 0,0,200 2 103,74, 0,14675,

    // into

    chr27 bed2gtf gene 17266470 17285418 . + . gene_id "ENSG00000151743";

    chr27 bed2gtf transcript 17266470 17281218 . + . gene_id "ENSG00000151743"; transcript_id "ENST00000541931.8";

    chr27 bed2gtf exon 17266470 17266572 . + . gene_id "ENSG00000151743"; transcript_id "ENST00000541931.8"; exon_number "1"; exon_id "ENST00000541931.8.1";

    ...
  </p>

</p>


Converts
- *Homo sapiens* GRCh38 GENCODE 44 (252,835 transcripts) in 3.25 seconds.
- *Mus musculus* GRCm39 GENCODE 44 (149,547 transcritps) in 1.99 seconds.
- *Canis lupus familiaris* ROS_Cfam_1.0 Ensembl 110 (55,335 transcripts) in 1.20 seconds.
- *Gallus galus* bGalGal1 Ensembl 110 (72,689 transcripts) in 1.36 seconds.

> What's new on v.1.9.3
>
> - Fixes a bug with .gz decoder
> - Implements reading .bed.gz files!
> - Fixes bug described in issue #11 with versioning

## Usage
``` rust
Usage: bed2gtf[EXE] --bed/-b <BED> --isoforms/-i <ISOFORMS> --output/-o <OUTPUT>

Arguments:
    -b, --bed <BED>: a .bed file
    -i, --isoforms <ISOFORMS>: a tab-delimited file
    -o, --output <OUTPUT>: path to output file
    -g, --gz[=<FLAG>]          Compress output file [default: false] [possible values: true, false]
    -n, --no-gene[=<FLAG>]     Flag to disable gene_id feature [default: false] [possible values: true, false]

Options:
    --help: print help
    --version: print version
    --threads/-t: number of threads (default: max ncpus)
    --gz: compress output .gtf
```

> [!WARNING]
>
> All the transcripts in .bed file should appear in the isoforms file.


> [!TIP]
>
> Here are some commands to get you started:
> ```bash
> # convert a .bed file to .gtf (if you have an isoforms file [gene -> transcript names] and want gene_ids in the output .gtf)
> bed2gtf -b file.bed -i isoforms.txt -o file.gtf
>
> # convert a ,.bed file to .gtf without isoforms [same things as UCSC bedToGtf]
> bed2gtf -b file.bed -o file.gtf --no-gene
>
> # convert a .bed.gz file to a .gtf [with or without isoforms]
> bed2gtf -b file.bed.gz -i isoforms.txt -o file.gtf --gz
> bed2gtf -b file.bed.gz -o file.gtf --gz --no-gene
>
> # convert a .bed.gz to a .gtf.gz [with or without isoforms]
> bed2gtf -b file.bed.gz -i isoforms.txt -o file.gtf --gz
> bed2gtf -b file.bed.gz -o file.gtf --gz --no-gene
>

#### crate: [https://crates.io/crates/bed2gtf](https://crates.io/crates/bed2gtf)

<details>
<summary>click for detailed formats</summary>
<p>
bed2gtf just needs two files:

1. a .bed file

    tab-delimited files with 3 required and 9 optional fields:

    ```
    chrom   chromStart  chromEnd      name    ...
      |         |           |           |
    chr20   50222035    50222038    ENST00000595977    ...
    ```

    see [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) for more information

2. a tab-delimited .txt/.tsv/.csv/... file with genes/isoforms (all the transcripts in .bed file should appear in the isoforms file):

    ```
    > cat isoforms.txt

    ENSG00000198888 ENST00000361390
    ENSG00000198763 ENST00000361453
    ENSG00000198804 ENST00000361624
    ENSG00000188868 ENST00000595977
    ```

    you can build a custom file for your preferred species using [Ensembl BioMart](https://www.ensembl.org/biomart/martview).

</p>
</details>

## Installation
to install bed2gtf on your system follow this steps:
1. get rust: `curl https://sh.rustup.rs -sSf | sh` on unix, or go [here](https://www.rust-lang.org/tools/install) for other options
2. run `cargo install bed2gtf` (make sure `~/.cargo/bin` is in your `$PATH` before running it)
4. use `bed2gtf` with the required arguments
5. enjoy!

## Build
to build bed2gtf from this repo, do:

1. get rust (as described above)
2. run `git clone https://github.com/alejandrogzi/bed2gtf.git && cd bed2gtf`
3. run `cargo run --release -- -b <BED> -i <ISOFORMS> -o <OUTPUT>`

## Container image
to build the development container image:
1. run `git clone https://github.com/alejandrogzi/bed2gtf.git && cd bed2gtf`
2. initialize docker with `start docker` or `systemctl start docker`
3. build the image `docker image build --tag bed2gtf .`
4. run `docker run --rm -v "[dir_where_your_gtf_is]:/dir" bed2gtf -b /dir/<BED> -i /dir/<ISOFORMS> -o /dir/<OUTPUT>`

## Conda
to use bed2gtf through Conda just:
1. `conda install bed2gtf -c bioconda` or `conda create -n bed2gtf -c bioconda gtfsort`

## Output

bed2gtf will send the output directly to the same .bed file path if you specify so

```
bed2gtf annotation.bed isoforms.txt output.gtf

.
├── ...
├── isoforms.txt
├── annotation.bed
└── output.gtf
```
where `output.gtf` is the result.

## FAQ
### Why?
UCSC offers a fast way to convert BED into GTF files through KentUtils or specific binaries (1) + several other bioinformaticians have shared scripts trying to replicate a similar solution (2,3,4).

A GTF file is a 9-column tab-delimited file that holds gene annotation data for a specific assembly (5). The 9th column defines the attributes of each entry. This field is important, as some post-processing tools that handle GTF files need them to extract gene information (e.g. STAR, arriba, etc). An incomplete GTF attribute field would probably lead to annotation-related errors in these software.

Of the available tools/scripts mentioned above, none produce a fully functional attribute GTF file conversion. (1) uses a two-step approach (bedToGenePred | genePredToGtf) written in C, which is extremely fast. Since a .bed file does not preserve any gene-related information, this approach fails to a) include correct gene_id attributes (duplicated transcript_ids) if no refTable is included b) append 3rd column gene features.

This is an example:

```
chr27 stdin transcript 17266470 17281218 . + . gene_id "ENST00000541931.8"; transcript_id "ENST00000541931.8";

chr27 stdin exon 17266470 17266572 . + . gene_id "ENST00000541931.8"; transcript_id "ENST00000541931.8"; exon_number "1"; exon_id "ENST00000541931.8.1";
```


On the other hand, available scripts (2,3,4) fall into bad-formatted outputs unable to be used as input to other tools. Some of them show a very customed format, far from a complete GTF file (2):

```
chr20 ---- peak 50222035 50222038 . + . peak_id "chr20_50222035_50222038";

chr20 ---- peak 50188548 50189130 . + . peak_id "chr20_50188548_50189130";
```
and others (4) just provide exon-related information:

```
chr20 ensembl exon 50222035 50222038 . + . gene_id "ENST00000595977.1735"; transcript_id "ENST00000595977.1735"; exon_number "0

chr20 ensembl exon 50188548 50188930 . + . gene_id "ENST00000595977.3403"; transcript_id "ENST00000595977.3403"; exon_number "0
```

This is where bed2gtf comes in: a fast and memory efficient BED-to-GTF converter written in Rust. In ~4 seconds this tool produces a fully functional GTF converted file with all the needed features needed for post-processing tools.

### How?
bed2gtf is basically the reimplementation of C binaries merged in 1 step. This tool evaluates the position of k exons in j transcript, calculates start/stop/codon/UTR positions preserving reading frames and adjust the index + 1 (to be compatible with GTF convention). The isoforms file works as the refTable in C binaries to map each transcript to their respective gene; however, bed2gtf takes advantage of this and adds an additional "gene" line (to be compatible with other tools).

## References

1. http://hgdownload.soe.ucsc.edu/admin/exe/
2. https://github.com/pfurio/bed2gtf
3. https://rdrr.io/github/wyguo/RTDBox/src/R/gtf2bed.R
4. https://github.com/MikeDacre/mike_tools/blob/master/bin/bed2gtf.py
5. https://www.ensembl.org/info/website/upload/gff.html
