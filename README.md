![Crates.io](https://img.shields.io/crates/v/bed2gtf?color=green)
![GitHub](https://img.shields.io/github/license/alejandrogzi/bed2gtf?color=blue)

# <span style="color:green">**bed2gtf**</span>
A fast bed-to-gtf converter written in Rust.

translates
```
chr27 17266469 17281218 ENST00000541931.8 1000 + 17266469 17281218 0,0,200 2 103,74, 0,14675,
```
into
```
chr27 bed2gtf gene 17266470 17281218 . + . gene_id "ENSG00000151743"; gene_biotype "protein_coding";

chr27 bed2gtf transcript 17266470 17281218 . + . gene_id "ENSG00000151743"; transcript_id "ENST00000541931.8";

chr27 bed2gtf exon 17266470 17266572 . + . gene_id "ENSG00000151743"; transcript_id "ENST00000541931.8"; exon_number "1"; exon_id "ENST00000541931.8.1";

...
```

in a few seconds.

>[!IMPORTANT]
>
>Now bed2gtf uses a lexicograph-based algorithm to offer the user not only a .gtf file but a nicely sorted .gtf file. The algorithm was originally implemented in [gtfsort](https://github.com/alejandrogzi/gtfsort), and some parts have been coupled with bed2gtf code. 


## Usage
``` rust
Usage: bed2gtf[EXE] --bed <BED> --isoforms <ISOFORMS> --output <OUTPUT>

Arguments:
    --bed <BED>: a .bed file
    --isoforms <ISOFORMS>: a tab-delimited file
    --output <OUTPUT>: path to output file

Options:
    --help: print help
    --version: print version
```

>**Warning** 
>
>All the transcripts in .bed file should appear in the isoforms file.
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


## Library
to include bed2gtf as a library and use it within your project follow these steps:
1. include `bed2gtf = 1.8.0` under `[dependencies]` in the `Cargo.toml` file
2. the library name is `bed2gtf`, to use it just write:

    ``` rust
    use bed2gtf::bed2gtf; 
    ```
    or 
    ``` rust
    use bed2gtf::*;
    ```
3. invoke
    ``` rust
    let gtf = bed2gtf(bed: &String, isoforms: &String, output: &String)
    ```

## Build
to build bed2gtf from this repo, do:

1. get rust (as described above)
2. run `git clone https://github.com/alejandrogzi/bed2gtf.git && cd bed2gtf`
3. run `cargo run --release <BED> <ISOFORMS> <OUTPUT>`(arguments are positional, so you do not need to specify --bed/--isoforms)


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

### Annex

As part of this project, I developed a similar approach in python to benchmark computation times. I decided to include that script in this repository for people that maybe feel more comfortable using python over Rust. The usage is practically the same:

```
./bed2gtf.py -bed <BED> -gtf <OUTPUT> -iso <ISOFORMS>
```

## References

1. http://hgdownload.soe.ucsc.edu/admin/exe/
2. https://github.com/pfurio/bed2gtf
3. https://rdrr.io/github/wyguo/RTDBox/src/R/gtf2bed.R
4. https://github.com/MikeDacre/mike_tools/blob/master/bin/bed2gtf.py
5. https://www.ensembl.org/info/website/upload/gff.html
