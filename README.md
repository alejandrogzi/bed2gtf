# bed2gtf
A fast bed-to-gtf converter written in Rust.

translates
```
chr27 17266469 17281218 ENST00000541931.8 1000 + 17266469 17281218 0,0,200 2 103,74, 0,14675,
```
into
```
chr27 stdin gene 17266470 17281218 . + . gene_id "ENSG00000151743"; gene_biotype "protein_coding";

chr27 stdin transcript 17266470 17281218 . + . gene_id "ENSG00000151743"; transcript_id "ENST00000541931.8";

chr27 stdin exon 17266470 17266572 . + . gene_id "ENSG00000151743"; transcript_id "ENST00000541931.8"; exon_number "1"; exon_id "ENST00000541931.8.1";

...
```

in a few seconds.

## Usage
``` rust
Usage: bed2gtf[EXE] --bed <BED> --isoforms <ISOFORMS> --verbose [verbose]

Arguments:
    --bed <BED>: a .bed file
    --isoforms <ISOFORMS>: a tab-delimited file
    --verbose [verbose]: print verbose [default: true]

Options:
    --help: print help
    --version: print version
```


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
    chr20   50222035    50222038    ENST00000595977.1735    ...
    ```

    see [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) for more information

2. a tab-delimited .txt/.tsv/.csv/... file with genes/isoforms:

    ```
    > cat isoforms.txt

    ENSG00000198888 ENST00000361390
    ENSG00000198763 ENST00000361453
    ENSG00000198804 ENST00000361624
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
1. include `bed2gtf = 1.1.0` under `[dependencies]` in the `Cargo.toml` file
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
    let gtf = bed2gtf(bed: PathBuf, isoforms: PathBuf)
    ```

## Build
to build bed2gtf from this repo, do:

1. get rust (as described above)
2. run `git clone https://github.com/alejandrogzi/bed2gtf.git && cd bed2gtf`
3. run `cargo run --release <BED> <ISOFORMS>`(arguments are positional, so you do not need to specify --bed/--isoforms)


## Output

bed2gtf will send the output directly to the same .bed file path

```
.
├── ...
├── isoforms.txt
├── annotation.bed
└── output.gtf
```
where `output.gtf` is the result. Any intermediate files are deleted, since there is no need to keep them. 

## FAQ
### Why?
UCSC offers a fast way to convert BED into GTF files through KentUtils or specific binaries (1) + several other bioinformaticians have shared scripts trying to replicate a similar solution (2,3,4).

A GTF file is a 9-column tab-delimited file that holds gene annotation data for a specific assembly (5). The 9th column defines the attributes of each entry. This field is important, as some post-processing tools that handle GTF files need them to extract gene information (e.g. STAR, arriba, etc). An incomplete GTF attribute field would probably lead to annotation-related errors in these software. 

Of the available tools/scripts mentioned above, none produce a fully functional attribute GTF file conversion. (1) uses a two-step approach (bedToGenePred | genePredToGtf) written in C, which is extremely fast. Since a .bed file does not preserve any gene-related information, this approach fails to a) include correct gene_id attributes (duplicated transcript_ids) b) append 3rd column gene features.

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

This is where bed2gtf comes in: a fast and memory efficient BED-to-GTF converter written in Rust. In ~15 seconds this tool produces a fully functional GTF converted file with all the needed features needed for post-processing tools. 

### How?
bed2gtf takes advantage of UCSC binaries. These executables show slightly lower computation times than those seen in a similar rust implementation. This tool uses bedToGenePred + genePredToGtf as a first step, automates their download and updates permissions. After that, this tool overcomes the limitations of gene_id and gene-feature fields by expecting a tab-delimited isoforms file. Each unique transcript seen in the new gtf file is mapped to their respective gene_id. With this information, bed2gtf simply constructs gene-feature lines and corrects gene_ids.

### Limitations
At the time of bed2gtf being publicly available some gaps have not been covered yet. 

1. Sorting. This tool does not implement a sorting step. While it is not mandatory, and probably would not affect post-processing tools, it could increase computation time. A sorting step should consider a) chromosome, b) boundaries (start-end) and c) features (gene > transcript > exon, start_codon, CDS, stop_codon, UTR). 
2. Time. Although bed2gtf is pretty fast compared to other implementations (e.g python3; unpublished data), I personally think that this tool could achieve even faster computation times with more detailed algorithms and structures.
3. Biotype. As you may know (or not), GTF files specify the gene_biotype of each entry (e.g. protein_coding, processed_pseudogene, snoRNA, etc). This is probably the biggest limitation in this release. Currently, bed2gtf assumes all genes/transcripts provided are protein_coding. In future releases will probably be an option to specify the gene_biotype [-b/--biotype]. This feature could also be extracted from Ensembl BioMart and append it to the isoforms file; however, map them to each line will probably cost efficiency without a better framework.


## References

1. http://hgdownload.soe.ucsc.edu/admin/exe/
2. https://github.com/pfurio/bed2gtf
3. https://rdrr.io/github/wyguo/RTDBox/src/R/gtf2bed.R
4. https://github.com/MikeDacre/mike_tools/blob/master/bin/bed2gtf.py
5. https://www.ensembl.org/info/website/upload/gff.html
