use crate::bed::BedRecord;
use crate::codon::*;
use std::cmp::{max, min};
use std::collections::HashMap;
use std::fmt::Write;

pub fn build_gene_line(
    gene: &String,
    record: &BedRecord,
    coords: &HashMap<String, (u32, u32)>,
) -> String {
    assert!(gene.len() > 0);

    let (tx_start, tx_end) = coords.get(gene).unwrap();

    let gene_line = format!(
        "{}\t{}\tgene\t{}\t{}\t.\t{}\t.\tgene_id \"{}\";",
        record.chrom,
        "bed2gtf",
        tx_start + 1,
        tx_end,
        record.strand,
        gene
    );
    gene_line
}

pub fn build_gtf_line(
    record: &BedRecord,
    gene: &String,
    gene_type: &str,
    exon_start: u32,
    exon_end: u32,
    frame: u32,
    exon: i16,
    result: &mut Vec<(String, String, u32, u32, String, String, String)>,
) {
    assert!(record.tx_start < record.tx_end);

    let phase = match frame {
        0 => "0",
        1 => "2",
        2 => "1",
        _ => ".",
    };

    let mut attr = format!("gene_id \"{}\"; transcript_id \"{}\";", gene, record.name);

    if exon >= 0 {
        let (exon_id, nexon) = if record.strand == "+" {
            let exon_id = exon + 1;
            (exon_id as u16, exon + 1)
        } else {
            let exon_id = record.exon_count - exon as u16;
            (exon_id, exon_id as i16)
        };

        write!(
            attr,
            " exon_number \"{}\"; exon_id \"{}.{}\";",
            nexon, record.name, exon_id
        )
        .expect("Failed to write exon information");
    }

    result.push((
        record.chrom.clone(),
        gene_type.to_string(),
        exon_start + 1,
        exon_end,
        record.strand.clone(),
        phase.to_string(),
        attr,
    ));
}

pub fn write_features(
    i: usize,
    record: &BedRecord,
    gene: &String,
    // first_utr_end: u32,
    cds_start: u32,
    cds_end: u32,
    // last_utr_start: u32,
    frame: u32,
    result: &mut Vec<(String, String, u32, u32, String, String, String)>,
) {
    let exon_start = record.exon_start[i];
    let exon_end = record.exon_end[i];

    // if exon_start < first_utr_end {
    //     let end = min(exon_end, first_utr_end);
    //     let utr_type = if record.strand == "+" {
    //         "five_prime_utr"
    //     } else {
    //         "three_prime_utr"
    //     };
    //     build_gtf_line(record, gene, utr_type, exon_start, end, frame, -1, result);
    // }

    if record.cds_start < exon_end && exon_start < record.cds_end {
        let start = max(exon_start, cds_start);
        let end = min(exon_end, cds_end);
        if start < end {
            build_gtf_line(record, gene, "CDS", start, end, frame, i as i16, result);
        }
    }

    // if exon_end > last_utr_start {
    //     let start = max(exon_start, last_utr_start);
    //     let utr_type = if record.strand == "+" {
    //         "three_prime_utr"
    //     } else {
    //         "five_prime_utr"
    //     };
    //     build_gtf_line(record, gene, utr_type, start, exon_end, frame, -1, result);
    // }
}

pub fn write_codon(
    record: &BedRecord,
    gene: &String,
    gene_type: &str,
    codon: Codon,
    result: &mut Vec<(String, String, u32, u32, String, String, String)>,
) {
    build_gtf_line(
        record,
        gene,
        gene_type,
        codon.start,
        codon.end,
        0,
        codon.index as i16,
        result,
    );

    if codon.start2 < codon.end2 {
        build_gtf_line(
            record,
            gene,
            gene_type,
            codon.start,
            codon.end,
            codon.start2,
            (codon.end - codon.start) as i16,
            result,
        );
    }
}
