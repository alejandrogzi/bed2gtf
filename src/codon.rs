use crate::bed::BedRecord;
use std::cmp::{max, min};

#[derive(Debug, Clone)]
pub struct Codon {
    pub start: u32,
    pub end: u32,
    pub index: u32,
    pub start2: u32,
    pub end2: u32,
}

impl Codon {
    pub fn new() -> Codon {
        Codon {
            start: 0,
            end: 0,
            index: 0,
            start2: 0,
            end2: 0,
        }
    }
}

pub fn first_codon(record: &BedRecord) -> Option<Codon> {
    let exon_frames = record.get_frames();
    record
        .exon_start
        .iter()
        .zip(record.exon_end.iter())
        .enumerate()
        .find_map(|(mut index, (&start, &end))| {
            let frame = exon_frames.get(index)?;
            let mut codon = Codon::new();

            if *frame < 0 {
                return Some(codon);
            }

            let cds_start = max(start, record.cds_start);
            let cds_end = min(end, record.cds_end);

            let frame = if record.strand == "+" {
                *frame
            } else {
                (*frame + (cds_end - cds_start) as i16) % 3
            };

            if frame == 0 {
                codon.start = cds_start;
                codon.end = cds_start + 3;
                codon.index = index as u32;
                let diff = cds_end - cds_start;

                if diff >= 3 {
                    Some(codon)
                } else {
                    index += 1;
                    if index >= exon_frames.len() {
                        Some(codon)
                    } else {
                        let need = 3 - diff;
                        if diff < need {
                            Some(codon)
                        } else {
                            codon.start2 = cds_start;
                            codon.end2 = cds_start + need;
                            Some(codon)
                        }
                    }
                }
            } else {
                Some(Codon::new())
            }
        })
}

pub fn last_codon(record: &BedRecord) -> Option<Codon> {
    let exon_frames = record.get_frames();
    record
        .exon_start
        .iter()
        .zip(record.exon_end.iter())
        .enumerate()
        .rev() // Reverse the iterator to start from the last exon
        .find_map(|(mut index, (&start, &end))| {
            let mut codon = Codon::new();
            let frame = exon_frames.get(index)?;
            let cds_start = max(start, record.cds_start);
            let cds_end = min(end, record.cds_end);

            let frame = if record.strand == "+" {
                (*frame + (cds_end - cds_start) as i16) % 3
            } else {
                *frame
            };

            if frame == 0 {
                codon.start = max(cds_start, cds_end - 3); // Find the last 3 bases of the CDS
                codon.end = cds_end;
                codon.index = index as u32;
                let diff = cds_end - cds_start;

                if diff >= 3 {
                    Some(codon)
                } else {
                    index += 1;
                    if index >= exon_frames.len() {
                        Some(codon)
                    } else {
                        let need = 3 - diff;
                        if diff < need {
                            Some(codon)
                        } else {
                            codon.start2 = cds_start;
                            codon.end2 = cds_start + need;
                            Some(codon)
                        }
                    }
                }
            } else {
                Some(Codon::new())
            }
        })
}

pub fn codon_complete(codon: &Codon) -> bool {
    ((codon.end - codon.start) + (codon.end2 - codon.start2)) == 3
}
