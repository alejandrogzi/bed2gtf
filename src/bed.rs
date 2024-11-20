use std::cmp::{max, min};

#[derive(Debug, PartialEq)]
pub struct BedRecord {
    pub chrom: String,
    pub tx_start: u32,
    pub tx_end: u32,
    pub name: String,
    pub strand: String,
    pub cds_start: u32,
    pub cds_end: u32,
    pub exon_count: u16,
    pub exon_start: Vec<u32>,
    pub exon_end: Vec<u32>,
}

impl BedRecord {
    pub fn parse(line: &str) -> Result<BedRecord, String> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            return Err(format!("Expected at least 12 fields, found {}", line));
        }

        let chrom = fields[0].to_string();
        let name = fields[3].to_string();
        let strand = fields[5].to_string();

        let get = |field: &str| field.parse::<u32>().map_err(|_| "Cannot parse field");
        let tx_start = get(fields[1])?;
        let tx_end = get(fields[2])?;
        let cds_start = get(fields[6])?;
        let cds_end = get(fields[7])?;
        let exon_count = get(fields[9])? as u16;

        let group = |field: &str| -> Result<Vec<u32>, &'static str> {
            field
                .split(',')
                .filter_map(|num| {
                    if !num.is_empty() {
                        Some(num.parse::<u32>().expect("Cannot parse number"))
                    } else {
                        None
                    }
                })
                .map(|num| Ok(num))
                .collect()
        };

        let exon_start = group(fields[11])?;
        let exon_end = group(fields[10])?;

        if exon_start.len() != exon_end.len() {
            return Err("Exon start and end vectors have different lengths".to_string());
        }

        let exon_starts: Vec<u32> = exon_start.iter().map(|&s| s + tx_start).collect();
        let exon_ends: Vec<u32> = exon_end
            .iter()
            .enumerate()
            .map(|(i, &s)| s + exon_starts[i])
            .collect();

        Ok(BedRecord {
            chrom: chrom.to_string(),
            tx_start: tx_start,
            tx_end: tx_end,
            name: name.to_string(),
            strand: strand.to_string(),
            cds_start: cds_start,
            cds_end: cds_end,
            exon_count: exon_count,
            exon_start: exon_starts,
            exon_end: exon_ends,
        })
    }

    pub fn get_frames(&self) -> Vec<i16> {
        let mut exon_frames: Vec<i16> = vec![0; self.exon_count as usize];
        let mut cds: u32 = 0;

        let exon_range = if self.strand == "+" {
            (0..(self.exon_count as usize)).collect::<Vec<_>>()
        } else {
            (0..(self.exon_count as usize)).rev().collect::<Vec<_>>()
        };

        for exon in exon_range {
            let cds_exon_start = max(self.exon_start[exon], self.cds_start);
            let cds_exon_end = min(self.exon_end[exon], self.cds_end);

            if cds_exon_start < cds_exon_end {
                exon_frames[exon] = (cds % 3) as i16;
                cds += cds_exon_end - cds_exon_start;
            } else {
                exon_frames[exon] = -1;
            }
        }

        exon_frames
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_record() {
        let line =
            "chr15\t81000922\t81005788\tENST00000267984\t0\t+\t81002271\t81003360\t0\t1\t4866,\t0,";
        let record = BedRecord::parse(line).unwrap();

        assert_eq!(record.chrom, "chr15");
        assert_eq!(record.tx_start, 81000922);
        assert_eq!(record.tx_end, 81005788);
        assert_eq!(record.name, "ENST00000267984");
        assert_eq!(record.strand, "+");
        assert_eq!(record.cds_start, 81002271);
        assert_eq!(record.cds_end, 81003360);
        assert_eq!(record.exon_count, 1);
        assert_eq!(record.exon_start, vec![81000922]);
        assert_eq!(record.exon_end, vec![81005788]);
    }

    #[test]
    fn get_exon_frames() {
        let line = "chr11\t13934505\t13958243\tENST00000674667\t1000\t-\t13934505\t13958243\t0,0,200\t9\t224,217,228,198,149,142,115,157,49,\t0,1305,2811,5576,10085,14837,18016,19498,23689,";
        let record = BedRecord::parse(line).unwrap();

        assert_eq!(record.get_frames(), vec![1, 0, 0, 0, 1, 0, 2, 1, 0]);
    }

    #[test]
    fn invalid_record() {
        let line =
            "chr15\t81000922\t81005788\tENST00000267984\t0\t+\t81002271\t81003360\t0\t1\t4866,";
        let record = BedRecord::parse(line);

        assert_eq!(
            record,
            Err(
                "Bed line has less than 12 fields and cannot be parsed into a BedRecord"
                    .to_string()
            )
        );
    }

    #[test]
    fn empty_record() {
        let line = "";
        let record = BedRecord::parse(line);

        assert_eq!(
            record,
            Err(
                "Bed line has less than 12 fields and cannot be parsed into a BedRecord"
                    .to_string()
            )
        );
    }
}
