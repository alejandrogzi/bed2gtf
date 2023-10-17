use crate::ParseError;
use colored::Colorize;
use std::cmp::{max, min};

#[derive(Clone, Debug, PartialEq)]
pub struct BedRecord {
    chrom: String,
    tx_start: i32,
    tx_end: i32,
    name: String,
    strand: String,
    cds_start: i32,
    cds_end: i32,
    exon_count: i16,
    exon_start: Vec<i32>,
    exon_end: Vec<i32>,
}

impl BedRecord {
    pub fn new(line: &str) -> Result<BedRecord, ParseError> {
        if line.is_empty() {
            eprintln!(
                "{}",
                "bed2gtf found an empty line! Check your input file.".bright_red()
            );
            return Err(ParseError::Empty);
        }

        let fields = splitb(line, "\t")?;

        if fields.len() < 12 {
            eprintln!(
                "{}",
                "bed2gtf found an invalid BED line! Check your input file.".bright_red()
            );
            return Err(ParseError::Invalid);
        }

        let chrom: String = fields[0].clone();
        let tx_start: i32 = fields[1].parse().unwrap();
        let tx_end: i32 = fields[2].parse().unwrap();
        let name: String = fields[3].clone();
        let strand: String = fields[5].clone();
        let cds_start: i32 = fields[6].parse().unwrap();
        let cds_end: i32 = fields[7].parse().unwrap();
        let exon_count: i16 = fields[9].parse().unwrap();
        let mut exon_start: Vec<i32> = splitb(fields[11].as_str(), ",")?
            .iter()
            .map(|x| x.parse().unwrap())
            .collect();
        let mut exon_end: Vec<i32> = splitb(fields[10].as_str(), ",")?
            .iter()
            .map(|x| x.parse().unwrap())
            .collect();

        for x in 0..exon_count as usize {
            exon_start[x] += tx_start;
            exon_end[x] += exon_start[x];
        }

        Ok(BedRecord {
            chrom,
            tx_start,
            tx_end,
            name,
            strand,
            cds_start,
            cds_end,
            exon_count,
            exon_start,
            exon_end,
        })
    }

    pub fn chrom(&self) -> &str {
        &self.chrom
    }

    pub fn tx_start(&self) -> i32 {
        self.tx_start
    }

    pub fn tx_end(&self) -> i32 {
        self.tx_end
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn strand(&self) -> &str {
        &self.strand
    }

    pub fn cds_start(&self) -> i32 {
        self.cds_start
    }

    pub fn cds_end(&self) -> i32 {
        self.cds_end
    }

    pub fn exon_count(&self) -> i16 {
        self.exon_count
    }

    pub fn exon_end(&self) -> &Vec<i32> {
        &self.exon_end
    }

    pub fn exon_start(&self) -> &Vec<i32> {
        &self.exon_start
    }

    pub fn layer(line: &str) -> Result<(String, i32, String), ParseError> {
        if line.is_empty() {
            eprintln!(
                "{}",
                "bed2gtf found an empty line! Check your input file.".bright_red()
            );
            return Err(ParseError::Empty);
        }

        let fields = splitb(line, "\t")?;

        if fields.len() < 12 {
            eprintln!(
                "{}",
                "bed2gtf found an invalid BED line! Check your input file.".bright_red()
            );
            return Err(ParseError::Invalid);
        }

        let chrom: String = fields[0].clone();
        let tx_start: i32 = fields[1].parse().unwrap();
        let line: String = fields.join("\t");

        Ok((chrom, tx_start, line))
    }

    pub fn get_exon_frames(&self) -> Vec<i32> {
        let mut exon_frames: Vec<i32> = vec![0; self.exon_count as usize];
        let mut cds: i32 = 0;
        let (start, end) = (0i16, self.exon_count);

        if self.strand == "+" {
            for exon in (start..end).step_by(1) {
                let cds_exon_start = max(self.exon_start[exon as usize], self.cds_start);
                let cds_exon_end = min(self.exon_end[exon as usize], self.cds_end);

                if cds_exon_start < cds_exon_end {
                    exon_frames[exon as usize] = cds % 3;
                    cds += cds_exon_end - cds_exon_start;
                } else {
                    exon_frames[exon as usize] = -1;
                }
            }
            exon_frames
        } else {
            for exon in (start..end).step_by(1).rev() {
                let cds_exon_start = max(self.exon_start[exon as usize], self.cds_start);
                let cds_exon_end = min(self.exon_end[exon as usize], self.cds_end);

                if cds_exon_start < cds_exon_end {
                    exon_frames[exon as usize] = cds % 3;
                    cds += cds_exon_end - cds_exon_start;
                } else {
                    exon_frames[exon as usize] = -1;
                }
            }
            exon_frames
        }
    }
}

pub fn splitb(line: &str, sep: &str) -> Result<Vec<String>, ParseError> {
    let bytes = line.as_bytes().iter().enumerate();
    let mut start = 0;
    let mut entries = Vec::new();

    for (i, byte) in bytes {
        if *byte == sep.as_bytes()[0] {
            let word = line[start..i].to_string();
            if !word.is_empty() {
                entries.push(word);
            }
            start = i + 1;
        }
    }

    let last = line[start..].to_string();
    if !last.is_empty() {
        entries.push(last);
    }

    Ok(entries)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_splitb() {
        let line = "chr11\t13934505\t13958243\tENST00000674667\t1000\t-\t13934505\t13958243\t0,0,200\t9\t224,217,228,198,149,142,115,157,49,\t0,1305,2811,5576,10085,14837,18016,19498,23689,";
        let list = splitb(line, "\t").unwrap();
        assert_eq!(list.len(), 12);
        assert_eq!(
            list,
            [
                "chr11",
                "13934505",
                "13958243",
                "ENST00000674667",
                "1000",
                "-",
                "13934505",
                "13958243",
                "0,0,200",
                "9",
                "224,217,228,198,149,142,115,157,49,",
                "0,1305,2811,5576,10085,14837,18016,19498,23689,"
            ]
        );
    }

    #[test]
    fn new_record() {
        let line =
            "chr15\t81000922\t81005788\tENST00000267984\t0\t+\t81002271\t81003360\t0\t1\t4866,\t0,";
        let record = BedRecord::new(line).unwrap();

        assert_eq!(record.chrom(), "chr15");
        assert_eq!(record.tx_start(), 81000922);
        assert_eq!(record.tx_end(), 81005788);
        assert_eq!(record.name(), "ENST00000267984");
        assert_eq!(record.strand(), "+");
        assert_eq!(record.cds_start(), 81002271);
        assert_eq!(record.cds_end(), 81003360);
        assert_eq!(record.exon_count(), 1);
        assert_eq!(record.exon_start(), &vec![81000922]);
        assert_eq!(record.exon_end(), &vec![81005788]);
    }

    #[test]
    fn get_exon_frames() {
        let line = "chr11\t13934505\t13958243\tENST00000674667\t1000\t-\t13934505\t13958243\t0,0,200\t9\t224,217,228,198,149,142,115,157,49,\t0,1305,2811,5576,10085,14837,18016,19498,23689,";
        let record = BedRecord::new(line).unwrap();

        assert_eq!(record.get_exon_frames(), vec![1, 0, 0, 0, 1, 0, 2, 1, 0]);
    }

    #[test]
    fn invalid_record() {
        let line =
            "chr15\t81000922\t81005788\tENST00000267984\t0\t+\t81002271\t81003360\t0\t1\t4866,";
        let record = BedRecord::new(line);

        assert_eq!(record, Err(ParseError::Invalid));
    }

    #[test]
    fn empty_record() {
        let line = "";
        let record = BedRecord::new(line);

        assert_eq!(record, Err(ParseError::Empty));
    }
}
