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
    exon_end: Vec<i32>
}


impl BedRecord {
    pub fn new(fields: Vec<&str>) -> BedRecord {
        let chrom = fields[0].to_string();
        let tx_start = fields[1].parse::<i32>().unwrap();
        let tx_end = fields[2].parse::<i32>().unwrap();
        let name = fields[3].to_string();
        let strand = fields[5].to_string();
        let cds_start = fields[6].parse::<i32>().unwrap();
        let cds_end = fields[7].parse::<i32>().unwrap();
        let exon_count: i16 = fields[9].parse::<i16>().unwrap();
        let mut exon_start: Vec<i32> = fields[11].split(',').filter(|&x| !x.is_empty()).map(|x| x.parse::<i32>().unwrap()).collect();
        let mut exon_end: Vec<i32> = fields[10].split(',').filter(|&x| !x.is_empty()).map(|x| x.parse::<i32>().unwrap()).collect();

        for x in 0..exon_count as usize {
            exon_start[x] += tx_start;
            exon_end[x] += exon_start[x];
        }
        
        BedRecord {
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
        }
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


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_record() {
        let fields = vec![
            "chr1",
            "11869",
            "14409",
            "DDX11L1",
            "0",
            "+",
            "11869",
            "11869",
            "0",
            "3",
            "354,120,247",
            "11869,12612,13220",
        ];

        let record = BedRecord::new(fields);

        assert_eq!(record.chrom(), "chr1");
        assert_eq!(record.tx_start(), 11869);
        assert_eq!(record.tx_end(), 14409);
        assert_eq!(record.name(), "DDX11L1");
        assert_eq!(record.strand(), "+");
        assert_eq!(record.cds_start(), 11869);
        assert_eq!(record.cds_end(), 11869);
        assert_eq!(record.exon_count(), 3);
        assert_eq!(record.exon_start(), &vec![23738, 24481, 25089]);
        assert_eq!(record.exon_end(), &vec![24092, 24601, 25336]);
    }

    #[test]
    fn get_exon_frames() {
        let fields = vec![
            "chr1",
            "164680187",
            "164841823",
            "ABC",
            "0",
            "+",
            "164685945",
            "164841196",
            "0",
            "14",
            "5845,177,1568,95,235,75,154,237,98,229,202,189,91,114",
            "0,12033,14081,19217,20334,24256,24764,41887,42237,47781,49726,63499,160968,161522",
        ];

        let record = BedRecord::new(fields);

        assert_eq!(record.get_exon_frames(), vec![0, 0, 0, 2, 1, 2, 2, 0, 0, 2, 0, 1, 1, -1]);
    }
}