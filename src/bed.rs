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