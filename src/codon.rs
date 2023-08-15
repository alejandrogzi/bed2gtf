
#[derive(Debug, Clone)]
pub struct Codon {
    pub start: i32,
    pub end: i32,
    pub index: i32,
    pub start2: i32,
    pub end2: i32,
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