use thiserror::Error;

#[derive(Error, Debug, PartialEq)]
pub enum ParseError {
    #[error("Empty line")]
    Empty,
    #[error("Invalid GTF line")]
    Invalid,
}

impl From<std::io::Error> for ParseError {
    fn from(error: std::io::Error) -> Self {
        ParseError::Invalid;
        ParseError::Empty
    }
}
