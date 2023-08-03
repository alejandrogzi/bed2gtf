use std::process::Command;


pub struct Shell;

impl Shell {
    pub fn run(cmd: &str) -> Result<String, std::io::Error> {
        let output = Command::new("sh")
            .arg("-c")
            .arg(cmd)
            .output()?;

        if output.status.success() {
            Ok(String::from_utf8_lossy(&output.stdout).into_owned())
        } else {
            Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!("Command exited with non-zero status code:\n{}", cmd),
            ))
        }
    }
}
