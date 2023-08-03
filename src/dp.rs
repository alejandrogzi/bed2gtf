use std::path::PathBuf;
use std::fs;
use std::fs::Permissions;
use std::os::unix::fs::PermissionsExt;

use crate::sh::Shell;


const BED_TO_GENEPRED: &str = "bedToGenePred";
const GENEPRED_TO_GTF: &str = "genePredToGtf";
const HTTP: &str = "http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/";


pub struct Dependencies;


impl Dependencies {

    // creates a directory called "dependencies" one level above cwd
    fn make_dir() -> Result<PathBuf, std::io::Error> {
        let mut path = std::env::current_dir().unwrap();
        path.push("dependencies");

        if path.exists() {
            println!("{} already exists", path.display());
            return Ok(path);
        } else {
            std::fs::create_dir_all(&path)?;
            println!("Dependencies will be download at: {}", path.display());
            return Ok(path);
        }
    }

    // downloads a given binary in the "dependencies" directory
    fn get_binary(path: &PathBuf, binary: &str) {
        let url: String = [HTTP, binary].join("");
        let cmd: String = format!("wget -P {} {}", path.display(), url);
        let _ = Shell::run(&cmd);
    }


    // checks if binaries are already installed
    fn check_existence(path: &PathBuf, binary: &str) -> bool {
        let mut x = path.clone();
        x.push(binary);
        return x.exists();
    }


    // makes a binary executable
    fn make_executable(file: &PathBuf) {
        fs::set_permissions(file, Permissions::from_mode(0o755)).unwrap();
        println!("changed {:?} mode", file.file_name().unwrap());
    }


    // runs the whole script
    pub fn get() -> Vec<PathBuf> {
        let path: PathBuf = Dependencies::make_dir().unwrap();
        let binaries: Vec<&str> = vec![BED_TO_GENEPRED, GENEPRED_TO_GTF];
        let mut output = Vec::new();

        for binary in binaries {
            if Dependencies::check_existence(&path, binary) {
                println!("{} already exists", binary);

                let file = path.join(&binary);            
                Dependencies::make_executable(&file);
                output.push(file);
            } else {
                Dependencies::get_binary(&path, binary);

                let file = path.join(&binary);            
                Dependencies::make_executable(&file);
                output.push(file);
            }
        }

        println!("\nDependencies installed successfully!");
        return output;
    }

}