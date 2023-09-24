use polars::prelude::*;

fn main() {
    use std::fs;
    use std::path::{Display, Path};
    println!("Hello, world!");

    use std::env;

    let current_exe = env::current_exe();
    let current_dir: Display = match current_exe {
        Ok(ref exe_path) => exe_path.display(),
        Err(e) => panic!("{e}"),
    };

    println!("{}", current_dir.to_string());

    let paths = fs::read_dir("./").unwrap();

    for path in paths {
        println!("Name: {}", path.unwrap().path().display())
    }

    // let file: File = File::open("./Wine Dataset 2023.csv").expect("Could not open this file.");

    let csv_reader = CsvReader::from_path(Path::new("./Wine Dataset 2023.csv"));
    let df_reader = match csv_reader {
        Ok(data) => data.truncate_ragged_lines(true).finish(),
        Err(e) => panic!("{e}"),
    };

    let df = match df_reader {
        Ok(data) => data,
        Err(e) => panic!("{e}"),
    };

    let exploded = df.explode(["TH_T_50"]).unwrap();

    // println!("{df}");
    println!("{exploded}");
}
