fn main() {
    println!("Hello, world!");
}

use polars_core::prelude::*;

let df = CsvReader::from_path("../Wine Dataset 2023.csv").unwrap().finish().unwrap();

