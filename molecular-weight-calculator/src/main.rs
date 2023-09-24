/*
The Goal: Produce a CLI to calculate molecular mass from a formula.

Options:
    -e, --exact (default) (NOT REQUIRED!)
    -a, --average
    -p, --protein, --peptide => Calculate as a polypeptide chain. (IUPAC Letters)
    -d, --dna => Calculate as DNA. (exact mode requires no ambiguous IUPAC bases, average mode allows them and averages their mass).
    -r, --rna => Calculate as RNA. (same as above).
    -F, --file => Input file (read based on a .fasta, or .csv file (column must be named "formula")).
    -h, --help => Show this
    -f, --format => Show help for file formatting.
*/

use clap::{arg, error::ErrorKind, ArgGroup, Args, Command};
use regex::Regex;
use std::path::PathBuf;

#[derive(Debug, Args)]
#[command(author, version, about, long_about = None)]
struct Arguments {
    #[arg(
        help = "The molecular or IUPAC formula to calculate the mass from.",
        required = false
    )]
    formula: String,
}

fn main() {
    let cli = Command::new("calcm")
        .author("Trevor Fox, voxelmc2@student.ubc.ca")
        .version("1.0.0")
        .about("Calculates the molecular mass of a given formula! NOTE: CURRENTLY DEFAULTED TO CALCULATE AVERAGE VALUES. EXACT VALUES WILL BE AVAILABLE IN THE FUTURE.")
        .arg(arg!(-a --average "Requests the average mass. (Enables ambiguous bases when calculating nucleotides)"))
        .arg(arg!(-d --dna "Notify the calculator to parse the formula as IUPAC bases. Accepts [A, T, G, C]"))
        .arg(arg!(-r --rna "Notify the calculator to parse the formula as IUPAC bases. Accepts [A, U, G, C]."))
        .arg(arg!(-p --protein "Notify the calculator to parse the formula as single-letter IUPAC amino acids."))
        .arg(arg!(-f --format "Shows information about the required file format for file-based input."))
        .arg(arg!(-F --file <FILE> "Specify a file for the calculator to read the formula from."))
        .arg(arg!(--debug "Launch in verbose debugging mode."))
        .group(
            ArgGroup::new("molecule").args(["dna", "rna", "protein"])
        );

    let err = cli.clone().error(
        ErrorKind::MissingRequiredArgument,
        "Please specify a molecule! (or use the \"--file <FILE>\" flag).",
    );

    let cli = Arguments::augment_args(cli);

    let matches = cli.get_matches();
    let molecule = match matches.get_raw("formula") {
        Some(mut raw) => raw.next().to_owned(),
        None => {
            let _ = err.print();
            std::process::exit(1);
        }
    }
    .unwrap();

    println!("Calculating mass for: {:?}", molecule);

    let formula_vec: Vec<String> = parse_formula(
        molecule
            .to_str()
            .expect("Could not cast molecule into String")
            .to_owned(),
    );

    let output = mass_from_formula(formula_vec, matches.get_flag("debug"));
    println!("{:?}", output);
}

fn parse_formula<'a>(formula: String) -> Vec<String> {
    let reg = Regex::new(r"[A-Za-z][a-z]?\d*|(<!\([^)])\(.*\)\d+(![^(]*\))")
        .expect("RegEx parsing error.");
    let binding = reg.to_owned();
    let out = binding.find_iter(formula.as_str());

    let mut out_vec: Vec<String> = vec![];
    out.for_each(|val| {
        let in_vec: Vec<String> = vec![val.as_str().to_owned()];
        out_vec.append(in_vec.to_vec().as_mut());
    });

    out_vec
}

fn mass_from_formula<'ass>(formula: Vec<String>, debug: bool) -> f32 {
    let mut aggregate_mass: f32 = 0.0;

    for atom in formula {
        let reg = Regex::new(r"(\D+)|(\d+)").expect("RegEx parsing error.");
        let matches: Vec<String> = reg
            .find_iter(&atom)
            .map(|val| val.as_str().to_owned())
            .collect();

        let elements_csv_path: PathBuf = PathBuf::from("./data/elements2.csv");
        let elements_csv_stream = csv::Reader::from_path(elements_csv_path);
        let mut deserialize_binding = elements_csv_stream.unwrap();
        let mut elements_csv_deserialized = deserialize_binding.deserialize::<ElementRow>();

        if debug {
            println!("{:?}", matches[0].to_uppercase());
            println!("{:?}", matches);
        }

        let element: Option<Result<ElementRow, csv::Error>> = elements_csv_deserialized
            .find(|val| val.as_ref().unwrap().element == matches[0].to_uppercase());

        let element_weight = element
            .expect("1")
            .expect("2")
            .weight
            .parse::<f32>()
            .expect("Could not parse element mass from dataset as a float32.");

        if matches[1].is_empty() {
            println!("here");
            aggregate_mass += element_weight;
        } else {
            let element_count = matches[1]
                .parse::<f32>()
                .expect("Could not parse element count into a float32.");
            aggregate_mass += element_weight * element_count;
        }
    }

    aggregate_mass
}

#[derive(Debug, serde::Deserialize, Clone)]
struct ElementRow {
    element: String,
    weight: String,
    // uncertainty: String,
    // charge: i32,
}

// let element_abbr = &matches[0];

// println!("{:?}", matches);
// // println!("element abbr {:?}", element_abbr);
// println!(
//     "MANUAL H {:?}",
//     elements_csv_deserialized.find(|val| val.as_ref().expect("Couldn't.").element == "K")
// );
// // println!(
// //     "{:?}",
// //     elements_csv_deserialized
// //         .find(|val| val.as_ref().expect("Couldn't.").element == element_abbr.as_ref())
// // );

// let element: Option<Result<ElementRow, csv::Error>> =
//     elements_csv_deserialized.find(|val| val.as_ref().unwrap().element == matches[0]);

// let a = match element {
//     Some(aa) => {
//         // println!("{:?}", aa);
//         aa
//     }
//     None => panic!("asdddd"),
// };

// let ab = match a {
//     Ok(dd) => {
//         // println!("{:?}", &dd);
//         dd
//     }
//     Err(e) => panic!("asd"),
// };

// let mut found: Vec<ElementRow> = vec![];
// for element in elements_csv_tosort {
//     if element.element == matches[0] {
//         found.insert(found.len(), element);
//         break;
//     } else {
//         continue;
//     }
// }
