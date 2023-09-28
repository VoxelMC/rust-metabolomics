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
use colored::Colorize;
use regex::Regex;
use std::{env, path::PathBuf};
use std::num::ParseFloatError;

#[derive(Debug, serde::Deserialize, Clone)]
struct ElementRow {
    element: String,
    weight: String,
    // uncertainty: String,
    // charge: i32,
}

#[derive(Debug, serde::Deserialize, Clone)]
struct AbbreviationRow {
    abbreviation: String,
    formula: String,
    // charge: i32,
    // name: String,
}

#[derive(Debug, serde::Deserialize, Clone)]
struct AminoAcidRow {
    letter: String,
    formula: String,
    // ,charge,letter,name
    // charge: i32,
    // name: String,
}

#[derive(Debug, serde::Deserialize, Clone)]
struct NucleicAcidRow {
    letter: String,
    formula: String,
    // ,charge,letter,name
    // charge: i32,
}

#[derive(Debug, Args)]
#[command(author = "Trevor Fox, voxelmc2@student.ubc.ca", version = "1.0.0", about, long_about = None)]
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
        .about(
            format!(
                "{} | {}\n{}\n{}", 
                "cmo".bright_cyan().bold(), 
                "v1.0.0".bright_green(), 
                "Calculates the molecular mass of a given formula!", 
                "Developed by Trevor Fox.".bold()
            )
        )
        .arg(arg!(-a --average "Requests the average mass (enables ambiguous bases when calculating nucleotides)."))
        .arg(arg!(-d --dna "Notify the calculator to parse the formula as IUPAC bases. Accepts [A, T, G, C]"))
        .arg(arg!(-r --rna "Notify the calculator to parse the formula as IUPAC bases. Accepts [A, U, G, C]."))
        .arg(arg!(-p --protein "Notify the calculator to parse the formula as single-letter IUPAC amino acids."))
        .arg(arg!(-f --format "Shows information about the required file format for file-based input."))
        .arg(arg!(-F --file <FILE> "Specify a file for the calculator to read the formula from."))
        .arg(arg!(--debug "Launch in verbose debugging mode."))
        .arg(arg!(-s --silent "Print only the numeric mass."))
        .group(
            ArgGroup::new("molecule").args(["dna", "rna", "protein"])
        );

    let err = cli.clone().error(
        ErrorKind::MissingRequiredArgument,
        "Please specify a molecule! (or use the \"--file <FILE>\" flag).",
    );

    let cli = Arguments::augment_args(cli);

    let h_mass: f32 = 1.00794;
    let matches = cli.get_matches();
    let molecule_string = match matches.get_raw("formula") {
        Some(mut raw) => raw
            .next()
            .to_owned()
            .expect("Could not cast molecule into OsString")
            .to_str()
            .expect("Could not cast molecule into String")
            .to_owned(),
        None => {
            let _ = err.print();
            std::process::exit(1);
        }
    };

    let is_average = matches.get_flag("average");
    let is_debug = matches.get_flag("debug");

    if !matches.get_flag("silent") {
        println!(
            "Calculating {} mass for: {}",
            if matches.get_flag("average") {
                "average"
            } else {
                "exact"
            },
            molecule_string.bright_yellow().bold()
        );
    }

    if matches.get_flag("protein") {
        let formula_vec = parse_protein_formula(molecule_string);
        let mut number_of_hydrogens = -2.0;
        let _ = formula_vec.iter().for_each(|val| {
            if val.contains("H") {
                let reg = Regex::new(r"(\d+)").expect("RegEx parsing error.");
                let matched = reg
                    .find(&val).expect("").as_str().to_owned();
                let parsed: Result<f32, ParseFloatError> = matched.parse::<f32>();
                number_of_hydrogens += parsed.expect("Could not parse # of hydrogens as a float.");
            }
        });
            // .collect::<Vec<&String>>().len() as i32;
        let mass_to_subtract: f32 = number_of_hydrogens * h_mass;
        let output = mass_from_formula(formula_vec, is_debug, is_average) - mass_to_subtract;
        println!("{:?}", output);
        return;
    }

    if matches.get_flag("dna") {
        let formula_vec = parse_nucleic_formula(molecule_string, false, is_average);
        let output = mass_from_formula(formula_vec, is_debug, is_average);
        println!("{:?}", output);
        return;
    }

    if matches.get_flag("rna") {
        let formula_vec = parse_nucleic_formula(molecule_string, true, is_average);
        let output = mass_from_formula(formula_vec, is_debug, is_average);
        println!("{:?}", output);
        return;
    }

    let formula_vec: Vec<String> = parse_molecular_formula(molecule_string);

    let output = mass_from_formula(
        formula_vec,
        matches.get_flag("debug"),
        matches.get_flag("average"),
    );
    println!("{:?}", output);
}

/// CH3 -> ["C", "H3"]
fn parse_molecular_formula<'a>(formula: String) -> Vec<String> {
    let reg = Regex::new(r"[A-Za-z][a-z]{0,2}\d*|(<!\([^)])\(.*\)\d+(![^(]*\))")
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

fn mass_from_formula<'ass>(parsed_formula: Vec<String>, is_debug: bool, is_average: bool) -> f32 {
    let mut aggregate_mass: f32 = 0.0;

    for atom in parsed_formula {
        let reg = Regex::new(r"(\D+)|(\d+)").expect("RegEx parsing error.");
        let matches: Vec<String> = reg
            .find_iter(&atom)
            .map(|val| val.as_str().to_owned())
            .collect();

        let current_exe_res = env::current_exe();
        let mut current_exe_path: PathBuf =
            current_exe_res.expect("Could not read executable path.");
        current_exe_path.pop();

        let elements_csv_path: PathBuf = if is_average {
            current_exe_path.join("../../data/elements.csv")
        } else {
            current_exe_path.join("../../data/absolute.csv")
        }
        .canonicalize()
        .expect("Canonicalization of executable path failed.");

        // Try mapping to clones to remake each time, instead of new reader.
        let elements_csv_stream = csv::Reader::from_path(elements_csv_path);
        let mut elements_deserialize_binding = elements_csv_stream.unwrap();
        let mut elements_csv_deserialized =
            elements_deserialize_binding.deserialize::<ElementRow>();

        let abbr_csv_path: PathBuf = current_exe_path
            .join("../../data/abbreviations.csv")
            .canonicalize()
            .expect("Canonicalization of executable path failed.");
        let abbr_csv_stream = csv::Reader::from_path(abbr_csv_path);
        let mut abbr_deserialize_binding = abbr_csv_stream.unwrap();
        let mut abbr_csv_deserialized = abbr_deserialize_binding.deserialize::<AbbreviationRow>();

        let _expanded_abbr = match abbr_csv_deserialized.find(|row| {
            row.as_ref()
                .expect("Could not get AbbreviationRow")
                .abbreviation
                == matches[0]
        }) {
            Some(res) => match res {
                Ok(found) => {
                    let parsed = parse_molecular_formula(found.formula);
                    if is_debug {
                        println!("Expanded Abbreviation: {:?}", parsed)
                    };
                    let weight = mass_from_formula(parsed, is_debug, is_average);

                    if matches.len().eq(&1) {
                        aggregate_mass += weight;
                    } else if matches.len().eq(&2) {
                        let element_count = matches[1]
                            .parse::<f32>()
                            .expect("Could not parse element count into a float32.");
                        aggregate_mass += weight * element_count;
                    }
                    continue;
                }
                Err(_e) => (),
            },
            None => (),
        };

        if is_debug {
            println!("Element: {:?}", matches[0]);
            println!("Matches: {:?}", matches);
            println!("Matches Len: {:?}", matches.len());
            println!("Aggregated Mass: {:?}", aggregate_mass);
        }

        let element: Option<Result<ElementRow, csv::Error>> =
            elements_csv_deserialized.find(|val| val.as_ref().unwrap().element == matches[0]);

        let element_weight = element
            .expect("1")
            .expect("2")
            .weight
            .parse::<f32>()
            .expect("Could not parse element mass from dataset as a float32.");
        if matches.len().eq(&1) {
            aggregate_mass += element_weight;
        } else if matches.len().eq(&2) {
            let element_count = matches[1]
                .parse::<f32>()
                .expect("Could not parse element count into a float32.");
            aggregate_mass += element_weight * element_count;
        }
    }
    aggregate_mass
}

/// Gives deprotonated and protonated (M+2).
fn parse_protein_formula<'a>(formula: String) -> Vec<String> {
    let reg = Regex::new(r"[Aa]|[C-Yc-y]").expect("RegEx parsing error.");
    let binding = reg.to_owned();
    let out = binding.find_iter(formula.as_str());

    let current_exe_res = env::current_exe();
    let mut current_exe_path: PathBuf = current_exe_res.expect("Could not read executable path.");
    current_exe_path.pop();

    let mut out_vec: Vec<String> = vec![];
    out.for_each(|val| {
        let aa_csv_path: PathBuf = current_exe_path
            .join("../../data/amino.csv")
            .canonicalize()
            .expect("Canonicalization of executable path failed.");
        let aa_csv_stream = csv::Reader::from_path(aa_csv_path);
        let aa_deserialize_binding = aa_csv_stream.unwrap();
        let mut aa_csv_deserialized = aa_deserialize_binding.into_deserialize::<AminoAcidRow>();

        let aa_formula = aa_csv_deserialized
            .find(|aa_row| aa_row.as_ref().unwrap().letter == val.as_str().to_owned());
        let in_vec: Vec<String> = parse_molecular_formula(aa_formula.unwrap().unwrap().formula);

        out_vec.append(in_vec.to_vec().as_mut());
    });
    out_vec
}

fn parse_nucleic_formula<'a>(formula: String, is_rna: bool, is_average: bool) -> Vec<String> {
    let reg = if is_rna {
        Regex::new(r"[AUGCaugc]")
    } else {
        Regex::new(r"[ATGCatgc]")
    }
    .expect("RegEx parsing error.");

    let binding = reg.to_owned();
    let out = binding.find_iter(formula.as_str());

    let current_exe_res = env::current_exe();
    let mut current_exe_path: PathBuf = current_exe_res.expect("Could not read executable path.");
    current_exe_path.pop();

    let mut out_vec: Vec<String> = vec![];
    out.for_each(|val| {
        let na_csv_path: PathBuf = if is_average {
            current_exe_path.join("../../data/nucleic.csv")
        } else {
            current_exe_path.join("../../data/deoxyribonucleic.csv")
        }
        .canonicalize()
        .expect("Canonicalization of executable path failed.");

        let na_csv_stream = csv::Reader::from_path(na_csv_path);
        let na_deserialize_binding = na_csv_stream.unwrap();
        let mut na_csv_deserialized = na_deserialize_binding.into_deserialize::<NucleicAcidRow>();

        let na_formula = na_csv_deserialized.find(|na_row| {
            na_row.as_ref().unwrap().letter == val.as_str().to_owned().to_uppercase()
        });
        let in_vec: Vec<String> = parse_molecular_formula(na_formula.unwrap().unwrap().formula);

        out_vec.append(in_vec.to_vec().as_mut());
    });
    out_vec
}
