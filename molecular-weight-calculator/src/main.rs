use clap::Command;
use clap::Parser;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Name of the person to greet
    #[arg(short, long)]
    name: String,

    /// Number of times to greet
    #[arg(short, long, default_value_t = 1)]
    count: u8,
}

fn main() {
    let command = Command::new("test");
}

/*
The Goal: Produce a CLI to calculate molecular mass from a formula.

Options:
    -e, --exact (default)
    -a, --average
    -p, --protein, --peptide => Calculate as a polypeptide chain. (IUPAC Letters)
    -d, --dna => Calculate as DNA. (exact mode requires no ambiguous IUPAC bases, average mode allows them and averages their mass).
    -r, --rna => Calculate as RNA. (same as above).
    -F, --file => Input file (read based on a .fasta, or .csv file (column must be named "formula")).
    -h, --help => Show this
    -f, --format => Show help for file formatting.
*/
