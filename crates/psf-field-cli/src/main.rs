//! `psf-field` CLI (C1B.7). Subcommands are stubs until later PRs.

const SUBCOMMANDS: &[&str] = &[
    "stage1",
    "stage2",
    "eval",
    "score",
    "check-jacobian",
    "report",
];

fn main() {
    let mut args = std::env::args().skip(1);
    match args.next().as_deref() {
        None | Some("-h") | Some("--help") => {
            println!("Usage: psf-field <subcommand>");
            println!("Subcommands: {}", SUBCOMMANDS.join(", "));
        }
        Some(cmd) if SUBCOMMANDS.contains(&cmd) => {
            eprintln!("psf-field {cmd}: not implemented");
            std::process::exit(2);
        }
        Some(other) => {
            eprintln!("unknown subcommand: {other}");
            eprintln!("Subcommands: {}", SUBCOMMANDS.join(", "));
            std::process::exit(2);
        }
    }
}
