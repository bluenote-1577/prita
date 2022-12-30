use clap::{Args, Parser, Subcommand};

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub mode: Mode,
}

#[derive(Subcommand)]
pub enum Mode {
    /// Adds files to myapp
    Sketch(Sketch),
    Load(Load),
    Strain(Strain)
}

#[derive(Args)]
pub struct Sketch {
    #[arg(short)]
    pub query: Option<String>,
    #[arg(short,long, value_parser, num_args = 0..)]
    pub references : Option<Vec<String>>,
    #[arg(long="rl")]
    pub reference_list: Option<String>,
    #[arg(short,long, default_value_t = 0.99)]
    pub derep: f64,
    #[arg(long="ro", default_value = "reference_sketch.prs")]
    pub reference_sketch_output: String,
    #[arg(long="qo", default_value = "query_sketch.pqs")]
    pub query_sketch_output: String,
    #[arg(short, default_value_t = 31)]
    pub k: usize,
    #[arg(short, default_value_t = 1000)]
    pub c: usize,
}

#[derive(Args)]
pub struct Load {
    #[arg(short)]
    pub query_sketch : String,
    #[arg(short)]
    pub references_sketch:  String,
    #[arg(short, default_value_t = 0.95)]
    pub min_ani: f64,
    #[arg(long, default_value_t = 0.95)]
    pub cluster: f64,

}

#[derive(Args)]
pub struct Strain {
    #[arg(short)]
    pub sequences_sketch : String,
    #[arg(short)]
    pub genomes_sketch:  String,
}
