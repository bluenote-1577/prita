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
    #[arg(short, num_args = 0..)]
    pub queries: Option<Vec<String>>,
    #[arg(short,long, value_parser, num_args = 0..)]
    pub references : Option<Vec<String>>,
    #[arg(long="rl")]
    pub reference_list: Option<String>,
    #[arg(short,long, default_value_t = 0.99)]
    pub derep: f64,
    #[arg(long="ref-output", default_value = "ref_sketch.prs")]
    pub reference_sketch_output: String,
    #[arg(long="query-prefix", default_value = "")]
    pub query_sketch_output_prefix: String,
    #[arg(short, default_value_t = 31)]
    pub k: usize,
    #[arg(short, default_value_t = 100)]
    pub c: usize,
}

#[derive(Args)]
pub struct Load {
    #[arg(short, num_args = 0..)]
    pub query_sketch : Vec<String>,
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
