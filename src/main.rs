use clap::Parser;
use rayon::prelude::*;
use std::path::Path;
use linfa::prelude::*;
use linfa_elasticnet::{ElasticNet, Result};
use linfa_elasticnet::{ElasticNetError, ElasticNetParams};
use log::LevelFilter;
use log::*;
use ndarray::{array, s, Array, Array2};
use prita::cmdline::*;
use prita::load::*;
use prita::sketch::*;
use prita::types::*;
use std::collections::HashMap;
use std::fs::read;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};

pub fn sketch(args: Sketch) {
    if !args.queries.is_none() {
        let queries = args.queries.as_ref().unwrap();
        let iter_vec : Vec<usize> = (0..queries.len()).into_iter().collect();
        iter_vec.into_par_iter().for_each(|i| {
            let query_file = &queries[i];
            let pref = Path::new(&args.query_sketch_output_prefix);
            let query_file_name = Path::new(query_file).file_name().unwrap();
            let file_path = pref.join(query_file_name);
            let file_path_str = format!("{}.prs", file_path.to_str().unwrap());
            let mut query_sk_file = BufWriter::new(File::create(&file_path_str).expect(&format!("{} not valid", file_path_str)));

            let query_sketch = sketch_query(&args, query_file);
            info!("Sketching query done. Serializing...");
            let enc = SequencesSketchEncode::new(query_sketch);
            bincode::serialize_into(&mut query_sk_file, &enc).unwrap();
        });
    }
    if !args.references.is_none() || !args.reference_list.is_none() {
        if !args.references.is_none() && !args.reference_list.is_none() {
            panic!("Only one of --rl or -r can be specified.");
        }
        let references_sketch = sketch_references(&args);
        info!("Sketching reference done. Serializing...");

        let mut ref_sk_file =
            BufWriter::new(File::create(args.reference_sketch_output.clone()).unwrap());
        bincode::serialize_into(&mut ref_sk_file, &references_sketch).unwrap();
    }
    info!("Finished.");
}

//pub fn strain(args: Strain) {
//    let sequence_bytes = read(args.sequences_sketch).expect("Sequence sketch not a valid file");
//    let sequence_sketch: SequencesSketch = bincode::deserialize(&sequence_bytes).unwrap();
//
//    let genome_bytes = read(args.genomes_sketch).expect("Genome sketch not a valid file");
//    let genome_sketch: GenomesSketch = bincode::deserialize(&genome_bytes).unwrap();
//    dbg!(&genome_sketch.file_names);
//
//    let mut strain_kmer_map: HashMap<u64, usize> = HashMap::default();
//    let mut shared_kmer_map: HashMap<u64, usize> = HashMap::default();
//    let mut mean_cov = 0;
//    let mut cov_count = 0;
//
//    for kmer_count in genome_sketch.genome_kmer_counts.iter() {
//        for kmer in kmer_count.keys() {
//            let c = shared_kmer_map.entry(*kmer).or_insert(0);
//            *c += 1;
//            if !strain_kmer_map.contains_key(kmer) {
//                strain_kmer_map.insert(*kmer, strain_kmer_map.len());
//                if sequence_sketch.kmer_counts.contains_key(kmer) {
//                    mean_cov += sequence_sketch.kmer_counts[kmer];
//                    cov_count += 1;
//                }
//            }
//        }
//    }
//
//    let mean_cov = mean_cov as f64 / cov_count as f64;
//
//    let mut target_vec = vec![0.; strain_kmer_map.len()];
//    let mut strain_kmer_mat =
//        Array2::<f64>::zeros((strain_kmer_map.len(), genome_sketch.file_names.len()));
//
//    dbg!(shared_kmer_map.len());
//
//    let mut shared_count = 0;
//    let mut counts = vec![0; genome_sketch.file_names.len()];
//    let mut mean_counts = vec![0.; genome_sketch.file_names.len()];
//    for (i, kmer_count) in genome_sketch.genome_kmer_counts.iter().enumerate() {
//        for kmer in kmer_count.keys() {
//            if shared_kmer_map[kmer] != genome_sketch.genome_kmer_counts.len() {
//                strain_kmer_mat[[strain_kmer_map[kmer], i]] = 1.;
//                if sequence_sketch.kmer_counts.contains_key(kmer) {
//                    target_vec[strain_kmer_map[kmer]] =
//                        sequence_sketch.kmer_counts[kmer] as f64 / mean_cov as f64;
//                    counts[i] += 1;
//                    mean_counts[i] += sequence_sketch.kmer_counts[kmer] as f64;
//                } else {
//                    target_vec[strain_kmer_map[kmer]] = 0.;
//                }
//            } else {
//                shared_count += 1;
//            }
//        }
//    }
//    dbg!(shared_count);
//    dbg!(mean_counts);
//    dbg!(counts);
//
//    let ds = Dataset::new(strain_kmer_mat, Array::from_vec(target_vec));
//    let model1 = ElasticNetParams::new()
//        .penalty(1. / (genome_sketch.file_names.len() as f64))
//        .l1_ratio(1.0)
//        .with_intercept(false);
//
//    let result = model1.fit(&ds).unwrap();
//    dbg!(result.hyperplane());
//}

fn main() {
    simple_logging::log_to_stderr(LevelFilter::Info);
    rayon::ThreadPoolBuilder::new()
        .num_threads(10)
        .build_global()
        .unwrap();

    let cli = Cli::parse();
    match cli.mode {
        Mode::Sketch(Sketch) => sketch(Sketch),
        Mode::Load(Load) => load(Load),
        Mode::Strain(Strain) => panic!(),
    }
}
