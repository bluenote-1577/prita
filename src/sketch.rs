use crate::cmdline::*;
use crate::seeding::*;
use crate::types::*;
use flate2::read::GzDecoder;
use log::*;
use needletail::parse_fastx_file;
use nohash_hasher::IntMap;
use nohash_hasher::IntSet;
use rayon::prelude::*;
use seq_io::fasta::{Reader as ReaderA, Record as ReccordA};
use seq_io::fastq::{Reader as ReaderQ, Record as RecordQ};
use seq_io::parallel::read_parallel;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::Read;
use std::io::Write;
use std::io::{BufRead, BufReader};
use std::sync::Mutex;

pub fn sketch_query_needle(args: &Sketch) -> SequencesSketch {
    let read_file = &args.query.as_ref().unwrap();
    let mut kmer_map = HashMap::default();
    let ref_file = &read_file;
    let reader = parse_fastx_file(&ref_file);
    let mut vec = vec![];
    if !reader.is_ok() {
        warn!("{} is not a valid fasta/fastq file; skipping.", ref_file);
    } else {
        let mut reader = reader.unwrap();
        while let Some(record) = reader.next() {
            if record.is_ok() {
                let record = record.expect(&format!("Invalid record for file {}", ref_file));
                let seq = record.seq();
                unsafe {
                    extract_markers_avx2(&seq, &mut vec, args.k, args.c);
                }
            } else {
                warn!("File {} is not a valid fasta/fastq file", ref_file);
            }
        }
        for km in vec {
            let c = kmer_map.entry(km).or_insert(0);
            *c += 1;
        }
    }

    return SequencesSketch {
        kmer_counts: kmer_map,
        file_name: read_file.to_string(),
        c: args.c,
        k: args.k,
    };
}

pub fn sketch_query(args: &Sketch) -> SequencesSketch {
    let read_file = &args.query.as_ref().unwrap();
    let mut read_sketch = SequencesSketch::new(read_file.to_string(), args.c, args.k);
    if read_file.contains(".fq") || read_file.contains(".fastq") {
        let reader;
        if read_file.contains(".gz") || read_file.contains(".gzip") {
            let file = File::open(read_file).unwrap();
            let gz_decoder: Box<dyn Read + Send> = Box::new(BufReader::new(GzDecoder::new(file)));
            reader = ReaderQ::new(gz_decoder);
        } else {
            let file = File::open(read_file).unwrap();
            let decoder: Box<dyn Read + Send> = Box::new(BufReader::new(file));
            reader = ReaderQ::new(decoder);
        }
        read_parallel(
            reader,
            20,
            200,
            |record_set| {
                // this function does the heavy work
                let mut vec = vec![];
                unsafe {
                    for record in record_set.into_iter() {
                        extract_markers_avx2(record.seq(), &mut vec, args.c, args.k);
                    }
                }

                return vec;
            },
            |record_sets| {
                while let Some(result) = record_sets.next() {
                    let (_rec_set, vec) = result.unwrap();
                    for km in vec {
                        let c = read_sketch.kmer_counts.entry(km).or_insert(0);
                        *c += 1;
                    }
                }
            },
        );
    } else {
        let reader;
        if read_file.contains(".gz") || read_file.contains(".gzip") {
            let file = File::open(read_file).unwrap();
            let gz_decoder: Box<dyn Read + Send> = Box::new(BufReader::new(GzDecoder::new(file)));
            reader = ReaderA::new(gz_decoder);
        } else {
            let file = File::open(read_file).unwrap();
            let decoder: Box<dyn Read + Send> = Box::new(BufReader::new(file));
            reader = ReaderA::new(decoder);
        }
        read_parallel(
            reader,
            20,
            200,
            |record_set| {
                // this function does the heavy work
                let mut vec = vec![];
                unsafe {
                    for record in record_set.into_iter() {
                        extract_markers_avx2(record.seq(), &mut vec, args.c, args.k);
                    }
                }

                return vec;
            },
            |record_sets| {
                while let Some(result) = record_sets.next() {
                    let (_rec_set, vec) = result.unwrap();
                    for km in vec {
                        let c = read_sketch.kmer_counts.entry(km).or_insert(0);
                        *c += 1;
                    }
                }
            },
        );
    }

    return read_sketch;
}

pub fn sketch_references(args: &Sketch) -> GenomesSketch {
    let genome_files;
    let mut temp_vec = vec![];
    if !args.reference_list.is_none() {
        let file = File::open(&args.reference_list.as_ref().unwrap()).unwrap();
        let reader = BufReader::new(file);
        for line in reader.lines() {
            temp_vec.push(line.unwrap().trim().to_string());
        }
        genome_files = &temp_vec;
    } else {
        genome_files = args.references.as_ref().unwrap();
    }
    let kmer_maps: Mutex<Vec<_>> = Mutex::new(vec![IntSet::default(); genome_files.len()]);
    let index_vec = (0..genome_files.len()).collect::<Vec<usize>>();
    let count: Mutex<usize> = Mutex::new(0);
    index_vec.into_par_iter().for_each(|i| {
        let ref_file = &genome_files[i];
        let reader = parse_fastx_file(&ref_file);
        let mut vec = vec![];
        if !reader.is_ok() {
            warn!("{} is not a valid fasta/fastq file; skipping.", ref_file);
        } else {
            let mut reader = reader.unwrap();
            while let Some(record) = reader.next() {
                if record.is_ok() {
                    let record = record.expect(&format!("Invalid record for file {}", ref_file));
                    let seq = record.seq();
                    unsafe {
                        extract_markers_avx2(&seq, &mut vec, args.c, args.k);
                    }
                } else {
                    warn!("File {} is not a valid fasta/fastq file", ref_file);
                }
            }
            let mut kmer_set = IntSet::default();
            for km in vec {
                kmer_set.insert(km);
            }
            let mut locked = kmer_maps.lock().unwrap();
            locked[i] = kmer_set;
            let mut c = count.lock().unwrap();
            *c += 1;
            if *c % 100 == 0 && *c > 0 {
                info!("{} references sketched", *c);
            }
        }
    });

    let gn_sketch = GenomesSketch {
        genome_kmer_counts: kmer_maps.into_inner().unwrap(),
        file_names: genome_files.clone(),
        c: args.c,
        k: args.k,
    };

    let ret;
    if args.derep < 1.00 && args.derep > 0.00 {
        ret = derep_genome_sketch(gn_sketch, args.derep, args.c);
    } else {
        ret = gn_sketch;
    }

    return ret;
}

fn derep_genome_sketch(genome_sketch: GenomesSketch, derep_val: f64, c: usize) -> GenomesSketch {
    info!("Dereplicating...");
    let mut rep_genomes: MMHashSet<_> = MMHashSet::default();
    let mut kmer_to_genome_table: IntMap<u64, Vec<_>> = IntMap::default();
    let num_gn = genome_sketch.file_names.len();
    let genome_kmers = &genome_sketch.genome_kmer_counts;
    let thresh = u64::MAX / 1000;

    for (i, kmer_map) in genome_kmers.iter().enumerate() {
        for kmer in kmer_map.iter(){
            if *kmer < thresh{
                let genomes = kmer_to_genome_table.entry(*kmer).or_insert(vec![]);
                genomes.push(i);
            }
        }
    }

    info!("Dereplication - finished populating k-mer to genome index");

    let mut used_genomes: MMHashSet<usize> = MMHashSet::default();

    for i in 0..genome_kmers.len() {
        if used_genomes.contains(&i) {
            continue;
        }
        rep_genomes.insert(i);
        let mut genome_containment_count: IntMap<_, _> = IntMap::default();
        let kmers = genome_kmers[i].iter();
        for kmer in kmers {
            if *kmer < thresh{
                let corresp_genomes = &kmer_to_genome_table[kmer];
                for j in corresp_genomes {
                    if used_genomes.contains(j) {
                        continue;
                    }
                    let c = genome_containment_count.entry(*j).or_insert(0);
                    *c += 1;
                }
            }
        }

        for (j, count) in genome_containment_count {
            let mut ratio = count as f64
//                    / ((genome_kmers[i].len() as f64 + genome_kmers[j].len() as f64)/2.);
                    / f64::min(genome_kmers[i].len() as f64 , genome_kmers[j].len() as f64);
            if c < 1000{
                ratio = ratio / c as f64 * 1000.;
            }
            let ani = f64::powf(ratio, 1. / genome_sketch.k as f64);
            if ani > derep_val {
                used_genomes.insert(j);
            }
        }
    }

    let ret_kmer_map: Vec<IntSet<u64>> = genome_sketch
        .genome_kmer_counts
        .into_iter()
        .enumerate()
        .filter(|(i, _)| rep_genomes.contains(i))
        .map(|(_, x)| x)
        .collect();
    let ret_file_names = genome_sketch
        .file_names
        .into_iter()
        .enumerate()
        .filter(|(i, _)| rep_genomes.contains(i))
        .map(|(_, x)| x)
        .collect();

    info!(
        "Genomes before dereplication: {}. Genomes after dereplicaction: {}",
        num_gn,
        ret_kmer_map.len()
    );

    return GenomesSketch {
        genome_kmer_counts: ret_kmer_map,
        file_names: ret_file_names,
        c: genome_sketch.c,
        k: genome_sketch.k,
    };
}
