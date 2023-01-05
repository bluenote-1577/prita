use crate::cmdline::*;
use crate::types::*;
use log::*;
use rayon::prelude::*;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::read;
use std::fs::File;
use std::io::BufReader;
use std::sync::Mutex;

pub fn load(args: Load) {
    let genome_file =
        File::open(args.references_sketch.clone()).expect("Genome sketch not a valid file");
    let genome_reader = BufReader::with_capacity(10_000_000, genome_file);
    let genome_sketch: GenomesSketch = bincode::deserialize_from(genome_reader).unwrap();
    info!("Genome sketch loading complete.");

    let mut results = vec![];
    for query_file in args.query_sketch.iter() {
        let res_vec = work(&args, query_file, &genome_sketch);
        results.push(res_vec);
    }

    let mut ref_gn_set : HashSet<_> = HashSet::default();
    for res_vec in results.iter(){
        for res in res_vec.iter(){
            ref_gn_set.insert(res.file_name_ref_gn);
        }
    }

    let mut ref_gn_vec = ref_gn_set.into_iter().collect::<Vec<&str>>();
    ref_gn_vec.sort();
    let mut ref_to_index_map = HashMap::new();
    for (i, ref_gn) in ref_gn_vec.iter().enumerate(){
        ref_to_index_map.insert(*ref_gn,i);
    }

    let mut ani_matrix_rows = vec![];
    let mut cov_matrix_rows = vec![];
    for res_vec in results.iter(){
//        dbg!(&res_vec);
        let mut ani_row = vec![0.; ref_to_index_map.len()];
        let mut cov_row = vec![0.; ref_to_index_map.len()];
        let cov_sum = res_vec.iter().map(|x| x.cov).sum::<f32>();
        for res in res_vec{
            let ind = ref_to_index_map[res.file_name_ref_gn];
            ani_row[ind] = res.ani;
            cov_row[ind] = res.cov/cov_sum;
        }
        ani_matrix_rows.push(ani_row);
        cov_matrix_rows.push(cov_row);
    }
    if ref_gn_vec.len() == 0{
        return
    }
    print!("query");
    for i in 0..ref_gn_vec.len(){
        print!("\t{}", &ref_gn_vec[i]);
    }
    println!("");
    for (j,row) in cov_matrix_rows.into_iter().enumerate(){
        if row.len() == 0{
            return;
        }
        print!("{}\t{}", &args.query_sketch[j], row[0]);
        for i in 1..row.len(){
            print!("\t{}", row[i]);
        }
        println!("");
    }
}

fn work<'a>(
    args: &Load,
    query_file: &str,
    genome_sketch: &'a GenomesSketch,
) -> Vec<LoadResult<'a>> {
    let sequence_file = File::open(query_file).expect("Sequence sketch not a valid file");
    let seq_reader = BufReader::with_capacity(10_000_000, sequence_file);
    let sequence_sketch_enc: SequencesSketchEncode = bincode::deserialize_from(seq_reader).unwrap();
    let sequence_sketch = SequencesSketch::from_enc(sequence_sketch_enc);
    info!("Sequence sketch loading complete.");

    if genome_sketch.k != sequence_sketch.k {
        panic!(
            "k parameter for query {} != k parameter for reference {}",
            sequence_sketch.k, genome_sketch.k
        )
    }
    if genome_sketch.c != sequence_sketch.c {
        panic!(
            "c parameter for query {} != c parameter for reference {}",
            sequence_sketch.c, genome_sketch.c
        )
    }
    //    let used_kmersli : Mutex<HashSet<u64>> = Mutex::new(HashSet::default(););
    let ani_desc: Mutex<Vec<_>> = Mutex::new(vec![]);
    let anis: Mutex<MMHashMap<_, _>> = Mutex::new(MMHashMap::default());
    let useable_refs: Mutex<Vec<_>> = Mutex::new(vec![]);
    let used_kmers: Mutex<MMHashSet<_>> = Mutex::new(MMHashSet::default());
    let index_vec = (0..genome_sketch.file_names.len()).collect::<Vec<usize>>();
    let progress: Mutex<usize> = Mutex::new(0);
    index_vec.into_par_iter().for_each(|i| {
        let mut count = 0;
        let mut covs = vec![];
        let gn_kmers = &genome_sketch.genome_kmer_counts[i];
        for kmer in gn_kmers.iter() {
            if sequence_sketch.kmer_counts.contains_key(kmer) {
                used_kmers.lock().unwrap().insert(kmer);
                count += 1;
                covs.push(sequence_sketch.kmer_counts[kmer]);
            }
        }
        //        dbg!(count,gn_kmers.len(),sequence_sketch.kmer_counts.len());
        let ani = f64::powf(
            count as f64 / gn_kmers.len() as f64,
            1. / genome_sketch.k as f64,
        );
        covs.sort();
        //        dbg!(&covs);

        {
            let mut num_proc = progress.lock().unwrap();
            *num_proc += 1;
            if *num_proc % 1000 == 0 && *num_proc > 0 {
                info!("{} references processed", num_proc);
            }
        }
        if ani > args.min_ani {
            let ani = correct_ani_abund(ani, &covs, sequence_sketch.k as f64);
            let cov = covs
                .iter()
                .enumerate()
                .filter(|(i, _)| *i > covs.len() / 10 && *i * 10 < covs.len() * 9)
                .map(|(_, x)| *x)
                .sum::<u32>() as f32
                / (covs.len() as f32)
                * 0.9;
            anis.lock()
                .unwrap()
                .insert(i, (ani, cov, &genome_sketch.file_names[i]));
            useable_refs.lock().unwrap().push(i);
            ani_desc.lock().unwrap().push((
                ani,
                covs[covs.len() / 2],
                cov,
                genome_sketch.file_names[i].clone(),
            ));
        }
    });

    let mut used_kmer_c = 0;
    for kmer in used_kmers.into_inner().unwrap() {
        used_kmer_c += sequence_sketch.kmer_counts[kmer];
    }
    let mut all_kmers_sequence = 0;
    for val in sequence_sketch.kmer_counts.values() {
        all_kmers_sequence += val;
    }
    info!("Number of kmers found in references passing filter: {}. Number of total kmers in query: {}. Fraction: {}", used_kmer_c, all_kmers_sequence, used_kmer_c as f64 / all_kmers_sequence as f64);

    let mut ani_desc = ani_desc.into_inner().unwrap();
    let useable_refs = useable_refs.into_inner().unwrap();
    let anis = anis.into_inner().unwrap();
    ani_desc.sort_by(|x, y| x.partial_cmp(&y).unwrap());
    for ani in ani_desc {
        debug!("{}\t{}\t{}\t{}", ani.0, ani.1, ani.2, ani.3);
    }

    let mut ret_vec = vec![];
    let mut clusters = cluster_genomes(&genome_sketch, &useable_refs, args.cluster);
    for cluster in clusters.iter_mut() {
        cluster.sort();
        debug!("Cluster {:?}", &cluster);
        if cluster.len() > 1 {
            let anis_cluster: Vec<f64> = cluster.iter().map(|x| anis[x].0).collect();
            //        dbg!(&cluster);
            let res = em(&sequence_sketch, &genome_sketch, &cluster, &anis_cluster);
            ret_vec.extend(res);
        } else {
            let ani_res = &anis[&cluster[0]];
            let res = LoadResult {
                ani: ani_res.0 as f32,
                cov: ani_res.1 as f32,
                em_abund: 1.,
                file_name_ref_gn: ani_res.2,
            };
            ret_vec.push(res)
        }
    }
    return ret_vec;
}

fn cluster_genomes(
    genome_sketch: &GenomesSketch,
    useable_refs: &Vec<usize>,
    cluster_val: f64,
) -> Vec<Vec<usize>> {
    let mut clusters = vec![];
    let mut used_genomes: MMHashSet<_> = MMHashSet::default();
    let mut kmer_to_genome_table: MMHashMap<u64, Vec<_>> = MMHashMap::default();
    let genome_kmers = &genome_sketch.genome_kmer_counts;

    for i in useable_refs.iter() {
        let kmer_map = &genome_kmers[*i];
        //This may occur for eukaryotic contaminiation
        if kmer_map.len() > 10_000_000 {
            continue;
        }
        for kmer in kmer_map.iter() {
            let genomes = kmer_to_genome_table.entry(*kmer).or_insert(vec![]);
            genomes.push(*i);
        }
    }

    for i in useable_refs.iter() {
        if used_genomes.contains(i) {
            continue;
        }
        used_genomes.insert(*i);

        let mut cluster = vec![];
        cluster.push(*i);

        let mut genome_containment_count: MMHashMap<_, _> = MMHashMap::default();

        let kmers = genome_kmers[*i].iter();
        if kmers.len() > 10_000_00 {
            continue;
        }
        for kmer in kmers {
            let corresp_genomes = &kmer_to_genome_table[kmer];
            for j in corresp_genomes {
                if !used_genomes.contains(j) {
                    let c = genome_containment_count.entry(*j).or_insert(0);
                    *c += 1;
                }
            }
        }

        for (j, count) in genome_containment_count {
            let ratio = count as f64
                / f64::min(genome_kmers[*i].len() as f64, genome_kmers[j].len() as f64);
            let ani = f64::powf(ratio, 1. / genome_sketch.k as f64);
            if ani > cluster_val {
                cluster.push(j);
                used_genomes.insert(j);
            }
        }

        clusters.push(cluster);
    }
    return clusters;
}

pub fn em<'a>(
    sequence_sketch: &SequencesSketch,
    genome_sketch: &'a GenomesSketch,
    useable_refs: &Vec<usize>,
    anis: &Vec<f64>,
) -> Vec<LoadResult<'a>> {
    let thresh = 100000000;
    let anis_k: Vec<f64> = anis
        .iter()
        .map(|x| f64::powf(*x, sequence_sketch.k as f64))
        .collect();
    let mut genome_kmer_counts_in_sequence = vec![0.; useable_refs.len()];
    let full_lengths: Vec<f64> = genome_sketch
        .genome_kmer_counts
        .iter()
        .enumerate()
        .filter(|(i, _)| useable_refs.contains(i))
        .map(|(_, x)| x.len() as f64)
        .collect();
    let mut kmer_union: MMHashSet<_> = MMHashSet::default();
    let mut equiv_classes: HashMap<Vec<usize>, usize> = HashMap::default();
    let mut kmer_to_refs: HashMap<u64, Vec<usize>> = HashMap::default();
    for (iter_num, i) in useable_refs.iter().enumerate() {
        let kmer_map = &genome_sketch.genome_kmer_counts[*i];
        let mut useable_count = 0;
        for kmer in kmer_map.iter() {
            if sequence_sketch.kmer_counts.contains_key(kmer) {
                useable_count += 1;
                if sequence_sketch.kmer_counts[kmer] > thresh {
                    continue;
                }
                kmer_union.insert(kmer);
                let vec = kmer_to_refs.entry(*kmer).or_insert(vec![]);
                vec.push(iter_num);
            }
        }
        genome_kmer_counts_in_sequence[iter_num] = useable_count as f64;
    }
    let lengths = genome_kmer_counts_in_sequence;

    for (kmer, class) in kmer_to_refs {
        let c = equiv_classes.entry(class).or_insert(0);
        *c += sequence_sketch.kmer_counts[&kmer] as usize;
    }

    let mut total_kmers = 0.;
    for kmer in kmer_union.iter() {
        total_kmers += sequence_sketch.kmer_counts[kmer] as f64;
    }
    //    dbg!(useable_refs.len());
    //    dbg!(kmer_union.len(), genome_sketch.genome_kmer_counts[useable_refs[0]].len());

    let anis_k_sum: f64 = anis_k.iter().sum();
    let mut thetas: Vec<f64> = anis_k.iter().map(|x| x / anis_k_sum).collect();
    //    let mut thetas = vec![1. / useable_refs.len() as f64; useable_refs.len()];
    for _ in 0..100 {
        let mut new_thetas = vec![0.; useable_refs.len()];
        for (class, count) in equiv_classes.iter() {
            let mut denom = 0.;
            for j in class.iter() {
                denom += 1. / full_lengths[*j] * thetas[*j] * anis_k[*j];
                //                denom += 1. / lengths[*j] * thetas[*j];
            }
            for i in class.iter() {
                let num = (*count) as f64 / full_lengths[*i] * thetas[*i] * anis_k[*i];
                //                let num = (*count)  as f64 / lengths[*i] * thetas[*i];
                new_thetas[*i] += num / denom;
            }
        }
        for theta in new_thetas.iter_mut() {
            *theta /= total_kmers;
        }
        thetas = new_thetas;
    }
    let mut res = vec![];
    for i in 0..useable_refs.len() {
        let load_res = LoadResult {
            ani: anis[i] as f32,
            cov: (thetas[i] * total_kmers / lengths[i]) as f32,
            em_abund: thetas[i] as f32,
            file_name_ref_gn: &genome_sketch.file_names[useable_refs[i]],
        };
        res.push(load_res);
    }
    return res
}

fn correct_ani_abund(ani: f64, covs: &Vec<u32>, k: f64) -> f64 {
    let thresh = 5;
    let mut leq_threshold = 0;
    for cov in covs.iter() {
        if *cov < thresh {
            leq_threshold += 1;
        }
    }

    if (leq_threshold as f64) < 0.9 * covs.len() as f64 {
        return ani;
    }

    let mut counts: HashMap<_, _> = HashMap::default();
    for cov in covs.iter() {
        let c = counts.entry(cov).or_insert(0);
        *c += 1;
    }

    let mut lambda_ests = vec![];
    let mut weights = vec![];
    for i in 1..10 {
        if counts.contains_key(&i) && counts.contains_key(&(i + 1)) {
            lambda_ests.push(counts[&(i + 1)] as f64 / counts[&i] as f64 * (i + 1) as f64);
            weights.push(f64::min(counts[&(i + 1)] as f64, counts[&i] as f64));
        }
    }
    //    dbg!(&lambda_ests, &counts, &weights);
    let sum = weights.iter().sum::<f64>();
    if covs.len() < 30 || weights.len() == 0 {
        return ani;
    }

    let mut lambda = 0.;
    for i in 0..weights.len() {
        lambda += weights[i] * lambda_ests[i];
    }
    lambda /= sum;
    //    lambda = lambda_ests[0];
    let ret = f64::min(
        1.0,
        ani * f64::powf(1. / (1. - f64::powf(2.718, -lambda)), 1. / k),
    );
    trace!("{},{},{}", lambda, ret, ani);
    return ret;
}
