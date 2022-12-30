use std::collections::HashMap;
use rkyv::{Archive, Deserialize, Serialize};
// bytecheck can be used to validate your data if you want
use bytecheck::CheckBytes;
use nohash_hasher::IntMap;
use nohash_hasher::IntSet;
use std::hash::{BuildHasherDefault, Hash, Hasher};
use rustc_hash::FxHashMap;
use nohash_hasher::BuildNoHashHasher;
use std::collections::HashSet;

pub type Kmer = u64;
pub const BYTE_TO_SEQ: [u64; 256] = [
    0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

#[inline]
pub fn mm_hash(bytes: &[u8]) -> usize {
    let mut key = usize::from_ne_bytes(bytes.try_into().unwrap()) as usize;
    key = !key.wrapping_add(key << 21); // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8); // key * 265
    key = key ^ key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4); // key * 21
    key = key ^ key >> 28;
    key = key.wrapping_add(key << 31);
    return key;
}

pub struct MMHasher {
    hash: usize,
}

impl Hasher for MMHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        self.hash = mm_hash(bytes);
    }
    #[inline]
    fn finish(&self) -> u64 {
        self.hash as u64
    }
}

impl Default for MMHasher {
    #[inline]
    fn default() -> MMHasher {
        MMHasher { hash: 0 }
    }
}

//Implement minimap2 hashing, will test later.
pub type MMBuildHasher = BuildHasherDefault<MMHasher>;
pub type MMHashMap<K, V> = HashMap<K, V, MMBuildHasher>;
pub type MMHashSet<K> = HashSet<K, MMBuildHasher>;

#[derive(Archive, Default, Deserialize, Serialize, Debug, PartialEq)]
pub struct SequencesSketch{
//    pub kmer_counts: HashMap<Kmer, u32,BuildNoHashHasher<Kmer>>,
    pub kmer_counts: HashMap<Kmer, u32>,
//    pub kmer_counts: FxHashMap<Kmer, u32>,
    pub c: usize,
    pub k: usize,
    pub file_name: String
}

impl SequencesSketch{
    pub fn new(file_name: String, c: usize, k: usize) -> SequencesSketch{
        return SequencesSketch{kmer_counts : HashMap::default(), file_name, c, k}
    }
}

#[derive(Archive, Deserialize, Serialize, Debug, PartialEq)]
// This will generate a PartialEq impl between our unarchived and archived types
// To use the safe API, you have to derive CheckBytes for the archived type
#[derive(Default)]
pub struct GenomesSketch{
    pub genome_kmer_counts: Vec<IntSet<Kmer>>,
    pub file_names: Vec<String>,
    pub c: usize,
    pub k: usize,
}

impl GenomesSketch{
    pub fn new(file_names: Vec<String>, c:usize, k: usize) -> GenomesSketch{
        return GenomesSketch{genome_kmer_counts : vec![], file_names, c, k}
    }
}

