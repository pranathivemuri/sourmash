#[cfg(all(target_arch = "wasm32", target_vendor = "unknown"))]
use wasm_bindgen::prelude::*;

use getset::{CopyGetters, Getters, Setters};
use typed_builder::TypedBuilder;

use crate::index::MHBT;
use crate::signature::Signature;
use crate::sketch::minhash::{max_hash_for_scaled, HashFunctions, KmerMinHashBTree};
use crate::sketch::Sketch;
use crate::Error;

pub fn prepare(index_path: &str) -> Result<(), Error> {
    let mut index = MHBT::from_path(index_path)?;

    // TODO equivalent to fill_internal in python
    //unimplemented!();

    index.save_file(index_path, None)?;

    Ok(())
}

impl Signature {
    pub fn from_params(params: &ComputeParameters) -> Signature {
        let template = build_template(&params);

        Signature::builder()
            .hash_function("0.murmur64")
            .name(params.merge.clone())
            .filename(None)
            .signatures(template)
            .build()
    }
}

#[allow(dead_code)]
#[cfg_attr(all(target_arch = "wasm32", target_vendor = "unknown"), wasm_bindgen)]
#[derive(TypedBuilder, CopyGetters, Getters, Setters)]
pub struct ComputeParameters {
    #[getset(get = "pub", set = "pub")]
    #[builder(default = vec![21, 31, 51])]
    ksizes: Vec<u32>,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = false)]
    check_sequence: bool,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = true)]
    dna: bool,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = false)]
    dayhoff: bool,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = false)]
    hp: bool,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = false)]
    singleton: bool,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = 0usize)]
    count_valid_reads: usize,

    #[getset(get = "pub", set = "pub")]
    #[builder(default = None)]
    barcodes_file: Option<String>, // TODO: check

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = 1500usize)]
    line_count: usize,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = None)]
    rename_10x_barcodes: Option<bool>, // TODO: check

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = None)]
    write_barcode_meta_csv: Option<bool>, // TODO: check

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = None)]
    save_fastas: Option<bool>, // TODO: check

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = 0u64)]
    scaled: u64,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = false)]
    force: bool,

    #[getset(get = "pub", set = "pub")]
    #[builder(default = None)]
    output: Option<String>, // TODO: check

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = 500u32)]
    num_hashes: u32,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = false)]
    protein: bool,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = false)]
    name_from_first: bool,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = 42u64)]
    seed: u64,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = false)]
    input_is_protein: bool,

    #[getset(get = "pub", set = "pub")]
    #[builder(default = None)]
    merge: Option<String>,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = false)]
    track_abundance: bool,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = false)]
    randomize: bool,

    #[getset(get = "pub", set = "pub")]
    #[builder(default = "CC0".into())]
    license: String,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = false)]
    input_is_10x: bool,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = 2usize)]
    processes: usize,
}

impl Default for ComputeParameters {
    fn default() -> Self {
        Self::builder().build()
    }
}

pub fn build_template(params: &ComputeParameters) -> Vec<Box<dyn Sketch>> {
    let max_hash = max_hash_for_scaled(params.scaled).unwrap_or(0);

    params
        .ksizes
        .iter()
        .flat_map(|k| {
            let mut ksigs = vec![];

            if params.protein {
                ksigs.push(Box::new(
                    KmerMinHashBTree::builder()
                        .num(params.num_hashes)
                        .ksize(*k)
                        .hash_function(HashFunctions::murmur64_protein)
                        .max_hash(max_hash)
                        .seed(params.seed)
                        .abunds(if params.track_abundance {
                            Some(Default::default())
                        } else {
                            None
                        })
                        .build(),
                ) as _);
            }

            if params.dayhoff {
                ksigs.push(Box::new(
                    KmerMinHashBTree::builder()
                        .num(params.num_hashes)
                        .ksize(*k)
                        .hash_function(HashFunctions::murmur64_dayhoff)
                        .max_hash(max_hash)
                        .seed(params.seed)
                        .abunds(if params.track_abundance {
                            Some(Default::default())
                        } else {
                            None
                        })
                        .build(),
                ) as _);
            }

            if params.hp {
                ksigs.push(Box::new(
                    KmerMinHashBTree::builder()
                        .num(params.num_hashes)
                        .ksize(*k)
                        .hash_function(HashFunctions::murmur64_hp)
                        .max_hash(max_hash)
                        .seed(params.seed)
                        .abunds(if params.track_abundance {
                            Some(Default::default())
                        } else {
                            None
                        })
                        .build(),
                ) as _);
            }

            if params.dna {
                ksigs.push(Box::new(
                    KmerMinHashBTree::builder()
                        .num(params.num_hashes)
                        .ksize(*k)
                        .hash_function(HashFunctions::murmur64_DNA)
                        .max_hash(max_hash)
                        .seed(params.seed)
                        .abunds(if params.track_abundance {
                            Some(Default::default())
                        } else {
                            None
                        })
                        .build(),
                ) as _);
            }

            ksigs
        })
        .collect()
}
