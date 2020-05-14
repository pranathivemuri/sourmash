//! # Compressed representations of genomic data
//!
//! A signature is a collection of sketches for a genomic dataset.

use std::fs::File;
use std::io;
use std::iter::Iterator;
use std::path::Path;
use std::str;

use cfg_if::cfg_if;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use typed_builder::TypedBuilder;

#[cfg(all(target_arch = "wasm32", target_vendor = "unknown"))]
use wasm_bindgen::prelude::*;

use crate::index::storage::ToWriter;
use crate::sketch::minhash::HashFunctions;
use crate::sketch::Sketch;
use crate::Error;

#[derive(Serialize, Deserialize, Debug, Clone, TypedBuilder)]
pub struct Signature {
    #[serde(default = "default_class")]
    #[builder(default = default_class())]
    class: String,

    #[serde(default)]
    #[builder(default)]
    email: String,

    #[builder(setter(into))]
    hash_function: String,

    #[builder(default)]
    filename: Option<String>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) name: Option<String>,

    #[serde(default = "default_license")]
    #[builder(default = default_license())]
    license: String,

    pub(crate) signatures: Vec<Box<dyn Sketch>>,

    #[serde(default = "default_version")]
    #[builder(default = default_version())]
    version: f64,
}

fn default_license() -> String {
    "CC0".to_string()
}

fn default_class() -> String {
    "sourmash_signature".to_string()
}

fn default_version() -> f64 {
    0.4
}

impl Signature {
    pub fn name(&self) -> String {
        if let Some(name) = &self.name {
            name.clone()
        } else if let Some(filename) = &self.filename {
            filename.clone()
        } else {
            self.md5sum()
        }
    }

    pub fn set_name(&mut self, name: &str) {
        self.name = Some(name.into())
    }

    pub fn filename(&self) -> String {
        if let Some(filename) = &self.filename {
            filename.clone()
        } else {
            "".into()
        }
    }

    pub fn set_filename(&mut self, name: &str) {
        self.filename = Some(name.into())
    }

    pub fn size(&self) -> usize {
        self.signatures.len()
    }

    pub fn sketches(&self) -> Vec<Box<dyn Sketch>> {
        self.signatures.clone()
    }

    pub fn reset_sketches(&mut self) {
        self.signatures = vec![];
    }

    pub fn push(&mut self, sketch: Box<dyn Sketch>) {
        self.signatures.push(sketch);
    }

    pub fn license(&self) -> String {
        self.license.clone()
    }

    pub fn class(&self) -> String {
        self.class.clone()
    }

    pub fn hash_function(&self) -> String {
        self.hash_function.clone()
    }

    pub fn email(&self) -> String {
        self.email.clone()
    }

    pub fn md5sum(&self) -> String {
        if self.signatures.len() == 1 {
            self.signatures[0].md5sum()
        } else {
            // TODO: select the correct signature
            unimplemented!()
        }
    }

    pub fn select_sketch(&self, template: &Box<dyn Sketch>) -> Option<&Box<dyn Sketch>> {
        for sk in &self.signatures {
            if sk.check_compatible(template).is_ok() {
                return Some(sk);
            }
        }
        None
    }

    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Vec<Signature>, Error> {
        let mut reader = io::BufReader::new(File::open(path)?);
        Ok(Signature::from_reader(&mut reader)?)
    }

    pub fn from_reader<R>(rdr: R) -> Result<Vec<Signature>, Error>
    where
        R: io::Read,
    {
        let (rdr, _format) = niffler::get_reader(Box::new(rdr))?;

        let sigs: Vec<Signature> = serde_json::from_reader(rdr)?;
        Ok(sigs)
    }

    pub fn load_signatures<R>(
        buf: R,
        ksize: Option<usize>,
        moltype: Option<HashFunctions>,
        _scaled: Option<u64>,
    ) -> Result<Vec<Signature>, Error>
    where
        R: io::Read,
    {
        let orig_sigs = Signature::from_reader(buf)?;

        let flat_sigs = orig_sigs.into_iter().flat_map(|s| {
            s.signatures
                .iter()
                .map(|mh| {
                    let mut new_s = s.clone();
                    new_s.signatures = vec![(*mh).clone()];
                    new_s
                })
                .collect::<Vec<Signature>>()
        });

        let filtered_sigs = flat_sigs.filter_map(|mut sig| {
            let good_mhs: Vec<Box<dyn Sketch>> = sig
                .signatures
                .into_iter()
                .filter(|sketch| {
                    if let Some(k) = ksize {
                        if k != sketch.ksize() as usize {
                            return false;
                        }
                    };

                    match moltype {
                        Some(x) => {
                            if sketch.hash_function() == x {
                                return true;
                            }
                        }
                        None => return true, // TODO: match previous behavior
                    };
                    false
                })
                .collect();

            if good_mhs.is_empty() {
                return None;
            };

            sig.signatures = good_mhs;
            Some(sig)
        });

        Ok(filtered_sigs.collect())
    }

    pub fn add_sequence(&mut self, seq: &[u8], force: bool) -> Result<(), Error> {
        cfg_if! {
        if #[cfg(feature = "parallel")] {
            self.signatures
                .par_iter_mut()
                .for_each(|sketch| {
                    sketch.add_sequence(&seq, force).unwrap(); }
                );
        } else {
            self.signatures
                .iter_mut()
                .for_each(|sketch| {
                    sketch.add_sequence(&seq, force).unwrap(); }
                );
        }
        }

        Ok(())
    }

    pub fn add_protein(&mut self, seq: &[u8]) -> Result<(), Error> {
        cfg_if! {
        if #[cfg(feature = "parallel")] {
            self.signatures
                .par_iter_mut()
                .try_for_each(|sketch| {
                    sketch.add_protein(&seq) }
                )?;
        } else {
            self.signatures
                .iter_mut()
                .try_for_each(|sketch| {
                    sketch.add_protein(&seq) }
                )?;
        }
        }

        Ok(())
    }
}

impl ToWriter for Signature {
    fn to_writer<W>(&self, writer: &mut W) -> Result<(), Error>
    where
        W: io::Write,
    {
        serde_json::to_writer(writer, &vec![&self])?;
        Ok(())
    }
}

impl Default for Signature {
    fn default() -> Signature {
        Signature {
            class: default_class(),
            email: "".to_string(),
            hash_function: "0.murmur64".to_string(),
            license: default_license(),
            filename: None,
            name: None,
            signatures: vec![],
            version: default_version(),
        }
    }
}

impl PartialEq for Signature {
    fn eq(&self, other: &Signature) -> bool {
        let metadata = self.class == other.class
            && self.email == other.email
            && self.hash_function == other.hash_function
            && self.filename == other.filename
            && self.name == other.name;

        // TODO: find the right signature
        // as long as we have a matching
        let sk = &self.signatures[0];
        let other_sk = &other.signatures[0];
        return metadata && (sk == other_sk);
    }
}

#[cfg(test)]
mod test {
    use std::convert::TryInto;
    use std::fs::File;
    use std::io::{BufReader, Read};
    use std::path::PathBuf;

    use needletail::parse_fastx_reader;

    use crate::cmd::ComputeParameters;

    use super::Signature;

    #[test]
    fn load_sig() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/.sbt.v3/60f7e23c24a8d94791cc7a8680c493f9");

        let mut reader = BufReader::new(File::open(filename).unwrap());
        let sigs = Signature::load_signatures(
            &mut reader,
            Some(31),
            Some("DNA".try_into().unwrap()),
            None,
        )
        .unwrap();
        let _sig_data = sigs[0].clone();
        // TODO: check sig_data
    }

    #[test]
    fn load_signature() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/genome-s10+s11.sig");

        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let sigs: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");

        assert_eq!(sigs.len(), 1);

        let sig = sigs.get(0).unwrap();
        assert_eq!(sig.class, "sourmash_signature");
        assert_eq!(sig.email, "");
        if let Some(ref filename) = sig.filename {
            assert_eq!(filename, "-");
        }
        assert_eq!(sig.hash_function, "0.murmur64");
        if let Some(ref name) = sig.name {
            assert_eq!(name, "s10+s11");
        }
        assert_eq!(sig.signatures.len(), 4);
    }

    #[test]
    fn signature_from_computeparams() {
        let params = ComputeParameters::builder()
            .ksizes(vec![2, 3, 4])
            .num_hashes(3u32)
            .build();

        let mut sig = Signature::from_params(&params);
        sig.add_sequence(b"ATGC", false).unwrap();

        assert_eq!(sig.signatures.len(), 3);
        dbg!(&sig.signatures);
        assert_eq!(sig.signatures[0].size(), 3);
        assert_eq!(sig.signatures[1].size(), 2);
        assert_eq!(sig.signatures[2].size(), 1);
    }

    #[test]
    fn signature_slow_path() {
        let params = ComputeParameters::builder()
            .ksizes(vec![2, 3, 4, 5])
            .num_hashes(3u32)
            .build();

        let mut sig = Signature::from_params(&params);
        sig.add_sequence(b"ATGCTN", true).unwrap();

        assert_eq!(sig.signatures.len(), 4);
        dbg!(&sig.signatures);
        assert_eq!(sig.signatures[0].size(), 3);
        assert_eq!(sig.signatures[1].size(), 3);
        assert_eq!(sig.signatures[2].size(), 2);
        assert_eq!(sig.signatures[3].size(), 1);
    }

    #[test]
    fn signature_add_sequence_protein() {
        let params = ComputeParameters::builder()
            .ksizes(vec![3, 6])
            .num_hashes(3u32)
            .protein(true)
            .dna(false)
            .build();

        let mut sig = Signature::from_params(&params);
        sig.add_sequence(b"ATGCAT", false).unwrap();

        assert_eq!(sig.signatures.len(), 2);
        dbg!(&sig.signatures);
        assert_eq!(sig.signatures[0].size(), 3);
        assert_eq!(sig.signatures[1].size(), 1);
    }

    #[test]
    fn signature_add_protein() {
        let params = ComputeParameters::builder()
            .ksizes(vec![3, 6])
            .num_hashes(3u32)
            .protein(true)
            .dna(false)
            .build();

        let mut sig = Signature::from_params(&params);
        sig.add_protein(b"AGY").unwrap();

        assert_eq!(sig.signatures.len(), 2);
        dbg!(&sig.signatures);
        assert_eq!(sig.signatures[0].size(), 3);
        assert_eq!(sig.signatures[1].size(), 2);
    }

    #[test]
    fn signature_add_sequence_cp() {
        let mut cp = ComputeParameters::default();
        cp.set_dayhoff(true);
        cp.set_protein(true);
        cp.set_hp(true);
        cp.set_dna(true);

        let mut sig = Signature::from_params(&cp);

        let mut data: Vec<u8> = vec![];
        let mut f = File::open("../../tests/test-data/ecoli.genes.fna").unwrap();
        let _ = f.read_to_end(&mut data);

        let mut parser = parse_fastx_reader(&data[..]).unwrap();
        while let Some(record) = parser.next() {
            let record = record.unwrap();
            sig.add_sequence(&record.seq(), false).unwrap();
        }

        let sketches = sig.sketches();
        assert_eq!(sketches.len(), 12);
        for sk in sketches {
            assert_eq!(sk.size(), 500);
        }
    }
}
