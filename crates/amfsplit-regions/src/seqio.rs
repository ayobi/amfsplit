//! Reads, in the prototype's way: FASTA or FASTQ, gzipped when the path ends in `.gz`.
//!
//! As in `read_sequences()` of the prototype: the read ID is the first whitespace-separated
//! token of the header, FASTA sequence lines are stripped and joined, FASTQ records are four
//! lines, and sequences are upper-cased. Lines are split on LF; CR is stripped as whitespace.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;

use crate::pycompat::{py_first_token, py_strip};

/// Open a file for reading, transparently decompressing it when the name ends in `.gz`.
pub fn open_maybe_gz(path: &Path) -> Result<Box<dyn BufRead>> {
    let f = File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    let is_gz = path.to_string_lossy().ends_with(".gz");
    let r: Box<dyn Read> = if is_gz {
        Box::new(MultiGzDecoder::new(f))
    } else {
        Box::new(f)
    };
    Ok(Box::new(BufReader::with_capacity(1 << 20, r)))
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Fasta,
    Fastq,
}

/// Streaming reader yielding `(read_id, upper-cased sequence)`.
pub struct SeqReader {
    inner: Box<dyn BufRead>,
    kind: Option<Kind>,
    pending_header: Option<Vec<u8>>,
    done: bool,
    path: String,
}

impl SeqReader {
    pub fn open(path: &Path) -> Result<Self> {
        let mut inner = open_maybe_gz(path)?;
        let mut first = Vec::new();
        inner.read_until(b'\n', &mut first)?;
        let (kind, done) = match first.first() {
            None => (None, true),
            Some(b'@') => (Some(Kind::Fastq), false),
            Some(b'>') => (Some(Kind::Fasta), false),
            _ => bail!("{} is neither FASTA nor FASTQ", path.display()),
        };
        Ok(SeqReader {
            inner,
            kind,
            pending_header: if done { None } else { Some(first) },
            done,
            path: path.display().to_string(),
        })
    }

    fn line(&mut self) -> Result<Option<Vec<u8>>> {
        let mut buf = Vec::new();
        let n = self.inner.read_until(b'\n', &mut buf)?;
        Ok(if n == 0 { None } else { Some(buf) })
    }

    fn id_of(&self, header: &[u8]) -> Result<String> {
        match py_first_token(&header[1..]) {
            Some(tok) => Ok(String::from_utf8_lossy(tok).into_owned()),
            None => bail!("{}: empty header line", self.path),
        }
    }

    fn next_record(&mut self) -> Result<Option<(String, Vec<u8>)>> {
        if self.done {
            return Ok(None);
        }
        let header = match self.pending_header.take() {
            Some(h) => h,
            None => {
                self.done = true;
                return Ok(None);
            }
        };
        let id = self.id_of(&header)?;
        match self.kind.expect("kind set when not done") {
            Kind::Fastq => {
                let seq = self.line()?.unwrap_or_default();
                let mut seq = py_strip(&seq).to_vec();
                seq.make_ascii_uppercase();
                let _plus = self.line()?;
                let _qual = self.line()?;
                match self.line()? {
                    Some(h) => self.pending_header = Some(h),
                    None => self.done = true,
                }
                Ok(Some((id, seq)))
            }
            Kind::Fasta => {
                let mut seq = Vec::new();
                loop {
                    match self.line()? {
                        None => {
                            self.done = true;
                            break;
                        }
                        Some(l) if l.first() == Some(&b'>') => {
                            self.pending_header = Some(l);
                            break;
                        }
                        Some(l) => seq.extend_from_slice(py_strip(&l)),
                    }
                }
                seq.make_ascii_uppercase();
                Ok(Some((id, seq)))
            }
        }
    }
}

impl Iterator for SeqReader {
    type Item = Result<(String, Vec<u8>)>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.next_record() {
            Ok(Some(r)) => Some(Ok(r)),
            Ok(None) => None,
            Err(e) => {
                self.done = true;
                Some(Err(e))
            }
        }
    }
}

/// IUPAC complement as in the prototype's `revcomp()`; other characters are kept.
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'D' => b'H',
        b'H' => b'D',
        b'V' => b'B',
        b'N' => b'N',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        b'r' => b'y',
        b'y' => b'r',
        b's' => b's',
        b'w' => b'w',
        b'k' => b'm',
        b'm' => b'k',
        b'b' => b'v',
        b'd' => b'h',
        b'h' => b'd',
        b'v' => b'b',
        b'n' => b'n',
        other => other,
    }
}

/// Reverse complement (IUPAC-aware).
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn tmp(name: &str, content: &[u8]) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("amfsplit-seqio-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join(name);
        std::fs::File::create(&p)
            .unwrap()
            .write_all(content)
            .unwrap();
        p
    }

    #[test]
    fn fasta_multiline_lowercase_and_empty_record() {
        let p = tmp("a.fa", b">r1 desc\nacg\nTN\r\n>r2\n\n>r3\nAC");
        let recs: Vec<_> = SeqReader::open(&p).unwrap().map(|r| r.unwrap()).collect();
        assert_eq!(
            recs,
            vec![
                ("r1".into(), b"ACGTN".to_vec()),
                ("r2".into(), b"".to_vec()),
                ("r3".into(), b"AC".to_vec()),
            ]
        );
    }

    #[test]
    fn fastq_records() {
        let p = tmp("a.fq", b"@q1 x\nacgt\n+\nIIII\n@q2\nGG\n+\nII\n");
        let recs: Vec<_> = SeqReader::open(&p).unwrap().map(|r| r.unwrap()).collect();
        assert_eq!(
            recs,
            vec![
                ("q1".into(), b"ACGT".to_vec()),
                ("q2".into(), b"GG".to_vec())
            ]
        );
    }

    #[test]
    fn revcomp_iupac() {
        assert_eq!(revcomp(b"ACGTRYKMBDHVNX"), b"XNBDHVKMRYACGT".to_vec());
    }
}
