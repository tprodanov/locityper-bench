use std::{
    io::{self, Write, BufWriter, Read, BufRead, BufReader},
    fs::File,
};
use flate2::{
    bufread::{GzDecoder, MultiGzDecoder},
    write::GzEncoder,
    Compression,
};

#[inline(always)]
fn read_line(mut f: impl BufRead, buf: &mut Vec<u8>) -> io::Result<usize> {
    f.read_until(b'\n', buf)
}

fn split_and_interleave(in_filename: &str, out_filename: &str) -> io::Result<()> {
    let mut f_in = BufReader::new(MultiGzDecoder::new(BufReader::new(File::open(in_filename)?)));
    let tmp_filename = format!("{}.tmp", out_filename);
    let mut f_out = BufWriter::with_capacity(131_072,
        GzEncoder::new(BufWriter::new(File::create(&tmp_filename)?), Compression::default()));

    let mut buf = Vec::with_capacity(1024);
    read_line(&mut f_in, &mut buf)?;
    read_line(&mut f_in, &mut buf)?;
    let mut count = 1;
    loop {
        buf.clear();
        let n = read_line(&mut f_in, &mut buf)?;
        if buf[n - 3] == b'.' && buf[n - 2] == b'1' {
            break;
        }
        read_line(&mut f_in, &mut buf)?;
        count += 1;
    }

    let old_buf = buf;
    let mut buf = Vec::with_capacity(1024);
    let mut f_in1 = BufReader::new(MultiGzDecoder::new(BufReader::new(File::open(in_filename)?)));
    read_line(&mut f_in1, &mut buf)?;
    read_line(&mut f_in1, &mut buf)?;
    buf.extend_from_slice(&old_buf);
    read_line(&mut f_in, &mut buf)?;
    f_out.write_all(&buf)?;

    for _ in 1..count {
        buf.clear();
        read_line(&mut f_in1, &mut buf)?;
        read_line(&mut f_in1, &mut buf)?;
        assert!(read_line(&mut f_in, &mut buf)? > 0);
        read_line(&mut f_in, &mut buf)?;
        f_out.write_all(&buf)?;
    }
    f_out.flush()?;
    std::mem::drop(f_out);

    let last_recs = String::from_utf8_lossy(&buf);
    let last_recs: Vec<_> = last_recs.split("\n").collect();
    assert_eq!(last_recs[0], last_recs[2]);
    buf.clear();
    buf.resize(16, 0);
    assert_eq!(f_in.read(&mut buf)?, 0);

    std::fs::rename(&tmp_filename, &out_filename)?;
    Ok(())
}

#[inline]
fn check_seq(seq: &[u8]) -> bool {
    seq.iter().all(|&b| b == b'A' || b == b'C' || b == b'G' || b == b'T' || b == b'N')
}

fn check(filename: &str) -> io::Result<()> {
    let mut f_in = BufReader::new(GzDecoder::new(BufReader::new(File::open(filename)?)));
    let mut buf1 = Vec::with_capacity(1024);
    let mut buf2 = Vec::with_capacity(1024);
    let mut buf3 = Vec::with_capacity(1024);
    let mut buf4 = Vec::with_capacity(1024);
    loop {
        buf1.clear();
        buf2.clear();
        buf3.clear();
        buf4.clear();
        let name1_len = read_line(&mut f_in, &mut buf1)?;
        if name1_len == 0 {
            break;
        }
        let seq1_len = read_line(&mut f_in, &mut buf2)?;
        let name2_len = read_line(&mut f_in, &mut buf3)?;
        let seq2_len = read_line(&mut f_in, &mut buf4)?;

        if name2_len < 5 || seq1_len < 50 || seq1_len != seq2_len || &buf1 != &buf3
                || !check_seq(&buf2[..seq1_len - 1])
                || !check_seq(&buf4[..seq1_len - 1]) {
            panic!("Malformed record in {}:\n    {}    {}    {}    {}```",
                filename,
                String::from_utf8_lossy(&buf1), String::from_utf8_lossy(&buf2),
                String::from_utf8_lossy(&buf3), String::from_utf8_lossy(&buf4));
        }
    }
    Ok(())
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(args.len(), 3);

    if args[2] == "check" {
        check(&args[1]).unwrap()
    } else {
        split_and_interleave(&args[1], &args[2]).unwrap()
    }
}
