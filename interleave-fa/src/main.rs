use std::{
    io::{self, Write, BufWriter, Read, BufRead, BufReader},
    fs::File,
};
use flate2::{
    bufread::MultiGzDecoder,
    write::GzEncoder,
    Compression,
};

#[inline(always)]
fn read_line(mut f: impl BufRead, buf: &mut Vec<u8>) -> io::Result<usize> {
    f.read_until(b'\n', buf)
}

fn split_and_interleave() -> io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(args.len(), 3);

    let in_filename = &args[1];
    let f_in = BufReader::new(File::open(in_filename)?);
    let mut f_in = BufReader::new(MultiGzDecoder::new(f_in));
    let out_filename = &args[2];
    let mut f_out = BufWriter::with_capacity(131_072,
        GzEncoder::new(File::create(out_filename)?, Compression::default()));
    let tmp_filename = format!("{}.tmp.fasta", out_filename);
    let mut f_tmp = BufWriter::with_capacity(131_072, File::create(&tmp_filename)?);

    let mut buf = Vec::with_capacity(1024);
    read_line(&mut f_in, &mut buf)?;
    read_line(&mut f_in, &mut buf)?;
    f_tmp.write_all(&buf)?;
    let mut count1 = 1;
    loop {
        buf.clear();
        let n = read_line(&mut f_in, &mut buf)?;
        if buf[n - 3] == b'.' && buf[n - 2] == b'1' {
            break;
        }

        read_line(&mut f_in, &mut buf)?;
        f_tmp.write_all(&buf)?;
        count1 += 1;
    }
    f_tmp.flush()?;
    std::mem::drop(f_tmp);

    let old_buf = buf;
    let mut buf = Vec::with_capacity(1024);
    let mut f_tmp = BufReader::new(File::open(&tmp_filename)?);
    read_line(&mut f_tmp, &mut buf)?;
    read_line(&mut f_tmp, &mut buf)?;
    buf.extend_from_slice(&old_buf);
    read_line(&mut f_in, &mut buf)?;
    f_out.write_all(&buf)?;

    for _ in 1..count1 {
        buf.clear();
        read_line(&mut f_tmp, &mut buf)?;
        read_line(&mut f_tmp, &mut buf)?;
        assert!(read_line(&mut f_in, &mut buf)? > 0);
        read_line(&mut f_in, &mut buf)?;
        f_out.write_all(&buf)?;
    }
    std::mem::drop(f_tmp);
    std::mem::drop(f_out);

    std::fs::remove_file(&tmp_filename)?;

    let last_recs = String::from_utf8_lossy(&buf);
    let last_recs: Vec<_> = last_recs.split("\n").collect();
    assert_eq!(last_recs[0], last_recs[2]);
    buf.clear();
    buf.resize(16, 0);
    assert_eq!(f_in.read(&mut buf)?, 0);
    Ok(())
}

fn main() {
    split_and_interleave().unwrap()
}