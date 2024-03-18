use std::{
    io::{self, Write, BufWriter, BufRead, BufReader},
    fs::File,
};
use flate2::{
    bufread::MultiGzDecoder,
    write::GzEncoder,
    Compression,
};

#[inline]
fn read_line(mut f: impl BufRead, buf: &mut Vec<u8>) -> io::Result<usize> {
    f.read_until(b'\n', buf)
}

fn split_fasta() -> io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(args.len(), 3);

    let in_filename = &args[1];
    let f_in = BufReader::new(File::open(in_filename)?);
    let mut f_in = io::BufReader::new(MultiGzDecoder::new(f_in));

    let out1_filename = args[2].replace("{}", "1");
    let mut f_out1 = BufWriter::new(GzEncoder::new(File::create(out1_filename)?, Compression::default()));
    let out2_filename = args[2].replace("{}", "2");
    let mut f_out2 = BufWriter::new(GzEncoder::new(File::create(out2_filename)?, Compression::default()));

    let mut buf = Vec::with_capacity(1024);
    read_line(&mut f_in, &mut buf)?;
    read_line(&mut f_in, &mut buf)?;
    f_out1.write_all(&buf)?;
    let mut count1 = 1;

    loop {
        buf.clear();
        let n = read_line(&mut f_in, &mut buf)?;
        if buf[n - 3] == b'.' && buf[n - 2] == b'1' {
            break;
        }
        read_line(&mut f_in, &mut buf)?;
        f_out1.write_all(&buf)?;
        count1 += 1;
    }

    std::mem::drop(f_out1);
    read_line(&mut f_in, &mut buf)?;
    f_out2.write_all(&buf)?;
    let mut count2 = 1;

    loop {
        buf.clear();
        if read_line(&mut f_in, &mut buf)? == 0 {
            break;
        }
        read_line(&mut f_in, &mut buf)?;
        f_out2.write_all(&buf)?;
        count2 += 1;
    }
    assert_eq!(count1, count2);
    Ok(())
}

fn main() {
    split_fasta().unwrap()
}
