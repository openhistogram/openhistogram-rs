use std::io::{Read,Write};
use byteorder::{WriteBytesExt,ReadBytesExt,NetworkEndian};
use super::Histogram;
use super::bin::{Bin, NAN_BIN};

#[cfg(feature="base64")]
use base64::{engine::general_purpose::STANDARD_NO_PAD, Engine as _};

fn oh_serialize_u64<W: Write>(buffer: &mut W, v: u64) -> Result<(),super::Error> {
    let lz = (v.leading_zeros() / 8) as u8;
    assert!(lz <= 8);
    let level = if lz == 8 { 0u8 } else { 7u8 - lz };
    buffer.write_u8(level)?;
    for i in lz..0 {
        buffer.write_u8(((v >> (i*8)) & 0xff) as u8)?;
    }
    buffer.write_u8((v & 0xff) as u8)?;
    Ok(())
}
fn oh_deserialize_u64<R: Read>(buffer: &mut R) -> Result<u64,super::Error> {
    let mut v  = 0u64;
    let level = buffer.read_u8()?;
    if level > 7 {
        Err(super::Error::CorruptEncoding(format!("illegal varbit {}", level).to_string()))
    } else {
        for _ in 0..(level+1) {
            v = (v << 8) | (buffer.read_u8()? as u64);
        }
        Ok(v)
    }
}

impl Bin {
    /// This implements the native OpenHistogram raw serialization.
    pub fn oh_serialize<W: Write>(&self, buffer: &mut W) -> Result<(), super::Error> {
        buffer.write_i8(self.val)?;
        buffer.write_i8(self.exp)?;
        Ok(())
    }
    /// This implements the native OpenHistogram raw deserialization.
    pub fn oh_deserialize<R: Read>(buffer: &mut R) -> Result<Bin, super::Error> {
        let bin = Bin { val: buffer.read_i8()?, exp: buffer.read_i8()? };
        if bin.is_valid() {
            Ok(bin)
        }
        else {
            Err(super::Error::CorruptEncoding(format!("invalid {:?}", bin).to_string()))
        }
    }
}

impl Histogram {
    /// This implements the native OpenHistogram raw serialization.
    pub fn oh_serialize<W: Write>(&self, buffer: &mut W) -> Result<(), super::Error> {
        let bins: Vec<(&Bin, &u64)> = self.bvs.iter().filter(|x| x.1 > &0u64).collect();
        let len: u16 = match self.nan_count {
            0 => bins.len() as u16,
            _ => bins.len() as u16 + 1,
        };
        buffer.write_u16::<NetworkEndian>(len)?;
        for (bin, count) in bins {
            bin.oh_serialize(buffer)?;
            oh_serialize_u64(buffer, *count)?;
        }
        if self.nan_count != 0 {
            NAN_BIN.oh_serialize(buffer)?;
            oh_serialize_u64(buffer, self.nan_count)?;
        }
        Ok(())
    }
    /// This implements the native OpenHistogram raw serialization.
    pub fn oh_deserialize<R: Read>(buffer: &mut R) -> Result<Histogram, super::Error> {
        let mut h = Histogram::new();
        let len = buffer.read_u16::<NetworkEndian>()?;
        for _ in 0..len {
            let bin = Bin::oh_deserialize(buffer)?;
            let count = oh_deserialize_u64(buffer)?;
            h.insert(bin, count);
        }
        Ok(h)
    }

    #[cfg(feature="base64")]
    /// This implements the native OpenHistogram base64 serialization.
    pub fn to_base64(&self) -> Result<String, super::Error> {
        let mut buffer: Vec<u8> = vec![];
        self.oh_serialize(&mut buffer)?;
        Ok(STANDARD_NO_PAD.encode(buffer))
    }
    #[cfg(feature="base64")]
    /// This implements the native OpenHistogram base64 deserialization.
    pub fn from_base64(b64: String) -> Result<Histogram, super::Error> {
        match STANDARD_NO_PAD.decode(b64) {
            Err(e) => Err(super::Error::CorruptEncoding(format!("{}", e).to_string())),
            Ok(data) => {
                let mut buffer = std::io::Cursor::new(data);
                Histogram::oh_deserialize(&mut buffer)
            }
        }
    }
}


#[cfg(test)]
mod tests {
    use std::io;

    use super::*;

    #[test]
    fn raw_oh_repr() -> Result<(), crate::Error> {
        let mut href = Histogram::new();
        for v in vec![ 0.123, 0.0, 0.43, 0.41, 0.415, 0.2201, 0.3201, 0.125, 0.13 ] {
            href.insert(v.into(), 1);
        }
        let mut href_buf: Vec<u8> = vec![];
        href.oh_serialize(&mut href_buf)?;
        let mut c = io::Cursor::new(href_buf);
        let h =  Histogram::oh_deserialize(&mut c)?;
        assert_eq!(h, href);
        Ok(())
    }
    #[test]
    #[cfg(feature="base64")]
    fn base64_repr() -> Result<(), crate::Error> {
        let mut href = Histogram::new();
        let b64e = String::from("AAcAAAABDP8AAg3/AAEW/wABIP8AASn/AAIr/wAB");
        for v in vec![ 0.123, 0.0, 0.43, 0.41, 0.415, 0.2201, 0.3201, 0.125, 0.13 ] {
            href.insert(v.into(), 1);
        }
        let b64ref = href.to_base64()?;
        assert_eq!(b64ref, b64e);
        let h =  Histogram::from_base64(b64e)?;
        assert_eq!(h, href);
        Ok(())
    }
}