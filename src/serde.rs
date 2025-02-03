#![cfg(feature="serde")]
use super::Histogram;
use super::bin::{Bin, NAN_BIN};

use std::marker::PhantomData;
use std::fmt;

use serde::{de, de::{Deserializer, MapAccess, Visitor}, ser::{Serialize, SerializeMap, Serializer}, Deserialize};

impl Serialize for Bin {
    fn serialize<S>(&self, serializer: S) -> Result<<S as Serializer>::Ok, <S as Serializer>::Error>
    where S: Serializer {
        serializer.serialize_str(self.to_string().as_str())
    }
}

struct BinVisitor;

impl<'de> Visitor<'de> for BinVisitor {
    type Value = Bin;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("an OpenHistogram Bin")
    }

    fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
    where E: de::Error,
    {
        if value == "NaN" {
            Ok(super::bin::NAN_BIN)
        }
        else {
            match value.parse::<f64>() {
                Ok(v) => Ok(v.into()),
                Err(_) => Err(E::custom(format!("bin {} invalid", value)))
            }
        }
    }
}

impl<'de> Deserialize<'de> for Bin {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_str(BinVisitor)
    }
}

impl Serialize for Histogram {
    fn serialize<S>(&self, serializer: S) -> Result<<S as Serializer>::Ok, <S as Serializer>::Error>
    where S: Serializer {
        let iter = self.bvs.iter().filter(|x| x.1 != &0u64);
        let mut map = serializer.serialize_map(Some(iter.clone().count()))?;
        for (bin, count) in iter {
            map.serialize_entry(bin, count)?;
        }
        if self.nan_count != 0 {
            map.serialize_entry(&NAN_BIN, &self.nan_count)?;
        }
        map.end()
    }
}

struct HistogramVisitor {
    marker: PhantomData<fn() -> Histogram>
}

impl HistogramVisitor {
    fn new() -> Self {
        HistogramVisitor {
            marker: PhantomData
        }
    }
}

impl<'de> Visitor<'de> for HistogramVisitor {
    type Value = Histogram;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("an OpenHistogram")
    }

    fn visit_map<M>(self, mut access: M) -> Result<Self::Value, M::Error>
    where
        M: MapAccess<'de>,
    {
        let mut hist = Histogram::new();
        while let Some((key, value)) = access.next_entry()? {
            hist.insert(key, value);
        }
        Ok(hist)
    }
}

impl<'de> Deserialize<'de> for Histogram {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_map(HistogramVisitor::new())
    }
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn histogram_serde_json() -> Result<(), super::super::Error> {
        use serde_json::json;
        type Error = super::super::Error;

        let mut h1 = Histogram::new();
        for v in vec![ 0.123, 0.0, 0.43, 0.41, 0.415, 0.2201, 3201.0, 0.125, 0.13, 345e234, -100e200 ] {
            h1.insert(v.into(), 1);
        }
        let json_data =  json!(h1).to_string();
        match serde_json::from_str(json_data.as_str()) {
            Ok(h2) => assert_eq!(h1, h2),
            Err(e) => {
                return Err(Error::CorruptEncoding(e.to_string()))
            }
        }
        Ok(())
    }
}