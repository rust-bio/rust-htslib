use std::fmt;

use serde::de::{self, Deserialize, Deserializer, MapAccess, SeqAccess, Visitor};
use serde::ser::SerializeStruct;
use serde::{Serialize, Serializer};

use bam::record::Record;

impl Serialize for Record {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let core = self.inner().core;
        let mut state = serializer.serialize_struct("Record", 12)?;
        state.serialize_field("tid", &core.tid)?;
        state.serialize_field("pos", &core.pos)?;
        state.serialize_field("bin", &core.bin)?;
        state.serialize_field("mapq", &core.qual)?;
        state.serialize_field("qname_len", &core.l_qname)?;
        state.serialize_field("flag", &core.flag)?;
        state.serialize_field("n_cigar", &core.n_cigar)?;
        state.serialize_field("seq_len", &core.l_qseq)?;
        state.serialize_field("mtid", &core.mtid)?;
        state.serialize_field("mpos", &core.mpos)?;
        state.serialize_field("isize", &core.isize)?;
        state.serialize_field("data", self.data())?;
        state.end()
    }
}

impl<'de> Deserialize<'de> for Record {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        enum Field {
            Tid,
            Pos,
            Bin,
            Mapq,
            QnameLen,
            Flag,
            NCigar,
            SeqLen,
            Mtid,
            Mpos,
            Isize,
            Data,
        };

        impl<'de> Deserialize<'de> for Field {
            fn deserialize<D>(deserializer: D) -> Result<Field, D::Error>
            where
                D: Deserializer<'de>,
            {
                struct FieldVisitor;

                impl<'de> Visitor<'de> for FieldVisitor {
                    type Value = Field;

                    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                        formatter.write_str("expecting a bam field")
                    }

                    fn visit_str<E>(self, value: &str) -> Result<Field, E>
                    where
                        E: de::Error,
                    {
                        match value {
                            "tid" => Ok(Field::Tid),
                            "pos" => Ok(Field::Pos),
                            "bin" => Ok(Field::Bin),
                            "mapq" => Ok(Field::Mapq),
                            "qname_len" => Ok(Field::QnameLen),
                            "flag" => Ok(Field::Flag),
                            "n_cigar" => Ok(Field::NCigar),
                            "seq_len" => Ok(Field::SeqLen),
                            "mtid" => Ok(Field::Mtid),
                            "mpos" => Ok(Field::Mpos),
                            "isize" => Ok(Field::Isize),
                            "data" => Ok(Field::Data),
                            _ => Err(de::Error::unknown_field(value, FIELDS)),
                        }
                    }
                }

                deserializer.deserialize_identifier(FieldVisitor)
            }
        }

        struct RecordVisitor;

        impl<'de> Visitor<'de> for RecordVisitor {
            type Value = Record;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("struct Record")
            }

            fn visit_seq<V>(self, mut seq: V) -> Result<Record, V::Error>
            where
                V: SeqAccess<'de>,
            {
                let tid = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(0, &self))?;
                let pos = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(0, &self))?;
                let bin = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(0, &self))?;
                let mapq = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(0, &self))?;
                let qname_len = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(0, &self))?;
                let flag = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(0, &self))?;
                let n_cigar = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(0, &self))?;
                let seq_len = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(0, &self))?;
                let mtid = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(0, &self))?;
                let mpos = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(0, &self))?;
                let isize = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(0, &self))?;
                let data: Vec<u8> = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(0, &self))?;

                let mut rec = Record::new();
                {
                    let _m = rec.inner_mut();
                    let m = &mut _m.core;
                    m.tid = tid;
                    m.pos = pos;
                    m.bin = bin;
                    m.qual = mapq;
                    m.l_qname = qname_len;
                    m.flag = flag;
                    m.n_cigar = n_cigar;
                    m.l_qseq = seq_len;
                    m.mtid = mtid;
                    m.mpos = mpos;
                    m.isize = isize;
                }

                //println!()
                rec.set_data(&data);
                Ok(rec)
            }

            fn visit_map<V>(self, mut map: V) -> Result<Record, V::Error>
            where
                V: MapAccess<'de>,
            {
                let mut tid = None;
                let mut pos = None;
                let mut bin = None;
                let mut mapq = None;
                let mut qname_len = None;
                let mut flag = None;
                let mut n_cigar = None;
                let mut seq_len = None;
                let mut mtid = None;
                let mut mpos = None;
                let mut isize = None;
                let mut data: Option<Vec<u8>> = None;

                while let Some(key) = map.next_key()? {
                    match key {
                        Field::Tid => {
                            if tid.is_some() {
                                return Err(de::Error::duplicate_field("tid"));
                            }
                            tid = Some(map.next_value()?);
                        }
                        Field::Pos => {
                            if pos.is_some() {
                                return Err(de::Error::duplicate_field("pos"));
                            }
                            pos = Some(map.next_value()?);
                        }
                        Field::Bin => {
                            if bin.is_some() {
                                return Err(de::Error::duplicate_field("bin"));
                            }
                            bin = Some(map.next_value()?);
                        }
                        Field::Mapq => {
                            if mapq.is_some() {
                                return Err(de::Error::duplicate_field("mapq"));
                            }
                            mapq = Some(map.next_value()?);
                        }
                        Field::QnameLen => {
                            if qname_len.is_some() {
                                return Err(de::Error::duplicate_field("qname_len"));
                            }
                            qname_len = Some(map.next_value()?);
                        }
                        Field::Flag => {
                            if flag.is_some() {
                                return Err(de::Error::duplicate_field("flag"));
                            }
                            flag = Some(map.next_value()?);
                        }
                        Field::NCigar => {
                            if n_cigar.is_some() {
                                return Err(de::Error::duplicate_field("n_cigar"));
                            }
                            n_cigar = Some(map.next_value()?);
                        }
                        Field::SeqLen => {
                            if seq_len.is_some() {
                                return Err(de::Error::duplicate_field("seq_len"));
                            }
                            seq_len = Some(map.next_value()?);
                        }
                        Field::Mtid => {
                            if mtid.is_some() {
                                return Err(de::Error::duplicate_field("mtid"));
                            }
                            mtid = Some(map.next_value()?);
                        }
                        Field::Mpos => {
                            if mpos.is_some() {
                                return Err(de::Error::duplicate_field("mpos"));
                            }
                            mpos = Some(map.next_value()?);
                        }
                        Field::Isize => {
                            if isize.is_some() {
                                return Err(de::Error::duplicate_field("isize"));
                            }
                            isize = Some(map.next_value()?);
                        }
                        Field::Data => {
                            if data.is_some() {
                                return Err(de::Error::duplicate_field("data"));
                            }
                            data = Some(map.next_value()?);
                        }
                    }
                }

                let tid = tid.ok_or_else(|| de::Error::missing_field("tid"))?;
                let pos = pos.ok_or_else(|| de::Error::missing_field("pos"))?;
                let bin = bin.ok_or_else(|| de::Error::missing_field("bin"))?;
                let mapq = mapq.ok_or_else(|| de::Error::missing_field("mapq"))?;
                let qname_len = qname_len.ok_or_else(|| de::Error::missing_field("qname_len"))?;
                let flag = flag.ok_or_else(|| de::Error::missing_field("flag"))?;
                let n_cigar = n_cigar.ok_or_else(|| de::Error::missing_field("n_cigar"))?;
                let seq_len = seq_len.ok_or_else(|| de::Error::missing_field("seq_len"))?;
                let mtid = mtid.ok_or_else(|| de::Error::missing_field("mtid"))?;
                let mpos = mpos.ok_or_else(|| de::Error::missing_field("mpos"))?;
                let isize = isize.ok_or_else(|| de::Error::missing_field("isize"))?;
                let data = data.ok_or_else(|| de::Error::missing_field("data"))?;

                let mut rec = Record::new();
                {
                    let _m = rec.inner_mut();
                    let m = &mut _m.core;
                    m.tid = tid;
                    m.pos = pos;
                    m.bin = bin;
                    m.qual = mapq;
                    m.l_qname = qname_len;
                    m.flag = flag;
                    m.n_cigar = n_cigar;
                    m.l_qseq = seq_len;
                    m.mtid = mtid;
                    m.mpos = mpos;
                    m.isize = isize;
                }

                rec.set_data(&data);
                Ok(rec)
            }
        }

        const FIELDS: &'static [&'static str] = &[
            "tid", "pos", "bin", "qual", "l_qname", "flag", "n_cigar", "seq_len", "mtid", "mpos",
            "isize", "data",
        ];
        deserializer.deserialize_struct("Record", FIELDS, RecordVisitor)
    }
}

#[cfg(test)]
mod tests {
    use bam::record::Record;
    use bam::Read;
    use bam::Reader;

    use std::path::Path;

    use bincode::{deserialize, serialize, Infinite};
    use serde_json;

    #[test]
    fn test_bincode() {
        let mut bam = Reader::from_path(&Path::new("test/test.bam"))
            .ok()
            .expect("Error opening file.");

        let mut recs = Vec::new();
        for record in bam.records() {
            recs.push(record.unwrap());
        }

        let encoded: Vec<u8> = serialize(&recs, Infinite).unwrap();
        let decoded: Vec<Record> = deserialize(&encoded[..]).unwrap();
        assert_eq!(recs, decoded);
    }

    #[test]
    fn test_serde_json() {
        let mut bam = Reader::from_path(&Path::new("test/test.bam"))
            .ok()
            .expect("Error opening file.");

        let mut recs = Vec::new();
        for record in bam.records() {
            recs.push(record.unwrap());
        }

        let encoded: String = serde_json::to_string(&recs).unwrap();
        println!("encoded: {}", encoded);
        let decoded: Vec<Record> = serde_json::from_str(&encoded).unwrap();
        assert_eq!(recs, decoded);
    }

}
