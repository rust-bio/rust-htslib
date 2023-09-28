use crate::{htslib, utils};

#[derive(Debug)]
pub struct BcfBuildError {
    pub msg: String,
}

/// Index type to build.
pub enum Type {
    /// Tabix index
    Tbx,
    /// CSI index, with given minimum shift
    Csi(u32),
}

impl Type {
    fn min_shift(&self) -> i32 {
        match self {
            Self::Tbx => 0,
            Self::Csi(x) => *x as i32,
        }
    }
}

impl BcfBuildError {
    pub fn error_message(error: i32) -> &'static str {
        match error {
            -1 => "indexing failed",
            -2 => "opening @fn failed",
            -3 => "format not indexable",
            -4 => "failed to create and/or save the index",
            _ => "unknown error",
        }
    }
}
impl std::error::Error for BcfBuildError {}

impl std::fmt::Display for BcfBuildError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "BcfBuildError{{msg: {}}}", self.msg)
    }
}

/// Build a bcf or vcf.gz index.
/// Builds tbi or csi depending on index_type.
pub fn build<P: AsRef<std::path::Path>>(
    bcf_path: P,
    idx_path: Option<P>,
    n_threads: u32,
    index_type: Type,
) -> Result<(), BcfBuildError> {
    let min_shift = index_type.min_shift();
    let idx_path_cstr = idx_path.map(|p| utils::path_to_cstring(&p).expect("path_to_cstring"));
    let ret = unsafe {
        htslib::bcf_index_build3(
            utils::path_to_cstring(&bcf_path).unwrap().as_ptr(),
            idx_path_cstr.map_or(std::ptr::null(), |p| p.as_ptr()),
            min_shift,
            n_threads as i32,
        )
    };
    match ret {
        0 => Ok(()),
        e => Err(BcfBuildError {
            msg: format!(
                "Failed to build  bcf index. Error: {e:?}/{}",
                BcfBuildError::error_message(e)
            ),
        }),
    }
}
