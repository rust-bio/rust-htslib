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
///
///```
/// // Index a sorted bcf file with csi.
/// let bcf_path = concat!(env!("CARGO_MANIFEST_DIR"), "/test/test_multi.bcf");
/// rust_htslib::bcf::index::build(&bcf_path, Some(&"built_test.csi"), /*n_threads=*/ 4u32, rust_htslib::bcf::index::Type::Csi(14)).expect("Failed to build csi index for bcf file.");
/// assert!(std::path::Path::new(&"built_test.csi").exists());
///
/// // Index a bgzip-compresed vcf file with tabix.
/// let vcf_path = concat!(env!("CARGO_MANIFEST_DIR"), "/test/test_left.vcf.gz");
/// rust_htslib::bcf::index::build(&vcf_path, Some(&"built_test_vcf.tbx"), /*n_threads=*/ 4u32, rust_htslib::bcf::index::Type::Tbx).expect("Failed to build tbx index for vcf file.");
/// assert!(std::path::Path::new(&"built_test_vcf.tbx").exists());
///
/// // Cannot build a tbi index for a bcf file: returns an Err(BcfBuildError).
/// assert!(std::panic::catch_unwind(|| rust_htslib::bcf::index::build(bcf_path, Some("built_test.tbi"), 4u32, rust_htslib::bcf::index::Type::Tbx).unwrap()).is_err());
///
/// // Cannot built a csi index for a vcf file: returns an Err(BcfBuildError).
/// let vcf_path = concat!(env!("CARGO_MANIFEST_DIR"), "/test/test_various.vcf");
/// assert!(std::panic::catch_unwind(|| rust_htslib::bcf::index::build(&vcf_path, Some(&"built_test_vcf.csi"), /*n_threads=*/ 4u32, rust_htslib::bcf::index::Type::Csi(14)).expect("Failed to build csi index for vcf file.")).is_err());
///```
///
pub fn build<P: AsRef<std::path::Path>>(
    bcf_path: P,
    idx_path: Option<P>,
    n_threads: u32,
    index_type: Type,
) -> Result<(), BcfBuildError> {
    let min_shift = index_type.min_shift();
    let idx_path_cstr = idx_path.and_then(|x| utils::path_to_cstring(&x));
    let bcf_path = utils::path_to_cstring(&bcf_path).ok_or(BcfBuildError {
        msg: format!(
            "Failed to format bcf_path to cstring: {:?}",
            bcf_path.as_ref().display()
        ),
    })?;
    let return_code = unsafe {
        htslib::bcf_index_build3(
            bcf_path.as_ptr(),
            idx_path_cstr
                .as_ref()
                .map_or(std::ptr::null(), |p| p.as_ptr()),
            min_shift,
            n_threads as i32,
        )
    };
    if return_code == 0 {
        Ok(())
    } else {
        Err(BcfBuildError {
            msg: format!(
                "Failed to build  bcf index. Error: {return_code:?}/{}",
                BcfBuildError::error_message(return_code)
            ),
        })
    }
}
