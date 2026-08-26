//! A module that provides an easier way to work with SAM flags.
//! It achieves this by providing a struct (`Flag`) with associated constants representing
//! flags. And the `is_in_flag()` and `is_not_in_flag()` functions that allow testing if specific flags are set or not set.
//!
//! ```
//! use rust_htslib::bam::flags{Flag, is_in_flag, is_not_in_flag};
//! let read_flag = record.flag(); // in general this is the way to obtain a flag.
//! let read_flag = 64;
//! assert_eq!(is_in_flag(read_flag, Flag::FIRST_IN_PAIR), true);
//! assert_eq!(is_not_in_flag(read_flag, Flag::MATE_UNMAPPED), true);
//! ```
///
/// This structure contains constants representing SAM flag values as u16.
/// Using this structure incurs no runtime cost.
///
/// ```
/// use rust_htslib::bam::flags::Flag;
/// // to get the value of a flag representing a read paired, and reversly mapped.
/// let flag = Flag::PAIRED + Flag::READ_RERVERSE;
///
/// ```
pub struct Flag;

impl Flag {
    pub const PAIRED: u16 = 1;
    pub const PROPERLY_PAIRED: u16 = 2;
    pub const READ_UNMAPPED: u16 = 4;
    pub const MATE_UNMAPPED: u16 = 8;
    pub const READ_RERVERSE: u16 = 16;
    pub const MATE_REVERSE: u16 = 32;
    pub const FIRST_IN_PAIR: u16 = 64;
    pub const SECOND_IN_PAIR: u16 = 128;
    pub const NOT_PRIMARY_ALN: u16 = 256;
    pub const FAIL_QC: u16 = 512;
    pub const DUPLICATE: u16 = 1024;
    pub const SUPPLEMENTARY: u16 = 2048;
}

pub fn is_not_in_flag(flag: u16, not_in: u16) -> bool {
    //! This function uses bitwise operations to test if flags are not set
    //! # Arguments
    //! * `flag`: u16 - The record flag you want to test
    //! * `not_in`: u16 - The flags you want to check if they are not set (use 0 for no test)
    //!
    //! # Usage:
    //! example: let test if a flag is primary alignment and did not fail QC
    //! ```
    //! use rust_htslib::bam::flags;
    //! use rust_htslib::bam::flags::Flag;
    //! let read_flag = 65;
    //! assert_eq!(flags::is_not_in_flag(read_flag, Flag::NOT_PRIMARY_ALN + Flag::FAIL_QC), true);
    //! ```
    //! let test that the read is mapped.
    //! ```
    //!
    //! use rust_htslib::bam::flags::{Flag, is_not_in_flag};
    //! let read_flag = 18;
    //! assert_eq!(is_not_in_flag(read_flag, Flag::READ_UNMAPPED), true);
    //! ```
    //!
    if (not_in & flag) != 0 {
        return false;
    }
    true
}
pub fn is_in_flag(flag: u16, in_: u16) -> bool {
    //! This function uses bitwise operations to test if flags are set
    //! # Arguments
    //! * `flag`: u16 - The record flag you want to test
    //! * `in_`: u16 - The flags you want to check if they are set (use 0 for no test)
    //!
    //! # Usage:
    //! example: let test if a flag is both paired and first in pair
    //! ```
    //! use rust_htslib::bam::flags::{Flag, is_in_flag};
    //! let read_flag = 65;
    //! assert_eq!(is_in_flag(read_flag, Flag::PAIRED + Flag::FIRST_IN_PAIR), true);
    //! ```

    if (in_ & flag) != in_ {
        return false;
    }
    true
}
