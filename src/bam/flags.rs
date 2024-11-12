//! A module that provides an easier way to work with SAM flags.
//! It achieves this by providing a struct (`Flag`) with associated constants representing
//! flags, and the `check_flag()` function that allows testing if specific flags are set or not set.
//!
//! For example, the following code tests if the read flag has the FIRST_IN_PAIR flag set and the MATE_UNMAPPED flag not set:
//! ```
//! use rust_htslib::bam::{Flag, check_flag};
//! # let read_flag = record.flag(); in general this is the way to obtian a flag.
//! let read_flag = 64;
//! assert_eq!(check_flag(read_flag, Flag::FIRST_IN_PAIR, Flag::MATE_UNMAPPED), true);
//! ```

///
/// This structure contains constants representing SAM flag values as u16.
/// Using this structure incurs no runtime cost.
///
/// ```
/// use rust_htslib::bam::{Flag};
/// to get the value of a flag representing a read paired, and reversly mapped.
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

pub fn check_flag(flag: u16, in_: u16, not_in: u16) -> bool {
    //! This function uses bitwise operations to test if flags are set or not.
    //!
    //! # Arguments
    //!
    //! * `flag`: u16 - The record flag you want to test
    //! * `in_`: u16 - The flags you want to check if they are set (use 0 for no test)
    //! * `not_in`: u16 - The flags you want to check if they are not set (use 0 for no test)
    //!
    //! # Usage:
    //! example: let test if a flag is both paired and fisrt in pair
    //! ```
    //! use rust_htslib::bam::{Flag, check_flag};
    //! let read_flag = 18
    //! assert_eq!(check_flag(read_flag, Flag::PAIRED +Flag::FIRST_IN_PAIR, 0), true);
    //! ```
    //! let test that the read is mapped. READ_UNMAPPED
    //! ```
    //! use rust_htslib::bam::{Flag, check_flag};
    //! let read_flag = 18
    //! assert_eq!(check_flag(read_flag, 0, Flag::READ_UNMAPPED), true);
    //! ```
    //!
    //! Finally let do a more complexe real example test:
    //! ```
    //! use rust_htslib::bam::{Flag, check_flag};
    //! let read_flag = 19
    //! assert_eq!(check_flag(read_flag, Flag::PAIRED + Flag::PROPERLY_PAIRED + Flag::READ_RERVERSE , Flag::READ_UNMAPPED + Flag::MATE_UNMAPPED), true);
    //! ```
    //!
    //binary flag check
    //assert that: - in_ is in n
    //             - not_in is not in n
    // bitwise operation
    if (not_in & flag) != 0 {
        return false;
    }
    if (in_ & flag) != in_ {
        return false;
    }
    return true;
}
