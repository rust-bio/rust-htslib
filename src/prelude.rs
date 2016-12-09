// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! The purpose of this module is to provide reexports of core traits so that they can be then glob-imported all at once:
//!
//! ```
//! use rust_htslib::prelude::*;
//! ```

pub use bam::Read;
pub use bcf::record::Numeric;
