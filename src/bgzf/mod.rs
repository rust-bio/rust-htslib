// Copyright 2020 Manuel Landesfeind, Evotec International GmbH
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//!
//! Module for working with bgzipped file.
//!

use std::ffi;
use std::path::Path;
use url::Url;

use crate::htslib;
use crate::tpool::ThreadPool;

pub mod errors;
pub use errors::{Error, Result};

fn path_as_bytes<'a, P: 'a + AsRef<Path>>(path: P, must_exist: bool) -> Result<Vec<u8>> {
    if path.as_ref().exists() || !must_exist {
        Ok(path
            .as_ref()
            .to_str()
            .ok_or(Error::NonUnicodePath)?
            .as_bytes()
            .to_owned())
    } else {
        Err(Error::FileNotFound {
            path: path.as_ref().to_owned(),
        })
    }
}

/// Test if a file is a Bgzip compressed file
///
/// # Arguments
///
/// * `path` - the path to test.
///
/// # Returns:
/// Will return `Ok(true)` or `Ok(false)` if the file at `path` is BGZIP compressed. Will return an `Err` in
/// cases where no testing is possible.
pub fn is_bgzip<P: AsRef<Path>>(path: P) -> Result<bool, Error> {
    let byte_path = path_as_bytes(path, true)?;
    let cpath = ffi::CString::new(byte_path).unwrap();
    let is_bgzf = unsafe { htslib::bgzf_is_bgzf(cpath.as_ptr()) == 1 };
    Ok(is_bgzf)
}

/// A reader that transparently reads uncompressed, gzip, and bgzip files.
#[derive(Debug)]
pub struct Reader {
    inner: *mut htslib::BGZF,
}

impl Reader {
    /// Create a new Reader to read from stdin.
    pub fn from_stdin() -> Result<Self, Error> {
        Self::new(b"-")
    }

    /// Create a new Reader from a path.
    ///
    /// # Arguments
    ///
    /// * `path` - the path to open.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, Error> {
        Self::new(&path_as_bytes(path, true)?)
    }

    /// Create a new Reader from an URL.
    ///
    /// # Arguments
    ///
    /// * `url` - the url to open
    pub fn from_url(url: &Url) -> Result<Self, Error> {
        Self::new(url.as_str().as_bytes())
    }

    /// Internal function to create a Reader from some sort of path (could be file path but also URL).
    /// The path or URL will be handled by the c-implementation transparently.
    ///
    /// # Arguments
    ///
    /// * `path` - the path or URL to open
    fn new(path: &[u8]) -> Result<Self, Error> {
        let mode = ffi::CString::new("r").unwrap();
        let cpath = ffi::CString::new(path).unwrap();
        let inner = unsafe { htslib::bgzf_open(cpath.as_ptr(), mode.as_ptr()) };
        Ok(Self { inner })
    }

    /// Set the thread pool to use for parallel decompression. 
    ///
    /// # Arguments
    ///
    /// * `tpool` - the thread-pool to use
		pub fn set_thread_pool(&mut self, tpool: &ThreadPool) -> Result<()> {
    	let b = tpool.handle.borrow_mut();
			let r = unsafe {
				htslib::bgzf_thread_pool(
					self.inner, 
					b.inner.pool as *mut _, 
					0) // let htslib decide on the queue-size
			};
    	
    	if r != 0 {
    	    Err(Error::ThreadPool)
    	} else {
    	    Ok(())
    	}
	}
}

impl std::io::Read for Reader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let nbytes = unsafe {
            htslib::bgzf_read(
                self.inner,
                buf.as_mut_ptr() as *mut libc::c_void,
                buf.len() as u64,
            )
        };
        if nbytes < 0 {
            Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "Can not read",
            ))
        } else {
            Ok(nbytes as usize)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    // Define paths to the test files
    const FN_PLAIN: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/test/bgzip/plain.vcf");
    const FN_GZIP: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/test/bgzip/gzip.vcf.gz");
    const FN_BGZIP: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/test/bgzip/bgzip.vcf.gz");

    const CONTENT: &str = include_str!("../../test/bgzip/plain.vcf");

    #[test]
    fn test_is_bgzip_plain() {
        assert!(
            !is_bgzip(FN_PLAIN).unwrap(),
            "Plain file not detected as BGZIP"
        );
        assert!(
            !is_bgzip(FN_GZIP).unwrap(),
            "Zip file not detected as BGZIP"
        );
        assert!(is_bgzip(FN_BGZIP).unwrap(), "Bgzip file detected as BGZIP");
    }

    #[test]
    fn test_open_plain() {
        let r_result = Reader::from_path(FN_PLAIN);
        assert!(r_result.is_ok(), "Open plain file with Bgzip reader");

        let mut my_content = String::new();
        let reading_result = r_result.unwrap().read_to_string(&mut my_content);
        assert!(
            reading_result.is_ok(),
            "Reading plain file into buffer is ok"
        );
        assert_eq!(
            reading_result.unwrap(),
            190,
            "Reading plain file into buffer is correct size"
        );
        assert_eq!(
            my_content, CONTENT,
            "Reading plain file with correct content"
        );
    }

    #[test]
    fn test_open_gzip() {
        let r_result = Reader::from_path(FN_GZIP);
        assert!(r_result.is_ok(), "Open gzip file with Bgzip reader");

        let mut my_content = String::new();
        let reading_result = r_result.unwrap().read_to_string(&mut my_content);
        assert!(
            reading_result.is_ok(),
            "Reading gzip file into buffer is ok"
        );
        assert_eq!(
            reading_result.unwrap(),
            190,
            "Reading gzip file into buffer is correct size"
        );
        assert_eq!(
            my_content, CONTENT,
            "Reading gzip file with correct content"
        );
    }

    #[test]
    fn test_open_bgzip() {
        let r_result = Reader::from_path(FN_BGZIP);
        assert!(r_result.is_ok(), "Open bgzip file with Bgzip reader");

        let mut my_content = String::new();
        let reading_result = r_result.unwrap().read_to_string(&mut my_content);
        assert!(
            reading_result.is_ok(),
            "Reading bgzip file into buffer is ok"
        );
        assert_eq!(
            reading_result.unwrap(),
            190,
            "Reading bgzip file into buffer is correct size"
        );
        assert_eq!(
            my_content, CONTENT,
            "Reading bgzip file with correct content"
        );
    }
    #[test]
    fn test_set_threadpool() {
        let r_result = Reader::from_path(FN_BGZIP);
        assert!(r_result.is_ok(), "Open bgzip file with Bgzip reader");
				let mut r = r_result.unwrap();

				let tpool_result = ThreadPool::new(5);
				assert!(tpool_result.is_ok(), "Creating thread pool");
				let tpool = tpool_result.unwrap();

				let set_result = r.set_thread_pool(&tpool);
				assert_eq!(set_result, Ok(()), "Setting thread pool okay");

        let mut my_content = String::new();
        let reading_result = r.read_to_string(&mut my_content);
        assert!(
            reading_result.is_ok(),
            "Reading bgzip file into buffer is ok - using a threadpool"
        );
        assert_eq!(
            reading_result.unwrap(),
            190,
            "Reading bgzip file into buffer is correct size using a threadpool"
        );
        assert_eq!(
            my_content, CONTENT,
            "Reading bgzip file with correct content using a threadpool"
        );
    }
}
