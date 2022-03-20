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

use crate::errors::{Error, Result};

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
            htslib::bgzf_thread_pool(self.inner, b.inner.pool as *mut _, 0) // let htslib decide on the queue-size
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

/// The CompressionLevel used by the underlying GZIP writer
/// Note that the special level NoCompression will not use the GZIP writer.
/// Compression levels in BGZF files
///
/// * Uncompressed: No compression, zlib level 0
/// * Fastest: Lowest compression level, zlib level 1
/// * Maximum: Highest compression level, zlib level 9
/// * Default: Default compression level, zlib level 6
/// * Level(i): Custom compression level in the range [0, 9]
/// * NoCompression: No compression, zlib not used. Output will be identical to input
#[derive(Debug, Clone, Copy)]
pub enum CompressionLevel {
    Default,
    NoCompression,
    Uncompressed,
    Fastest,
    Maximum,
    Level(i8),
}
impl CompressionLevel {
    // Convert and check the variants of the `CompressionLevel` enum to a numeric level
    fn convert(self) -> Result<i8> {
        match self {
            CompressionLevel::NoCompression => Ok(-2),
            CompressionLevel::Default => Ok(-1),
            CompressionLevel::Uncompressed => Ok(0),
            CompressionLevel::Fastest => Ok(1),
            CompressionLevel::Maximum => Ok(9),
            CompressionLevel::Level(i @ -2..=9) => Ok(i),
            CompressionLevel::Level(i) => Err(Error::BgzfInvalidCompressionLevel { level: i }),
        }
    }
}

/// A writer that writes uncompressed, gzip, and bgzip files.
#[derive(Debug)]
pub struct Writer {
    inner: *mut htslib::BGZF,
    tpool: Option<ThreadPool>,
}

impl Writer {
    /// Create a new Writer to write to stdout with default compression.
    pub fn from_stdout() -> Result<Self, Error> {
        Self::from_stdout_with_compression(CompressionLevel::Default)
    }

    /// Create a new Writer to write to stdout with specific compression
    ///
    /// # Arguments
    ///
    /// * `level` the compression level to use
    pub fn from_stdout_with_compression(level: CompressionLevel) -> Result<Self, Error> {
        Self::new(b"-", level)
    }

    /// Create a new Writer from a path with default compression.
    ///
    /// # Arguments
    ///
    /// * `path` - the path to open.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, Error> {
        Self::from_path_with_level(path, CompressionLevel::Default)
    }

    /// Create a new Writer from a path with a specific compression level.
    ///
    /// # Arguments
    ///
    /// * `path` - the path to open.
    pub fn from_path_with_level<P: AsRef<Path>>(
        path: P,
        level: CompressionLevel,
    ) -> Result<Self, Error> {
        Self::new(&path_as_bytes(path, false)?, level)
    }

    /// Internal function to create a Writer from a file path
    ///
    /// # Arguments
    ///
    /// * `path` - the path or URL to open
    fn new(path: &[u8], level: CompressionLevel) -> Result<Self, Error> {
        let mode = Self::get_open_mode(level)?;
        let cpath = ffi::CString::new(path).unwrap();
        let inner = unsafe { htslib::bgzf_open(cpath.as_ptr(), mode.as_ptr()) };
        Ok(Self { inner, tpool: None })
    }

    /// Internal function to convert compression level to "mode"
    /// bgzf.c expects mode for writers to be one of: 'w', 'wu', 'w#', where # is 0-9.
    /// # Arguments
    ///
    /// * `level` - the level of compression to use
    fn get_open_mode(level: CompressionLevel) -> Result<ffi::CString, Error> {
        let write_string = match level.convert() {
            Ok(-2) => "wu".to_string(),
            Ok(-1) => "w".to_string(),
            Ok(n @ 0..=9) => format!("w{}", n),
            Err(e) => return Err(e),
            // This should be unreachable
            Ok(i) => return Err(Error::BgzfInvalidCompressionLevel { level: i }),
        };
        return Ok(ffi::CString::new(write_string).unwrap());
    }

    /// Set the thread pool to use for parallel compression.
    ///
    /// # Arguments
    ///
    /// * `tpool` - the thread-pool to use
    pub fn set_thread_pool(&mut self, tpool: &ThreadPool) -> Result<()> {
        self.tpool = Some(tpool.clone());
        let b = tpool.handle.borrow_mut();
        let r = unsafe {
            htslib::bgzf_thread_pool(self.inner, b.inner.pool as *mut _, 0) // let htslib decide on the queue-size
        };

        if r != 0 {
            Err(Error::ThreadPool)
        } else {
            Ok(())
        }
    }
}

impl std::io::Write for Writer {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let nbytes = unsafe {
            htslib::bgzf_write(
                self.inner,
                buf.as_ptr() as *mut libc::c_void,
                buf.len() as u64,
            )
        };
        if nbytes < 0 {
            Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "Can not write",
            ))
        } else {
            Ok(nbytes as usize)
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        let exit_code: i32 = unsafe { htslib::bgzf_flush(self.inner) };
        if exit_code == 0 {
            Ok(())
        } else {
            Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "Can not flush",
            ))
        }
    }
}

impl std::ops::Drop for Writer {
    fn drop(&mut self) {
        unsafe {
            htslib::bgzf_close(self.inner);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use std::io::Write;

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

    #[test]
    fn test_write_plain() {
        let tmp = tempfile::Builder::new()
            .prefix("rust-htslib")
            .tempdir()
            .expect("Cannot create temp dir");
        let out_path = tmp.path().join("test.vcf");
        println!("{:?}", out_path);

        {
            let w_result = Writer::from_path_with_level(&out_path, CompressionLevel::NoCompression);
            if let Err(ref e) = w_result {
                println!("w_result is {}", e);
            }
            assert!(w_result.is_ok(), "Create plain file with Bgzip writer");
            assert!(out_path.exists(), "Plain file is created with Bgzip writer");
            let mut w = w_result.unwrap();
            let write_result = w.write_all(CONTENT.as_bytes());
            assert!(
                write_result.is_ok(),
                "Plain file can write with Bgzip writer"
            );
        } // let Writer fall out of scope and implicitly close
        assert!(
            !is_bgzip(&out_path).unwrap(),
            "NoCompression file should not be detected as BGZIP"
        );
        let my_content = std::fs::read_to_string(&out_path).unwrap();
        assert_eq!(
            my_content, CONTENT,
            "Writing bgzip file with no compression"
        );

        tmp.close().expect("Failed to delete temp dir");
    }

    #[test]
    fn test_write_default() {
        let tmp = tempfile::Builder::new()
            .prefix("rust-htslib")
            .tempdir()
            .expect("Cannot create temp dir");
        let out_path = tmp.path().join("test.vcf.bgzf");
        println!("{:?}", out_path);
        {
            let w_result = Writer::from_path(&out_path);
            if let Err(ref e) = w_result {
                println!("w_result is {}", e);
            }
            assert!(w_result.is_ok(), "Create bgzip file with Bgzip writer");
            assert!(
                std::path::Path::new(&out_path).exists(),
                "Bgzip file is created with Bgzip writer"
            );
            let mut w = w_result.unwrap();
            let write_result = w.write_all(CONTENT.as_bytes());
            assert!(
                write_result.is_ok(),
                "Bgzip file can write with Bgzip writer"
            );
        } // let Writer fall out of scope and implicitly close

        // Read in with bgzip reader
        let mut my_content = String::new();
        Reader::from_path(&out_path)
            .unwrap()
            .read_to_string(&mut my_content)
            .unwrap();
        assert_eq!(
            my_content, CONTENT,
            "Writing bgzip file with default compression"
        );

        assert!(
            is_bgzip(&out_path).unwrap(),
            "Default BGZIP file detected as BGZIP"
        );
        tmp.close().expect("Failed to delete temp dir");
    }

    #[test]
    fn test_write_compression_levels() {
        let tmp = tempfile::Builder::new()
            .prefix("rust-htslib")
            .tempdir()
            .expect("Cannot create temp dir");
        let out_path = tmp.path().join("test.vcf.bgzf");

        // Test all levels except NoCompression
        let compression_levels = vec![
            CompressionLevel::Fastest,
            CompressionLevel::Maximum,
            CompressionLevel::Uncompressed,
        ]
        .into_iter()
        .chain((-1..=9_i8).map(|n| CompressionLevel::Level(n)));

        for level in compression_levels {
            {
                let w_result = Writer::from_path_with_level(&out_path, level);
                if let Err(ref e) = w_result {
                    println!("w_result is {}", e);
                }
                assert!(w_result.is_ok(), "Create bgzip file with Bgzip writer");
                assert!(
                    std::path::Path::new(&out_path).exists(),
                    "Bgzip file is created with Bgzip writer"
                );
                let mut w = w_result.unwrap();
                let write_result = w.write_all(CONTENT.as_bytes());
                assert!(
                    write_result.is_ok(),
                    "Bgzip file can write with Bgzip writer"
                );
            } // let Writer fall out of scope and implicitly close

            // Read in with bgzip reader
            let mut my_content = String::new();
            Reader::from_path(&out_path)
                .unwrap()
                .read_to_string(&mut my_content)
                .unwrap();
            assert_eq!(
                my_content, CONTENT,
                "Writing bgzip file with {:?} compression",
                level
            );

            assert!(
                is_bgzip(&out_path).unwrap(),
                "Writing BGZIP file with {:?} compression detected as BGZIP",
                level
            );
        }
        tmp.close().expect("Failed to delete temp dir");
    }

    #[test]
    fn test_write_with_threadpool() {
        let tmp = tempfile::Builder::new()
            .prefix("rust-htslib")
            .tempdir()
            .expect("Cannot create temp dir");
        let out_path = tmp.path().join("test.vcf.bgzf");

        let content = CONTENT.as_bytes();
        println!("{:?}", out_path);
        {
            let w_result = Writer::from_path(&out_path);
            if let Err(ref e) = w_result {
                println!("w_result is {}", e);
            }
            assert!(w_result.is_ok(), "Create bgzip file with Bgzip threadpool");
            assert!(
                std::path::Path::new(&out_path).exists(),
                "Bgzip file is created with Bgzip threadpool"
            );

            let mut w = w_result.unwrap();
            let tpool_result = ThreadPool::new(5);
            assert!(tpool_result.is_ok(), "Creating thread pool");
            let tpool = tpool_result.unwrap();

            let set_tpool_result = w.set_thread_pool(&tpool);
            assert!(set_tpool_result.is_ok(), "Setting thread pool");

            let write_result = w.write_all(content);
            assert!(
                write_result.is_ok(),
                "Bgzip file can write with Bgzip threadpool"
            );
        } // let Writer fall out of scope and implicitly close

        // Read in with bgzip reader
        let mut my_content = String::new();
        Reader::from_path(&out_path)
            .unwrap()
            .read_to_string(&mut my_content)
            .unwrap();
        assert_eq!(my_content, CONTENT, "Writing bgzip file with threadpool");

        assert!(
            is_bgzip(&out_path).unwrap(),
            "Threadpool BGZIP file detected as BGZIP"
        );

        tmp.close().expect("Failed to delete temp dir");
    }
}
