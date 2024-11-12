use std::cell::RefCell;
use std::sync::Arc;

pub use crate::errors::{Error, Result};
use crate::htslib;

/// An HTSlib thread pool. Create a thread pool and use `set_thread_pool()` methods
/// to share a thread pool across multiple BAM readers & writers.
/// The Rust wrapper holds the htslib thread pool behind a Rc, and a Rc reference
/// to the thread pool is held by each reader / writer so you don't need to
/// explicitly manage the lifetime of the `ThreadPool`.
#[derive(Clone, Debug)]
pub struct ThreadPool {
    pub(crate) handle: Arc<RefCell<InnerThreadPool>>,
}

impl ThreadPool {
    /// Create a new thread pool with `n_threads` threads.
    pub fn new(n_threads: u32) -> Result<ThreadPool> {
        let ret = unsafe { htslib::hts_tpool_init(n_threads as i32) };

        if ret.is_null() {
            Err(Error::ThreadPool)
        } else {
            let inner = htslib::htsThreadPool {
                pool: ret,
                // this matches the default size
                // used in hts_set_threads.
                qsize: n_threads as i32 * 2,
            };
            let inner = InnerThreadPool { inner };

            let handle = Arc::new(RefCell::new(inner));
            Ok(ThreadPool { handle })
        }
    }
}

/// Internal htsThreadPool
#[derive(Clone, Debug)]
pub struct InnerThreadPool {
    pub(crate) inner: htslib::htsThreadPool,
}

impl Drop for InnerThreadPool {
    fn drop(&mut self) {
        if !self.inner.pool.is_null() {
            unsafe {
                htslib::hts_tpool_destroy(self.inner.pool);
            }
        }

        self.inner.pool = std::ptr::null_mut();
    }
}
