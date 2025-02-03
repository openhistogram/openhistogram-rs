//! An implementation of the OpenHistogram log-linear histogram specification.
//! 
//! Use this crate for high-performance, high-volume measurement recording and analysis.
//!
//! # Examples
//! 
//! ```ignore
//! // Create a new empty histogram
//! let mut perf = openhistogram::Histogram::new();
//! 
//! // Instrument a function
//! fn some_interesting_function() {
//!     let start = std::time::Instant::now();
//!     // do some work
//!     let finish = std::time::Instant::now();
//!     perf.insert(finish.duration_since(start), 1);
//! }
//! 
//! // Analyze
//! println!("Average call time: {}s, 99th-percentile: {}s", perf.mean(), perf.quantile1(&[0.99])?[0]);
//! ```
mod tables;
mod bin;
#[macro_use]
mod histogram;
mod serde;
mod serialize;
mod analysis;
use thiserror;

pub use bin::Bin;
pub use bin::NAN_BIN;
pub use histogram::Histogram;
pub use analysis::QuantileType;

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("Invalid Bin")]
    InvalidBin,
    #[error("Corrupt encoding: {}", .0)]
    CorruptEncoding(String),
    #[error("serialization error")]
    WriteError(#[from] std::io::Error),
    #[error("quantiles out of order")]
    QuantilesOutOfOrder,
    #[error("quantile outside of range [0,1]")]
    QuantileOutOfBounds,
}