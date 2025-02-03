use crate::bin::NAN_BIN;

use super::bin::Bin;
use std::collections::BTreeMap;

#[allow(unused_macros)]
macro_rules! bv {
    (($bin:expr, $count:expr)) => ( (Into::<Bin>::into($bin), $count) );
    (($bin:expr)) => ( (Into::<Bin>::into($bin), 1) );
    ($bin:expr) => ( (Into::<Bin>::into($bin), 1) )
}
#[allow(unused_macros)]
/// hist! macro will turn an list of values into a histogram.
/// 
/// The values can be of the form `value`, `(value)`, or `(value, sample_count)`.
/// In the first two forms, the omitted `sample_count` defaults to 1.
macro_rules! hist {
    () => ( Histogram::new() );
    ($($bv:tt),+ $(,)?) => ( {
        let mut h = Histogram::new();
        [$(bv![$bv]),+].iter().for_each(|x| { h.insert(x.0, x.1); });
        h
    } )
}

#[derive(Debug,Clone)]
/// A distribution histogram containing sample counts in base 10 log-linear [Bin]s.
/// 
/// The [Histogram] is designed to be a fast, low-memory-overhead data
/// structure for storing high-frequency performance data.  Several
/// statistical/analytical approximations are provided.
pub struct Histogram {
    pub(crate) bvs: BTreeMap<Bin,u64>,
    pub(crate) nan_count: u64,
}

impl Histogram {
    /// This function creates a new, empty histogram.
    pub fn new() -> Histogram {
        Histogram { bvs: BTreeMap::<Bin,u64>::new(), nan_count: 0u64 }
    }
    /// This function will remove any [Bin]s in the histogram with no
    /// samples.
    fn sweep(&mut self) {
        self.bvs.retain(|_, count| count != &0u64)
    }
    /// This function will record `count` samples in [Bin] bin.
    /// 
    /// If the samples in a bit exceed `u64::MAX`, they will saturate the
    /// bin and be considered `u64::MAX`.
    /// 
    /// ```
    /// use openhistogram::Histogram;
    /// use openhistogram::Bin;
    /// 
    /// let mut h = Histogram::new();
    /// h.insert(Bin::from_int_scale(1, -9).into(), 1);
    /// h.insert((1, -9).into(), 1);
    /// h.insert((1e-9).into(), 1);
    /// ```
    pub fn insert(&mut self, bin: Bin, count: u64) {
        if bin == NAN_BIN {
            self.nan_count = self.nan_count.saturating_add(count);
        } else {
            self.bvs.entry(bin)
                .and_modify(|v| { *v = (*v).saturating_add(count); })
                .or_insert(count);
        }
    }
    /// This function will add (or remove) samples in the specified [Bin].
    ///
    /// If the samples in a bit exceed `u64::MAX`, they will saturate the
    /// bin and be considered `u64::MAX`.  If `count` is negative and of
    /// greater magnitude than the current sample count of the bin in question,
    /// the bin will be reduced to 0 samples.
    pub fn insert_signed(&mut self, bin: Bin, count: i64) {
        if bin == NAN_BIN {
            self.nan_count = self.nan_count.saturating_add_signed(count);
        } else {
            if self.bvs.entry(bin)
                .and_modify(|v| { *v = (*v).saturating_add_signed(count); })
                .or_insert(0u64.saturating_add_signed(count)) == &0u64 {
                self.sweep()
            }
        }
    }
    /// This function will remove samples from the specifiec [Bin].
    ///
    /// If `count` is of /// greater magnitude than the current sample count
    /// of the bin in question, the bin will be reduced to 0 samples.
    pub fn remove(&mut self, bin: Bin, count: u64) {
        if bin == NAN_BIN {
            self.nan_count = self.nan_count.saturating_sub(count);
        } else {
            if let Some(v) = self.bvs.get_mut(&bin) {
                *v = (*v).saturating_sub(count);
                if *v == 0 {
                    self.sweep()
                }
            }
        }
    }
    /// This function will accumulate an array of [Histogram]s onto the caller.
    /// 
    /// The sample count will saturate a `u64` instead of wrapping.
    pub fn accumulate_many(&mut self, others: &[&Self]) {
        others.iter().for_each(|other| {
            other.bvs.iter().for_each(|(bin, count)| {
                self.insert(*bin, *count);
            });
            self.nan_count = self.nan_count.saturating_add(other.nan_count);
        });
    }
    /// This function will accumulate a [Histogram] onto the caller.
    /// 
    /// The sample count will saturate a `u64` instead of wrapping.
    pub fn accumulate(&mut self, other: &Self) {
        self.accumulate_many(&[other]);
    }

    /// This function will decumulate the aggregate of `others` [Histogram]s from the caller.
    /// 
    /// Each [Histogram] in `others` will be removed from the caller. If
    /// an attempt is made to remove more samples than exist a given [Bin] 
    /// in the calling histogram, that [Bin] will have its sample count set
    /// to 0.
    pub fn decumulate_many(&mut self, others: &[&Self]) {
        if others.iter().fold(false, |zeroed, other| {
            self.nan_count = self.nan_count.saturating_sub(other.nan_count);
            let z = other.bvs.iter().fold(false, |zeroed, (bin, count)| {
                if let Some(v) = self.bvs.get_mut(&bin) {
                    *v = (*v).saturating_sub(*count);
                    zeroed || (*v == 0)
                } else {
                    zeroed
                }
            });
            zeroed || z
        }) {
            self.sweep();
        }
    }
    /// This function will decumulate a [Histogram] from the caller.
    /// 
    /// The sample count will saturate at 0.
    pub fn decumulate(&mut self, other: &Self) {
        self.decumulate_many(&[&other]);
    }
    pub fn subtract(&self, other: &Self) -> Histogram {
        let mut copy = self.clone();
        copy.decumulate(other);
        copy
    }
    /// This function will return the total number of samples in all [Bin]s
    /// excluding the `NAN_BIN`.
    pub fn count(&self) -> u128 {
        self.bvs.iter().fold(0u128, |total, (_, count)| { total + (*count) as u128 })
    }
    /// This function will return the total number of samples in all [Bin]s
    /// including the `NAN_BIN`.
    pub fn count_including_nan(&self) -> u128 {
        self.count() + self.nan_count as u128
    }
    /// This function will zero the sample count in all [Bin]s less
    /// than `lower` and greater than `upper`.
    pub fn clamp(&mut self, lower: &Bin, upper: &Bin) {
        self.bvs.retain(|bin, _| {
            println!("{} >= {} ({}) && {} <= {} ({})", bin, lower, bin>=lower, bin, upper, bin<=upper);
            bin >= lower && bin <= upper
        });
    }
    /// This function will reduce the sample counts in all [Bin]s to
    /// approximately `factor` of their original magnitutdes.
    /// 
    /// A binomial reduction is used and the results will not be exact.
    pub fn downsample(&mut self, factor: f64) {
        let factor = factor.clamp(0.0, 1.0);
        self.bvs.iter_mut().for_each(|(_, count)| {
            *count = super::analysis::binomial_reduce(*count, factor);
        });
    }
    /// This function will empty the histogram of all recorded samples.
    pub fn clear(&mut self) {
        self.bvs.clear();
    }
}

impl PartialEq for Histogram {
    fn eq(&self, other: &Self) -> bool {
        let i1 = self.bvs.iter().filter(|x| x.1 != &0u64);
        let mut i2 = other.bvs.iter().filter(|x| x.1 != &0u64);
        if self.nan_count != other.nan_count {
            return false;
        }
        for left in i1 {
            match i2.next() {
                Some(right) => {
                    if left.0 != right.0 || left.1 != right.1 {
                        return false;
                    }
                },
                None => {
                    return false;
                }
            }
        }
        i2.next() == None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::prelude::*;

    #[test]
    fn basic() {
        let mut h1 = Histogram::new();
        let mut h2 = Histogram::new();
        assert_eq!(h1, h2);
        h1.insert(12.into(), 10);
        h1.insert(NAN_BIN, 1);
        h2.insert(12.into(), 10);
        h2.insert(NAN_BIN, 1);
        assert_eq!(h1, h2);
        assert_eq!(10, h1.count());
        h1.clear();
        assert_eq!(0, h1.count());
    }

    #[test]
    fn downsample() {
        let mut rng = rand::rng();
        let mut total = 0u128;
        let mut h = hist![];
        for i in 10 .. 100 {
            for j in -9 .. 10 {
                let cnt = rng.random::<u64>() % 4;
                h.insert(Bin::from_int_scale(i,j), cnt);
                total += cnt as u128;
            }
        }
        assert_eq!(total, h.count());
        h.downsample(0.2);
        let actual = h.count() as f64;
        let expected = total as f64 * 0.2;
        let err = ((actual - expected) / total as f64).abs();
        assert!(err < 0.05);
    }

    #[test]
    fn clone() {
        let h = hist![(13,12),(1.3e5,2),(-1.0),0,(f64::NAN,3)];
        let h2 = h.clone();
        assert_eq!(h, h2);
    }

    #[test]
    fn sweep() {
        let h1 = hist![(13,12),(1.3e5,2),0,(f64::NAN,3)];
        let mut h2= hist![(13,12),(1.3e5,2),(-1.0),0,(f64::NAN,3)];
        h2.remove((-1.0).into(), 1);
        assert_eq!(h1, h2);
    }

    #[test]
    fn clamp() {
        let mut h1 = hist![(-10.0), (-9.9), 0, 99, 110];
        h1.clamp(&(-9.9).into(), &(99.0).into());
        assert_eq!(h1, hist![(-9.9),0,99]);
    }

    #[test]
    fn sub_accum() {
        let mut rng = rand::rng();
        let mut hists: [Histogram; 10] = core::array::from_fn(|_| Histogram::new());
        let mut samples = 0u128;
        for i in 0..10 {
            if i == 8 {
                continue;
            }
            for _ in 0..100 {
                let v: Bin = ((rng.random::<u64>() % 100) as f64 + 10.0).into();
                hists[i].insert(v, 1);
                samples += 1;
            }
        }
        let mut h = Histogram::new();
        h.accumulate_many(&hists.each_ref());
        assert_eq!(h.count(), samples);
        for i in 0..9 {
            h.decumulate(&hists[i]);
        }
        assert_eq!(hists[9].count(), h.count());
    }

    #[test]
    fn sample_saturate() {
        let mut h = hist![];
        h.insert(1.into(), u64::MAX);
        h.insert(1.into(), 1);
        h.insert(2.into(), u64::MAX);
        h.insert(2.into(), 2);
        h.insert(3.into(), 1);
        h.remove(3.into(), 3);
        assert_eq!(h.count_nearby(1.0), u64::MAX);
        assert_eq!(h.count_nearby(2.0), u64::MAX);
        assert_eq!(h.count_nearby(3.0), 0);
    }

    #[test]
    fn hist_macro() {
        let h = hist![(12,1),(10),(1e13,3),14];
        assert_eq!(h.count(), 6);
    }
}