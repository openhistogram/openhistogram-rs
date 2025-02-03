use super::bin::Bin;
use super::Histogram;
use std::cmp::Ordering;
use probability::prelude::{Binomial, Distribution};
use rand::prelude::*;


/// The QuantileType represents various methodologies for calculating quantiles.
pub enum QuantileType {
    /// corresponds to Type=1 quantiles in the Hyndman-Fan list (Statistical Computing, 1996).
    Type1,
    /// corresponds to discretized Type=7 quantiles in the Hyndman-Fan list (Statistical Computing, 1996).
    Type7,
}

fn binomial_reduce_random(n: u64, pr: f64, tgt: f64) -> u64 {
    if pr <= 0.0 { 0 }
    else if pr >= 1.0 { n }
    else {
        let cumbin = Binomial::new(n as usize, pr);
        let (mut left, mut right) = (0, n);
        let mut i = (n as f64 * pr) as u64;
        let mut bias = match (n as f64 * 0.005) as u64 {
                x if x < 2 => 2,
                x => x,
        };
        while left < right {
            let s = i as f64;
            let cum = cumbin.distribution(s);
            if i == right {
                if right - left  == 1 {
                    if tgt <= cum {
                        return i;
                    }
                    else {
                        return i-1;
                    }
                }
            }
            if tgt <= cum {
                if i == 0 {
                    return 0;
                }
                right = i;
            }
            else {
                if i == n {
                    return n;
                }
                else {
                    left = i;
                }
            }

            let mut skip = (right - left) / bias;
            if bias != 2 && skip < 1 {
                skip = (right - left) / 2;
            }
            if skip < 1 {
                skip = 1;
            }

            if tgt <= cum {
                i = right - skip;
            }
            else {
                i = left + skip;
            }
            bias = 2;
        }
        0
    }
}
pub(crate) fn binomial_reduce(n: u64, pr: f64) -> u64 {
    let mut rng = rand::rng();
    let tgt = rng.random::<f64>();
    binomial_reduce_random(n, pr, tgt)
}

impl Histogram {
    fn accum_moments<const N: usize>(&self, ks: &[i32; N]) -> [f64; N] {
        self.bvs.iter().fold([0.0; N], |v, (bin, count)| {
            let midpoint = bin.midpoint();
            let cardinality = *count as f64;
            let mut newv = [0.0; N];
            for i in 0..N {
                newv[i] = v[i] + midpoint.powi(ks[i]) * cardinality;
            }
            newv
        })
    }
    /// Calculate an approximate mean across all samples.
    pub fn mean(&self) -> f64 {
        let r = self.accum_moments(&[0, 1]);
        if r[0] == 0.0 {
            f64::NAN
        } else {
            r[1] / r[0]
        }
    }
    /// Calculate an approximate sum across all samples.
    pub fn sum(&self) -> f64 {
        self.accum_moments(&[1])[0]
    }
    /// Calculate an approximate k-th moment across all samples.
    pub fn moment(&self, k: i32) -> f64 {
        let r = self.accum_moments(&[0, k]);
        if r[0] == 0.0 {
            f64::NAN
        } else {
            r[1] / r[0].powi(k)
        }
    }
    /// Calculate the approximate standard deviation across all samples.
    pub fn stddev(&self) -> f64 {
        let r = self.accum_moments(&[0, 1, 2]);
        if r[0] == 0.0 {
            f64::NAN
        } else {
            (r[2] / r[0] - (r[1] / r[0]).powi(2)).sqrt()
        }
    }

    /// Calculate an approximate number of samples relative to a `pivot` [Bin]
    pub fn count_relative(&self, pivot: &Bin, ord: Ordering, inclusive: bool) -> u128 {
        assert_ne!(ord, Ordering::Equal);
        self.bvs.iter().fold(0u128, |v, (bin, count)| {
            let rel = bin.cmp(pivot);
            if rel == Ordering::Equal {
                if inclusive {
                    v + (*count as u128)
                } else {
                    v
                }
            } else if rel == ord {
                v + *count as u128
            } else {
                v
            }
        })
    }
    /// Count all samples in all [Bin]s less than or equal to the [Bin] containing `v`.
    pub fn count_below_inclusive(&self, v: f64) -> u128 {
        self.count_relative(&v.into(), Ordering::Less, true)
    }
    /// Count all samples in all [Bin]s less than the [Bin] containing `v`.
    pub fn count_below_exclusive(&self, v: f64) -> u128 {
        self.count_relative(&v.into(), Ordering::Less, false)
    }
    /// Count all samples in all [Bin]s greater than or equal to the [Bin] containing `v`.
    pub fn count_above_inclusive(&self, v: f64) -> u128 {
        self.count_relative(&v.into(), Ordering::Greater, true)
    }
    /// Count all samples in all [Bin]s greater than the [Bin] containing `v`.
    pub fn count_above_exclusive(&self, v: f64) -> u128 {
        self.count_relative(&v.into(), Ordering::Greater, false)
    }
    /// Count the samples in the [Bin] containing `v`
    pub fn count_nearby(&self, v: f64) -> u64 {
        self.bvs
            .get(&v.into())
            .and_then(|x| Some(*x))
            .or(Some(0u64))
            .unwrap()
    }

    /// Calculate a set of approximate quantiles from the samples.
    pub fn quantile_ex<const N: usize>(
        &self,
        qs: &[f64; N],
        qt: QuantileType,
    ) -> Result<[f64; N], super::Error> {
        let mut qout = [f64::NAN; N];
        if N == 0 {
            return Ok([0.0; N]);
        }
        for i in 1..N {
            if qs[i - 1] > qs[i] {
                return Err(super::Error::QuantilesOutOfOrder);
            }
        }
        if qs[0] < 0.0 || qs[N - 1] > 1.0 {
            return Err(super::Error::QuantileOutOfBounds);
        }

        match self.count() {
            0 => Ok(qout),
            cu128 => {
                let total_count = cu128 as f64;
                /* We use q_out as temporary space to hold the count-normalized quantiles */
                match qt {
                    QuantileType::Type1 => {
                        for i in 0..N {
                            qout[i] = total_count * qs[i]
                        }
                    }
                    QuantileType::Type7 => {
                        for i in 0..N {
                            qout[i] = ((total_count - 1.0) * qs[i] + 1.0).floor()
                        }
                    }
                };

                /* this is to track */
                struct Track {
                    bin_width: f64,
                    bin_left: f64,
                    bin_count: f64,
                    lower_cnt: f64,
                    upper_cnt: f64,
                }
                impl Track {
                    fn update(&mut self, bin: &Bin, count: u64) {
                        self.bin_width = bin.width_signed();
                        self.bin_left = bin.left();
                        self.bin_count = count as f64;
                        self.lower_cnt = self.upper_cnt;
                        self.upper_cnt = self.lower_cnt + count as f64;
                    }
                }
                let mut track = Track {
                    bin_width: 0.0,
                    bin_left: 0.0,
                    bin_count: 0.0,
                    lower_cnt: 0.0,
                    upper_cnt: 0.0,
                };
                // iterate through the tree and stop to act once we've place
                let mut bv_iter = self.bvs.iter();
                if let Some((bin, count)) = bv_iter.next() {
                    track.update(bin, *count);
                }

                for i_q in 0..N {
                    while track.upper_cnt < qout[i_q] {
                        match bv_iter.next() {
                            Some((bin,count)) => track.update(bin, *count),
                            _ => break
                        }
                    }
                    if track.bin_width == 0.0 {
                        qout[i_q] = track.bin_left;
                    }
                    else {
                        match qt {
                            // For `QuantileType::Type1` quantiles:`
                            // A q-quantile for the bucket, will be represented by the sample number:
                            //
                            // (q = 0)  k = 1
                            // (q > 0)  k = ceil(q*n)
                            // 
                            // so that q=0 => k=1 and q=1 => k=n. This corresponds to Type=1 quantiles
                            // in the Hyndman-Fan list (Statistical Computing, 1996).
                            QuantileType::Type1 => {
                                let qn = qout[i_q];
                                assert!(qn >= track.lower_cnt);
                                assert!(qn <= track.upper_cnt);
                                let k = match (qn - track.lower_cnt).ceil() {
                                    0.0 => 1.0,
                                    k => k,
                                };
                                qout[i_q] = track.bin_left + k / (track.bin_count + 1.0) * track.bin_width;
                            },
                            // For Type 7 Quantiles, we consider samples at indices:
                            // 
                            // k = floor( q*(n-1) + 1 )
                            // 
                            // This corresponds to discretized Type=7 quantiles in
                            // the Hyndman-Fan list (Statistical Computing, 1996).
                            QuantileType::Type7 => {
                                let k = qout[i_q] - track.lower_cnt;
                                qout[i_q] = track.bin_left + k / (track.bin_count + 1.0) * track.bin_width;
                            },
                        }
                    }
                }
                Ok(qout)
            }
        }
    }
    /// Calculate a set of Type=1 quantiles from the sample set.
    pub fn quantile1<const N: usize>(
        &self,
        qs: &[f64; N],
    ) -> Result<[f64; N], super::Error> {
        self.quantile_ex(qs, QuantileType::Type1)
    }
    /// Calculate a set of Type=7 quantiles from the sample set.
    pub fn quantile7<const N: usize>(
        &self,
        qs: &[f64; N],
    ) -> Result<[f64; N], super::Error> {
        self.quantile_ex(qs, QuantileType::Type7)
    }
    /// Calculate a set of inverse quantiles from the sample set.
    pub fn inverse_quantile<const N: usize>(
        &self,
        vs: &[f64; N],
    ) -> Result<[f64; N], super::Error> {
        let mut out = [f64::NAN; N];
        if N == 0 {
            return Ok(out);
        }
        for i in 1..N {
            if vs[i] < vs[i-1] {
                return Err(super::Error::QuantilesOutOfOrder);
            }
        }
        let total_count = match self.count() {
            0 => { return Ok(out); },
            x => { x as f64 },
        };
        let mut vs_idx = 0usize;
        let mut threshold = vs[vs_idx];
        let mut count_below = 0u128;
        for (bin, count) in self.bvs.iter() {
            let (bin_lower, bin_upper, bin_width) = match bin.val {
                v if v < 0i8 => {
                    let (b, w) = (bin.absolute_min(), bin.width());
                    (b-w, b, w)
                },
                v if v > 0i8 => {
                    let (b, w) = (bin.absolute_min(), bin.width());
                    (b, b+w, w)
                },
                _ => (-1e-128,1e-128,0.0),
            };

            while threshold < bin_lower {
                out[vs_idx] = count_below as f64 / total_count as f64;
                vs_idx += 1;
                if vs_idx >= N {
                    return Ok(out);
                }
                threshold = vs[vs_idx];
            }
            while threshold < bin_upper {
                out[vs_idx] = if bin_width > 0.0 {
                        let pr = (threshold - bin_lower) / bin_width;
                        (count_below as f64 + pr * (*count) as f64) / total_count as f64
                    } else {
                        count_below as f64 / total_count as f64
                    };
                vs_idx += 1;
                if vs_idx >= N {
                    return Ok(out);
                }
                threshold = vs[vs_idx];
            }
            count_below += *count as u128;
        }
        while vs_idx < N {
            out[vs_idx] = 1.0;
            vs_idx += 1;
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts() {
        let h = hist![
            (-10, 2),
            (-7, 3),
            (-0.1, 5),
            (0, 1),
            (3),
            (7, 14),
            (f64::NAN)
        ];
        assert_eq!(h.count_including_nan(), 27);
        assert_eq!(h.count(), 26);
        assert_eq!(h.count_below_inclusive(-0.1), 10);
        assert_eq!(h.count_above_inclusive(-0.1), 21);
        assert_eq!(h.count_below_exclusive(-0.1), 5);
        assert_eq!(h.count_above_exclusive(-0.1), 16);
        assert_eq!(h.count_nearby(-0.1), 5);
        assert_eq!(h.count_nearby(0.0), 1);
    }
    #[test]
    fn quantiles() -> Result<(), super::super::Error> {
        let h = hist![0.123, 0, 0.43, 0.41, 0.415, 0.2201, 0.3201, 0.125, 0.13];
        let qs = [0.0,0.95,0.99,1.0];
        assert_eq!(h.quantile1(&qs)?, [0.0, 0.435, 0.435, 0.435]);
        assert_eq!(h.quantile7(&qs)?, [0.0, 0.41666666666666663, 0.41666666666666663, 0.435]);
        assert_eq!(h.mean(), 0.2443387586038558);
        assert_eq!(h.sum(), 2.199048827434702);

        let h = hist![1,1];
        let a = 1.0+0.1*1.0/3.0;
        let b = 1.0+0.1*2.0/3.0;
        let qs = [ 0.0, 0.25, 0.5, 0.75, 1.0 ];
        assert_eq!(h.quantile_ex(&qs, QuantileType::Type1)?, [a,a,a,b,b]);
        assert_eq!(h.quantile_ex(&qs, QuantileType::Type7)?, [a,a,a,a,b]);

        let h = hist![1];
        let qs = [ 0.0, 0.25, 0.5, 1.0 ];
        assert_eq!(h.quantile1(&qs)?, [1.05, 1.05, 1.05, 1.05]);
        assert_eq!(h.quantile7(&qs)?, [1.05, 1.05, 1.05, 1.05]);

        let h = hist![1,1,1];
        let qs = [ 0.0, 0.5, 1.0 ];
        assert_eq!(h.quantile1(&qs)?, [1.025, 1.05, 1.075]);
        assert_eq!(h.quantile7(&qs)?, [1.025, 1.05, 1.075]);

        let h = hist![1.0, 2.0];
        assert_eq!(h.quantile1(&[0.5])?, [1.05]);
        assert_eq!(h.quantile7(&[0.5])?, [1.05]);

        let h = hist![1.0, 1e200];
        assert_eq!(h.quantile1(&[0.0,1.0])?, [1.05,1.05]);
        assert_eq!(h.quantile7(&[0.0,1.0])?, [1.05,1.05]);
        assert_eq!(h.mean(), 1.0476190476190477);

        let h = hist![(1e200,3),(0.0,2),(1e-20,3),1e-10];
        assert_eq!(h.quantile1(&[0.0,1.0])?, [0.0, 1.05e-10]);
        assert_eq!(h.quantile7(&[0.0,1.0])?, [0.0, 1.05e-10]);

        let h = hist![0,1];
        let qs = [0.0,0.1,0.499,0.501,0.9,1.0];
        assert_eq!(h.quantile1(&qs)?, [0.0,0.0,0.0,1.05,1.05,1.05]);
        assert_eq!(h.quantile7(&qs)?, [0.0,0.0,0.0,0.0,0.0,1.05]);

        let h = hist![10,100];
        let qs = [0.0,0.1,0.499,0.501,0.9,1.0];
        assert_eq!(h.quantile1(&qs)?, [10.5,10.5,10.5,105.0,105.0,105.0]);
        assert_eq!(h.quantile7(&qs)?, [10.5,10.5,10.5,10.5,10.5,105.0]);
        Ok(())
    }

    #[test]
    fn inverse_quantiles() -> Result<(), super::super::Error> {
        use approx::*;

        let vs = [-200.0,-100.0,0.0,1.0,1.001,1.1,1.2,2.0,3.0,4.0];
        let h = hist![];
        let rq = h.inverse_quantile(&vs)?;
        assert!(rq.iter().fold(true, |v, x| v && x.is_nan()));

        let h = hist![(-100)];
        assert_eq!([0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0], h.inverse_quantile(&vs)?);

        let h = hist![0];
        assert_eq!([0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0], h.inverse_quantile(&vs)?);

        let h = hist![1,2,3];
        [0.0,0.0,0.0,0.0,(1.0/3.0)/100.0,1.0/3.0,1.0/3.0,1.0/3.0,2.0/3.0,1.0]
            .iter().zip(h.inverse_quantile(&vs)?
            .iter()).for_each(|x| {
                assert_relative_eq!(x.0, x.1, epsilon = 0.00000001);
            });
        Ok(())
    }
}
