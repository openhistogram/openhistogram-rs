use super::tables::POWER_OF_TEN;

use std::{fmt, time::Duration};

macro_rules! bin {
    ($v:expr,$e:expr) => {
        Bin{ val:$v, exp:$e }
    };
}

/// A [Bin] represent a log(10)-linear range in some unitless value domain.
///
/// A bin can represent:
///  * NAN (not a representable number),
///  * \[0\], and
///  * \[x.y * 10^z, (x.y + 0.1) * 10^z) where:
///    * x is in the positive and negative integral range \[1,9\]
///    * y is in the integral range \[0,9\]
///    * z is in the integral range \[-128,127\]
/// 
#[derive(Eq, PartialEq, Debug, Copy, Clone)]
pub struct Bin {
    /// `val` is represented by an i8, but excepting the NaN representation must
    /// adhere to the domain (-100,-10] or [10,100) which represt the positive and
    /// negative numbers 1.0 through 9.9 in 0.1 incremements.
    pub(crate) val: i8,
    /// `exp` is the base-10 exponent by which `val` is multiplied to represent the
    /// bucket boundary.  As an i8, the exponent is limited to the range [-128,127].
    pub(crate) exp: i8
}

/// `NAN_BIN` is a solitary [Bin] for accumulations on values that cannot be expressed
/// in the OpenHistogram log-linear binning system.  This can happen because the value
/// being inserted into this histogram falls outside of a recognized range.
pub const NAN_BIN: Bin = bin![-1,0];
/// `ZERO_BIN` is the special zero-width, single-value [Bin] representing the number zero.
pub const ZERO_BIN: Bin = bin![0,0];

impl fmt::Display for Bin {
    /// This function implements the `fmt` function for `fmt:Display` and
    /// formats the Bin as its absolute minimum in the exponential notation
    /// (e.g. -3.2e-3)
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_nan() {
            write!(f, "NaN")
        }
        else if !self.is_valid() {
            // It should be nigh impossible to experience this.
            write!(f, "Invalid")
        }
        else {
            write!(f, "{:.1}e{:}", self.val as f64 / 10.0, self.exp)
        }
    }
}

impl Ord for Bin {
    /// Implements ordering of Bins based on their absolute minimum boundary.
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if self == other {
            std::cmp::Ordering::Equal
        } else {
            if self.is_nan() {
                return std::cmp::Ordering::Greater;
            }
            if other.is_nan() {
                return std::cmp::Ordering::Less;
            }
            if self.val == 0 || other.val == 0 {
                return self.val.cmp(&other.val);
            }
            if (self.val < 0 && other.val > 0) || (self.val > 0 && other.val < 0) {
                return self.val.cmp(&other.val);
            }
            if self.exp == other.exp {
                return self.val.cmp(&other.val);
            }
            if self.exp > other.exp {
                if self.val < 0 {
                    std::cmp::Ordering::Less
                }
                else {
                    std::cmp::Ordering::Greater
                }
            }
            else {
                if self.val < 0 {
                    std::cmp::Ordering::Greater
                }
                else {
                    std::cmp::Ordering::Less
                }
            }
        }
    }
}

impl PartialOrd for Bin {
    /// This function partially orders such that NAN Bins do not have ordering.
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self.is_nan() || other.is_nan() {
            None
        }
        else {
            Some(self.cmp(&other))
        }
    }
}

impl Bin {
    /// Determines if the bin is the NAN bin.
    pub fn is_nan(&self) -> bool {
        self.val == NAN_BIN.val && self.exp == NAN_BIN.exp
    }
    /// This function validates the domain for the `val` field.
    pub fn is_valid(&self) -> bool {
        self.is_nan() ||
        (self.val == 0 && self.exp == 0) ||
        (self.val < 100 && self.val >= 10) ||
        (self.val > -100 && self.val <= -10)
    }
    /// This function returns a new [Bin] containing the supplied `f64`
    pub fn from_f64(v: f64) -> Bin {
        if v.is_infinite() || v.is_nan() {
            NAN_BIN
        }
        else if v == 0.0 {
            ZERO_BIN
        }
        else {
            let neg = v.is_sign_negative();
            let mut v = v.abs();
            let big_exp = v.log10().floor() as i32;
            let mut hb_exp = big_exp as i8;
            if hb_exp as i32 != big_exp {
                // rolled.
                if big_exp > 0 {
                    NAN_BIN
                }
                else {
                    ZERO_BIN
                }
            }
            else {
                let pot_idx = hb_exp as u8;
                v /= POWER_OF_TEN[pot_idx as usize];
                v *= 10.0;
                // avoid rounding problem at the bucket boundary
                // e.g. d=0.11 results in hb.val = 10 (should be 11)
                // by allowing an error margin (in the order or magnitude
                // of the expected rounding errors of the above transformations)
                let mut hb_val = ((v + 1e-13).floor() * if neg { -1.0 } else { 1.0 }) as i8;
                if hb_val == 100 || hb_val == -100 {
                    if hb_exp < 127 {
                        hb_val /= 10;
                        hb_exp +=1;
                    }
                    else {
                        return NAN_BIN;
                    }
                }
                if hb_val == 0 {
                    ZERO_BIN
                }
                else if !((hb_val >= 10 && hb_val < 100) || (hb_val <= -10 && hb_val > -100)) {
                    NAN_BIN
                }
                else {
                    Bin { val: hb_val, exp: hb_exp }
                }
            }
        }
    }
    /// This function creates a new Bin from a Duration (as seconds)
    pub fn from_duration(d: Duration) -> Self {
        Bin::from_int_scale(d.as_nanos() as i128, -9)
    }
    /// This function creates a new bin from a decomposition of integral values.
    /// 
    /// `value` * 10 ^ 'scale`, made negative is `positive` is `false`
    pub fn from_uint_scale(positive: bool, value: u128, scale: i16) -> Self {
        if value == 0 { return ZERO_BIN }
        let sign = if positive { 1i8 } else { -1i8 };
        let mut value = value;
        let mut exp = scale + 1; // this is b/c `Bin.val` is *10
        if value < 10 {
            value = value * 10;
            exp = exp - 1;
        }
        while value >= 100 {
            value = value / 10;
            exp = exp + 1;
        }
        if exp < -128 {
            ZERO_BIN
        }
        else if exp > 127 {
            NAN_BIN
        }
        else {
            assert!(value >= 10 && value < 100);
            if let Some(val) = sign.checked_mul(value as i8) {
                Bin { val: val, exp: exp as i8}
            } else {
                NAN_BIN
            }
        }
    }
    /// This function creates a new Bin from a scaled integer.
    /// 
    /// `value` * 10 ^ `scale`
    pub fn from_int_scale(value: i128, scale: i16) -> Self {
        if value == 0 { return ZERO_BIN }
        let (positive, unsigned_value) =
            if value >= 0 {
                (true, value as u128)
            } else {
                if value == std::i128::MIN {
                    (false, (std::i128::MAX - 1) as u128)
                } else {
                    (false, (-value) as u128)
                }
            };
        Bin::from_uint_scale(positive, unsigned_value, scale)
    }
    /// This function calculates the signed width of a Bin.
    /// 
    /// If the bin is exist in the negative, the width will be negative.
    /// If the width cannot be calculated, `f64::NAN` is returned.
    pub fn width_signed(&self) -> f64 {
        if self.is_nan() {
            f64::NAN
        } else {
            if self.val == 0 {
                0.0
            } else {
                let divisor = if self.val.is_negative() { -10.0 } else { 10.0 };
                let exp_idx = self.exp as u8;
                POWER_OF_TEN[exp_idx as usize]/divisor
            }
        }
    }
    /// This function calculates the width of a [Bin].
    /// 
    /// If the width cannot be calculated, `f64::NAN` is returned.
    pub fn width(&self) -> f64 {
        if self.is_nan() {
            f64::NAN
        } else {
            if self.val == 0 {
                0.0
            } else {
                let exp_idx = self.exp as u8;
                POWER_OF_TEN[exp_idx as usize]/10.0
            }
        }
    }
    /// This function returns a error-minimizing (not linear) midpoint of the [Bin].
    /// 
    /// Take the [B0,B1) bin, with a bin-width of B1-B0...
    /// as a sample S approaches B1, we see error ((B1-B0)(1-M))/B1
    /// and as S approaches B0, we see error ((B1-B0)M)/B0.
    ///
    /// M should be chosen such that:
    ///
    /// > ((B1-B0)(1-M))/B1 = ((B1-B0)M)/B0  
    /// > (B0)(B1-B0)(1-M) = (B1)(B1-B0)(M)  
    /// > B0 - B0*M = B1*M  
    /// > B0 = (B0 + B1)(M)  
    /// > M = (B0)/(B0 + B1)  
    pub fn midpoint(&self) -> f64 {
        if self.is_nan() {
            f64::NAN
        } else {
            if self.val == 0 {
                0.0
            } else {
                let bottom: f64 = self.into();
                let interval = if bottom < 0.0 {
                    self.width() * -1.0
                } else {
                    self.width()
                };
                let top = bottom + interval;
                let ratio = bottom / (bottom + top);
                bottom + interval * ratio
            }
        }
    }
    /// This function returns the absolute minimum boundary (closest to 0) of the [Bin]
    pub fn absolute_min(&self) -> f64 {
        if self.is_nan() {
            f64::NAN
        } else {
            if self.val == 0 {
                0.0
            } else {
                let exp_idx = self.exp as u8;
                (self.val as f64 / 10.0) * POWER_OF_TEN[exp_idx as usize]
            }
        }
    }
    /// This function returns the absolute maximum boundary (furthest from 0) of the [Bin]
    pub fn absolute_max(&self) -> f64 {
        if self.is_nan() {
            f64::NAN
        } else {
            self.absolute_min() + self.width_signed()
        }
    }
    /// This function returns the boundary of the [Bin] closest to -Inf
    pub fn left(&self) -> f64 {
        if self.is_nan() {
            f64::NAN
        }
        else if self.val == 0 {
            0.0
        }
        else if self.val > 0 {
            self.absolute_min()
        }
        else {
            self.absolute_min() - self.width()
        }
    }
    /// This function returns the boundary of the [Bin] closest to +Inf
    pub fn right(&self) -> f64 {
        if self.is_nan() {
            f64::NAN
        }
        else if self.val == 0 {
            0.0
        }
        else if self.val < 0 {
            self.absolute_min()
        }
        else {
            self.absolute_min() + self.width()
        }
    }
}

impl From<&Bin> for f64 {
    fn from(bin: &Bin) -> Self {
        bin.absolute_min()
    }
}

impl From<f64> for Bin {
    /// This function converts a f64 into the [Bin] that contains it, if possible.
    fn from(v: f64) -> Self {
        Bin::from_f64(v)
    }
}

impl From<f32> for Bin {
    /// This function converts a f32 into the [Bin] that contains it, if possible.
    fn from(v: f32) -> Self {
        Bin::from_f64(v as f64)
    }
}

impl From<(u128,i8)> for Bin {
    fn from((value,scale): (u128,i8)) -> Self {
        Bin::from_uint_scale(true, value, scale as i16)
    }
}
impl From<(u128,i16)> for Bin {
    fn from((value,scale): (u128,i16)) -> Self {
        Bin::from_uint_scale(true, value, scale)
    }
}
impl From<(u128,i32)> for Bin {
    fn from((value,scale): (u128,i32)) -> Self {
        Bin::from_uint_scale(true, value, scale as i16)
    }
}
impl From<u128> for Bin {
    fn from(v: u128) -> Self {
        Bin::from_uint_scale(true, v, 0)
    }
}

macro_rules! bin_from_int {
    ($itype:ty) => {
        impl From<($itype,i8)> for Bin {
            fn from((value,scale): ($itype,i8)) -> Self {
                Bin::from_int_scale(value as i128, scale as i16)
            }
        }
        impl From<($itype,i16)> for Bin {
            fn from((value,scale): ($itype,i16)) -> Self {
                Bin::from_int_scale(value as i128, scale)
            }
        }
        impl From<($itype,i32)> for Bin {
            fn from((value,scale): ($itype,i32)) -> Self {
                Bin::from_int_scale(value as i128, scale as i16)
            }
        }
        impl From<$itype> for Bin {
            fn from(v: $itype) -> Self {
                Bin::from_int_scale(v as i128, 0)
            }
        }
    };
}

bin_from_int!(i8);
bin_from_int!(u8);
bin_from_int!(i16);
bin_from_int!(u16);
bin_from_int!(i32);
bin_from_int!(u32);
bin_from_int!(i64);
bin_from_int!(u64);
bin_from_int!(i128);

#[cfg(test)]
mod tests {
    use std::time::Duration;

    use super::*;

    fn validate_bin<T>(is: T, v: f64, bin: Bin, s: &str)
    where T: Into::<Bin> {
        let ibin: Bin = is.into();
        assert_eq!(ibin, bin);
        let vbin: Bin = v.into();
        assert_eq!(vbin, bin);
        assert_eq!(vbin.to_string().as_str(), s);
    }

    #[test]
    fn construction() {
        validate_bin((i128::MAX,1), i128::MAX as f64 * 10.0, bin![17,39], "1.7e39");
        validate_bin((i128::MIN,1), i128::MIN as f64 * 10.0, bin![-17,39], "-1.7e39");
        validate_bin((u128::MAX,1), u128::MAX as f64 * 10.0, bin![34,39], "3.4e39");
        validate_bin((100,0), 100.0, bin![10,2], "1.0e2");
        validate_bin((i64::MAX,1), i64::MAX as f64 * 10.0, bin![92,19], "9.2e19");
        validate_bin((i64::MIN,1), i64::MIN as f64 * 10.0, bin![-92,19], "-9.2e19");
        validate_bin((i64::MIN,-127), i64::MIN as f64 * 1e-127, bin![-92,-109], "-9.2e-109");
        validate_bin((i64::MAX,-200), i64::MAX as f64 * 1e-200, ZERO_BIN, "0.0e0");
        validate_bin((10,-128), 1e-127, bin![10,-127], "1.0e-127");
        validate_bin((0,0), 0.0, ZERO_BIN, "0.0e0");
        validate_bin((10, 1), 100.0, bin![10,2], "1.0e2");
        validate_bin((2,0), 2.0, bin![20,0], "2.0e0");
        validate_bin((1,-9), 1e-9, bin![10,-9], "1.0e-9");
        validate_bin((1300000000,-9), 1.3, bin![13,0], "1.3e0");
        validate_bin((-2700,-9), -2.7e-6, bin![-27,-6], "-2.7e-6");
        validate_bin((7,-9), 7e-9, bin![70,-9], "7.0e-9");
        validate_bin((99999,-134), 9.9999e-129, ZERO_BIN, "0.0e0");
        validate_bin((10,-129), 1e-128, bin![10,-128], "1.0e-128");
        validate_bin((100001,-133), 1.00001e-128, bin![10,-128], "1.0e-128");
        validate_bin((109999,-133), 1.09999e-128, bin![10,-128], "1.0e-128");
        validate_bin((11,-129), 1.1e-128, bin![11,-128], "1.1e-128");
        validate_bin((10,126), 1e127, bin![10,127], "1.0e127");
        validate_bin((9999,124), 9.999e127, bin![99,127], "9.9e127");
        validate_bin((1,128), 1e128, NAN_BIN, "NaN");
        validate_bin((-99999,-134), -9.9999e-129, ZERO_BIN, "0.0e0");
        validate_bin((-1,-128), -1e-128, bin![-10,-128], "-1.0e-128");
        validate_bin((-100001,-133), -1.00001e-128, bin![-10,-128], "-1.0e-128");
        validate_bin((-11,-129), -1.1e-128, bin![-11,-128], "-1.1e-128");
        validate_bin((-1,-127), -1e-127, bin![-10,-127], "-1.0e-127");
        validate_bin((-9999,124), -9.999e127, bin![-99,127], "-9.9e127");
        validate_bin((-1,128), -1e128, NAN_BIN, "NaN");
        validate_bin((9999,124), 9.999e127, bin![99,127], "9.9e127");
        validate_bin((172,3), 172000.0, Bin { val: 17, exp: 5 }, "1.7e5");
    }

    fn test_boundary_width(v: f64, b: f64, w: f64) {
        let bin: Bin = v.into();
        assert_eq!(b, Into::<f64>::into(&bin));
        assert_eq!(w, bin.width_signed());
    }

    #[test]
    fn boundary_width() {
        test_boundary_width(43.3, 43.0, 1.0);
        test_boundary_width(99.9, 99.0, 1.0);
        test_boundary_width(10.0, 10.0, 1.0);
        test_boundary_width(1.0, 1.0, 0.1);
        test_boundary_width(0.0002, 0.0002, 0.00001);
        test_boundary_width(0.003, 0.003, 0.0001);
        test_boundary_width(0.3301, 0.33, 0.01);
        test_boundary_width(0.0035, 0.0035, 0.0001);
        test_boundary_width(-1.0, -1.0, -0.1);
        test_boundary_width(-0.00123, -0.0012, -0.0001);
        test_boundary_width(-997324.0, -990000.0, -10000.0);
    }
    #[test]
    fn ops() {
        let b: Bin = 12345.0.into();
        assert_eq!(b.midpoint(), 12480.0);
        assert_eq!(b.left(), 12000.0);
        assert_eq!(b.width(), 1000.0);
        let b: Bin = (-12345.0).into();
        assert_eq!(b.midpoint(), -12480.0);
        assert_eq!(b.left(), -13000.0);
        assert_eq!(b.width(), 1000.0);
    }
    #[test]
    fn from_duration() {
        let s = Duration::from_secs(42);
        assert_eq!(bin![42, 1], Bin::from_duration(s));
        let ns = Duration::from_nanos(42);
        assert_eq!(bin![42, -8], Bin::from_duration(ns));
    }
    #[test]
    fn display() {
        let b: Bin = 4200.0.into();
        assert_eq!(b.to_string(), "4.2e3")
    }
}
