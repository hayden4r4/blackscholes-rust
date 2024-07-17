use num_traits::{AsPrimitive, Float, FromPrimitive};
use once_cell::sync::Lazy;
use statrs::distribution::{ContinuousCDF, Normal};

static STANDARD_NORMAL: Lazy<Normal> = Lazy::new(|| Normal::new(0.0, 1.0).unwrap());

pub fn standard_normal_cdf<T>(x: T) -> T
where
    T: Float + FromPrimitive + AsPrimitive<f64>,
{
    T::from(STANDARD_NORMAL.cdf(x.as_())).unwrap()
}

pub fn inverse_normal_cdf<T>(p: T) -> T
where
    T: Float + FromPrimitive + AsPrimitive<f64>,
{
    T::from(STANDARD_NORMAL.inverse_cdf(p.as_())).unwrap()
}

pub fn inverse_f_upper_map<T>(f: T) -> T
where
    T: Float + FromPrimitive + AsPrimitive<f64>,
{
    T::from(-2.0 * inverse_normal_cdf(f.as_())).unwrap()
}
