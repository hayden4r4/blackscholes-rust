use once_cell::sync::Lazy;
use statrs::distribution::{Continuous, ContinuousCDF, Normal};

static STANDARD_NORMAL: Lazy<Normal> = Lazy::new(|| Normal::new(0.0, 1.0).unwrap());

pub fn standard_normal_pdf(x: f64) -> f64 {
    STANDARD_NORMAL.pdf(x)
}

pub fn standard_normal_cdf(x: f64) -> f64 {
    STANDARD_NORMAL.cdf(x)
}

pub fn standard_normal_inverse_cdf(p: f64) -> f64 {
    STANDARD_NORMAL.inverse_cdf(p)
}

pub fn inverse_norm_cdf(p: f64) -> f64 {
    STANDARD_NORMAL.inverse_cdf(p)
}

pub fn inverse_f_upper_map(f: f64) -> f64 {
    -2.0 * inverse_norm_cdf(f)
}
