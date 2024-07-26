use statrs::{consts::SQRT_2PI, function::erf};
use std::f64::consts::SQRT_2;

const MEAN: f64 = 0.0;
const STD_DEV: f64 = 1.0;

// Re-implementations of statrs::distribution::Normal::standard()
// with const values for compiler optimization
pub fn standard_normal_cdf(x: f64) -> f64 {
    0.5 * erf::erfc((MEAN - x) / (STD_DEV * SQRT_2))
}

pub fn inverse_normal_cdf(x: f64) -> f64 {
    if !(0.0..=1.0).contains(&x) {
        panic!("x must be in [0, 1]");
    } else {
        MEAN - (STD_DEV * SQRT_2 * erf::erfc_inv(2.0 * x))
    }
}

pub fn inverse_f_upper_map(f: f64) -> f64 {
    -2.0 * inverse_normal_cdf(f)
}

pub fn standard_normal_pdf(x: f64) -> f64 {
    let d = (x - MEAN) / STD_DEV;
    (-0.5 * d * d).exp() / (SQRT_2PI * STD_DEV)
}
