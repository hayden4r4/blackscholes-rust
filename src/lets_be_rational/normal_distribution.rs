use statrs::distribution::{Continuous, ContinuousCDF, Normal};

pub fn standard_normal_cdf(x: f64) -> f64 {
    Normal::standard().cdf(x)
}

pub fn inverse_normal_cdf(x: f64) -> f64 {
    Normal::standard().inverse_cdf(x)
}

pub fn inverse_f_upper_map(f: f64) -> f64 {
    -2.0 * inverse_normal_cdf(f)
}

pub fn standard_normal_pdf(x: f64) -> f64 {
    Normal::standard().pdf(x)
}
