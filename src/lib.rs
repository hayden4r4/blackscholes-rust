//! This library provides an simple, lightweight, and efficient (though not heavily optimized) implementation of the Black-Scholes-Merton model for pricing European options.
//!
//! Provides methods for pricing options, calculating implied volatility, and calculating the first, second, and third order Greeks.
//!
//! ### Example:
//! ```
//! use blackscholes::{Inputs, OptionType, Pricing};
//! let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
//! let price: f64 = inputs.calc_price().unwrap();
//! ```
//!
//! Criterion benchmark can be ran by running:
//! ```bash
//! cargo bench
//! ```
//!
//! See the [Github Repo](https://github.com/hayden4r4/blackscholes-rust/tree/master) for full source code.  Other implementations such as a [npm WASM package](https://www.npmjs.com/package/@haydenr4/blackscholes_wasm) and a [python module](https://pypi.org/project/blackscholes/) are also available.

pub use greeks::Greeks;
pub use implied_volatility::ImpliedVolatility;
pub use inputs::{Inputs, OptionType};
use lets_be_rational::normal_distribution::{standard_normal_cdf, standard_normal_pdf};
pub use pricing::Pricing;

mod greeks;
mod implied_volatility;
mod inputs;
pub mod lets_be_rational;
mod pricing;

pub(crate) const DAYS_PER_YEAR: f64 = 365.25;

pub(crate) const A: f64 = 4.626_275_3e-1;
pub(crate) const B: f64 = -1.168_519_2e-2;
pub(crate) const C: f64 = 9.635_418_5e-4;
pub(crate) const D: f64 = 7.535_022_5e-5;
pub(crate) const _E: f64 = 1.424_516_45e-5;
pub(crate) const F: f64 = -2.102_376_9e-5;

/// Calculates the d1 and d2 values for the option.
/// # Requires
/// s, k, r, q, t, sigma.
/// # Returns
/// Tuple (f64, f64) of (d1, d2)
pub(crate) fn calc_d1d2(inputs: &Inputs) -> Result<(f64, f64), String> {
    let sigma = inputs
        .sigma
        .ok_or("Expected Some(f64) for self.sigma, received None")?;
    // Calculating numerator of d1
    let part1 = (inputs.s / inputs.k).ln();

    if part1.is_infinite() {
        return Err("Log from s/k is infinity".to_string());
    }

    let part2 = (inputs.r - inputs.q + (sigma.powi(2)) / 2.0) * inputs.t;
    let numd1 = part1 + part2;

    // Calculating denominator of d1 and d2
    if inputs.t == 0.0 {
        return Err("Time to maturity is 0".to_string());
    }

    let den = sigma * (inputs.t.sqrt());

    let d1 = numd1 / den;
    let d2 = d1 - den;

    Ok((d1, d2))
}

/// Calculates the nd1 and nd2 values for the option.
/// # Requires
/// s, k, r, q, t, sigma
/// # Returns
/// Tuple (f64, f64) of (nd1, nd2)
pub(crate) fn calc_nd1nd2(inputs: &Inputs) -> Result<(f64, f64), String> {
    let (d1, d2) = calc_d1d2(inputs)?;

    // Calculates the nd1 and nd2 values
    // Checks if OptionType is Call or Put
    match inputs.option_type {
        OptionType::Call => Ok((standard_normal_cdf(d1), standard_normal_cdf(d2))),
        OptionType::Put => Ok((standard_normal_cdf(-d1), standard_normal_cdf(-d2))),
    }
}

/// # Returns
/// f64 of the derivative of the nd1.
pub fn calc_nprimed1(inputs: &Inputs) -> Result<f64, String> {
    let (d1, _) = calc_d1d2(inputs)?;

    // Get the standard n probability density function value of d1
    let nprimed1 = standard_normal_pdf(d1);
    Ok(nprimed1)
}

/// # Returns
/// f64 of the derivative of the nd2.
pub(crate) fn calc_nprimed2(inputs: &Inputs) -> Result<f64, String> {
    let (_, d2) = calc_d1d2(inputs)?;

    // Get the standard n probability density function value of d1
    let nprimed2 = standard_normal_pdf(d2);
    Ok(nprimed2)
}
