//! This library provides an simple, lightweight, and efficient (though not heavily optimized) implementation of the Black-Scholes-Merton model for pricing European options.
//!
//! Provides methods for pricing options, calculating implied volatility, and calculating the first, second, and third order Greeks.
//!
//! ### Example:
//! ```
//! use blackscholes::{Inputs, OptionType, Pricing};
//! let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
//! let price: f32 = inputs.calc_price().unwrap();
//! ```
//!
//! Criterion benchmark can be ran by running:
//! ```bash
//! cargo bench
//! ```
//!
//! See the [Github Repo](https://github.com/hayden4r4/blackscholes-rust/tree/master) for full source code.  Other implementations such as a [npm WASM package](https://www.npmjs.com/package/@haydenr4/blackscholes_wasm) and a [python module](https://pypi.org/project/blackscholes/) are also available.

pub use std::f32::consts::{E, PI};

use num_traits::NumCast;
use statrs::distribution::{ContinuousCDF, Normal};

pub use greeks::Greeks;
pub use implied_volatility::ImpliedVolatility;
pub use inputs::{Inputs, OptionType};
pub use pricing::Pricing;

mod greeks;
mod implied_volatility;
mod inputs;
pub mod lets_be_rational;
mod pricing;

// Static error messages to avoid string allocations
static ERR_MISSING_SIGMA: &str = "Expected Some(f32) for self.sigma, received None";
static ERR_INFINITY: &str = "Log from s/k is infinity";
static ERR_ZERO_MATURITY: &str = "Time to maturity is 0";
static ERR_NUM_CAST: &str = "Failed to cast f64 to f32";

pub(crate) const N_MEAN: f32 = 0.0;
pub(crate) const N_STD_DEV: f32 = 1.0;
pub(crate) const SQRT_2PI: f32 = 2.5066282;
pub(crate) const HALF: f32 = 0.5;
pub(crate) const DAYS_PER_YEAR: f32 = 365.25;

pub(crate) const A: f32 = 4.626_275_3e-1;
pub(crate) const B: f32 = -1.168_519_2e-2;
pub(crate) const C: f32 = 9.635_418_5e-4;
pub(crate) const D: f32 = 7.535_022_5e-5;
pub(crate) const _E: f32 = 1.424_516_45e-5;
pub(crate) const F: f32 = -2.102_376_9e-5;

/// Calculates the d1 and d2 values for the option.
/// # Requires
/// s, k, r, q, t, sigma.
/// # Returns
/// Tuple (f32, f32) of (d1, d2)
#[inline]
pub(crate) fn calc_d1d2(inputs: &Inputs) -> Result<(f32, f32), String> {
    let sigma = inputs
        .sigma
        .ok_or(ERR_MISSING_SIGMA)?;
    // Calculating numerator of d1
    let part1 = (inputs.s / inputs.k).ln();

    if part1.is_infinite() {
        return Err(ERR_INFINITY.to_string());
    }

    let part2 = (inputs.r - inputs.q + (sigma.powi(2)) / 2.0) * inputs.t;
    let numd1 = part1 + part2;

    // Calculating denominator of d1 and d2
    if inputs.t == 0.0 {
        return Err(ERR_ZERO_MATURITY.to_string());
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
/// Tuple (f32, f32) of (nd1, nd2)
#[inline]
pub(crate) fn calc_nd1nd2(inputs: &Inputs) -> Result<(f32, f32), String> {
    let nd1nd2 = {
        let d1d2 = calc_d1d2(inputs)?;

        let n: Normal = Normal::new(N_MEAN as f64, N_STD_DEV as f64).unwrap();

        // Calculates the nd1 and nd2 values
        // Checks if OptionType is Call or Put
        match inputs.option_type {
            OptionType::Call => (
                NumCast::from(n.cdf(NumCast::from(d1d2.0).ok_or(ERR_NUM_CAST)?))
                    .ok_or(ERR_NUM_CAST)?,
                NumCast::from(n.cdf(NumCast::from(d1d2.1).ok_or(ERR_NUM_CAST)?))
                    .ok_or(ERR_NUM_CAST)?,
            ),
            OptionType::Put => (
                NumCast::from(n.cdf(NumCast::from(-d1d2.0).ok_or(ERR_NUM_CAST)?))
                    .ok_or(ERR_NUM_CAST)?,
                NumCast::from(n.cdf(NumCast::from(-d1d2.1).ok_or(ERR_NUM_CAST)?))
                    .ok_or(ERR_NUM_CAST)?,
            ),
        }
    };
    Ok(nd1nd2)
}

/// Calculates the n probability density function (PDF) for the given input.
/// # Returns
/// f32 of the value of the n probability density function.
#[inline]
pub(crate) fn calc_npdf(x: f32) -> f32 {
    let d: f32 = (x - N_MEAN) / N_STD_DEV;
    (-HALF * d * d).exp() / (SQRT_2PI * N_STD_DEV)
}

/// # Returns
/// f32 of the derivative of the nd1.
#[inline]
pub fn calc_nprimed1(inputs: &Inputs) -> Result<f32, String> {
    let (d1, _) = calc_d1d2(inputs)?;

    // Get the standard n probability density function value of d1
    let nprimed1 = calc_npdf(d1);
    Ok(nprimed1)
}

/// # Returns
/// f32 of the derivative of the nd2.
#[inline]
pub(crate) fn calc_nprimed2(inputs: &Inputs) -> Result<f32, String> {
    let (_, d2) = calc_d1d2(inputs)?;

    // Get the standard n probability density function value of d1
    let nprimed2 = calc_npdf(d2);
    Ok(nprimed2)
}
