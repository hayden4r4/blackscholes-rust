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

use num_traits::{Float, FromPrimitive, NumCast};
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

pub(crate) const N_MEAN: f64 = 0.0;
pub(crate) const N_STD_DEV: f64 = 1.0;
pub(crate) const SQRT_2PI: f64 = 2.506628274630002;
pub(crate) const HALF: f64 = 0.5;
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
/// Tuple (f32, f32) of (d1, d2)
pub(crate) fn calc_d1d2<T>(inputs: &Inputs<T>) -> Result<(T, T), String>
where
    T: Float + FromPrimitive,
{
    let sigma = inputs
        .sigma
        .ok_or("Expected Some(float) for self.sigma, received None")?;
    // Calculating numerator of d1
    let part1 = (inputs.s / inputs.k).ln();

    if part1.is_infinite() {
        return Err("Log from s/k is infinity".to_string());
    }

    let part2 = (inputs.r - inputs.q + (sigma.powi(2)) * T::from(0.5).unwrap()) * inputs.t;
    let numd1 = part1 + part2;

    // Calculating denominator of d1 and d2
    if inputs.t == T::zero() {
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
/// Tuple (f32, f32) of (nd1, nd2)
pub(crate) fn calc_nd1nd2<T>(inputs: &Inputs<T>) -> Result<(T, T), String>
where
    T: Float + FromPrimitive,
{
    let nd1nd2 = {
        let d1d2 = calc_d1d2(inputs)?;

        let n: Normal = Normal::new(N_MEAN, N_STD_DEV).unwrap();

        let num_cast_err: String = "Failed to cast f64 to f32".into();
        // Calculates the nd1 and nd2 values
        // Checks if OptionType is Call or Put
        match inputs.option_type {
            OptionType::Call => (
                NumCast::from(n.cdf(NumCast::from(d1d2.0).ok_or(&num_cast_err)?))
                    .ok_or(&num_cast_err)?,
                NumCast::from(n.cdf(NumCast::from(d1d2.1).ok_or(&num_cast_err)?))
                    .ok_or(&num_cast_err)?,
            ),
            OptionType::Put => (
                NumCast::from(n.cdf(NumCast::from(-d1d2.0).ok_or(&num_cast_err)?))
                    .ok_or(&num_cast_err)?,
                NumCast::from(n.cdf(NumCast::from(-d1d2.1).ok_or(&num_cast_err)?))
                    .ok_or(&num_cast_err)?,
            ),
        }
    };
    Ok(nd1nd2)
}

/// Calculates the n probability density function (PDF) for the given input.
/// # Returns
/// f32 of the value of the n probability density function.
pub(crate) fn calc_npdf<T>(x: T) -> T
where
    T: Float + FromPrimitive,
{
    let n_mean = T::from(N_MEAN).unwrap();
    let n_std_dev = T::from(N_STD_DEV).unwrap();
    let sqrt_2pi = T::from(SQRT_2PI).unwrap();
    let half = T::from(HALF).unwrap();
    let d: T = (x - n_mean) / n_std_dev;
    (-half * d * d).exp() / (sqrt_2pi * n_std_dev)
}

/// # Returns
/// f32 of the derivative of the nd1.
pub fn calc_nprimed1<T>(inputs: &Inputs<T>) -> Result<T, String>
where
    T: Float + FromPrimitive,
{
    let (d1, _) = calc_d1d2(inputs)?;

    // Get the standard n probability density function value of d1
    let nprimed1 = calc_npdf(d1);
    Ok(nprimed1)
}

/// # Returns
/// f32 of the derivative of the nd2.
pub(crate) fn calc_nprimed2<T>(inputs: &Inputs<T>) -> Result<T, String>
where
    T: Float + FromPrimitive,
{
    let (_, d2) = calc_d1d2(inputs)?;

    // Get the standard n probability density function value of d1
    let nprimed2 = calc_npdf(d2);
    Ok(nprimed2)
}
