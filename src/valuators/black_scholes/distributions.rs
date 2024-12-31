use crate::lets_be_rational::normal_distribution::{standard_normal_cdf, standard_normal_pdf};
use crate::valuators::black_scholes::Inputs;
use crate::OptionType;

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
