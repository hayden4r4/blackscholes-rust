use crate::lets_be_rational::black::normalised_black;
use crate::OptionType;
use statrs::consts::SQRT_2PI;

mod black;

#[cfg(any(test, feature = "testing"))]
pub mod test_utils {
    pub use super::black::*;
}

mod intrinsic;
mod normal_distribution;
mod rational_cubic;
mod so_rational;

const IMPLIED_VOLATILITY_MAXIMUM_ITERATIONS: i32 = 2;
pub(crate) const DENORMALISATION_CUTOFF: f64 = 0.0;
pub(crate) const ONE_OVER_SQRT_TWO_PI: f64 = 1.0 / SQRT_2PI;

/// Calculates the price of a European option using the Black model.
///
/// This function computes the theoretical price of a European call or put option based on the Black model, which is an extension of the Black-Scholes model for futures contracts. The model assumes that the price of the underlying asset follows a geometric Brownian motion and that markets are frictionless.
///
/// # Arguments
/// * `forward_price` - The forward price of the underlying asset.
/// * `strike_price` - The strike price of the option.
/// * `sigma` - The volatility of the underlying asset's returns.
/// * `time_to_maturity` - The time to maturity of the option, in years.
/// * `option_type` - The type of the option (call or put), represented by `OptionType`.
///
/// # Returns
/// The theoretical price of the option as a `f64`.
///
/// # Examples
/// ```
/// use blackscholes::{OptionType, lets_be_rational::black};
///
/// let forward_price = 100.0;
/// let strike_price = 95.0;
/// let sigma = 0.2;
/// let time_to_maturity = 1.0;
/// let option_type = OptionType::Call; // For a call option
///
/// let price = black(forward_price, strike_price, sigma, time_to_maturity, option_type);
/// println!("The price of the option is: {}", price);
/// ```
///
/// # Note
/// The function uses the natural logarithm of the forward price over the strike price,
/// multiplies it by the square root of time to maturity, and applies the option type
/// to determine the final price. It's suitable for European options *only*.
pub fn black(
    forward_price: f64,
    strike_price: f64,
    sigma: f64,
    time_to_maturity: f64,
    option_type: OptionType,
) -> f64 {
    let q: f64 = option_type.into();
    let intrinsic = (if q < 0.0 {
        strike_price - forward_price
    } else {
        forward_price - strike_price
    })
    .max(0.0)
    .abs();
    // Map in-the-money to out-of-the-money
    if q * (forward_price - strike_price) > 0.0 {
        return intrinsic
            + black(
                forward_price,
                strike_price,
                sigma,
                time_to_maturity,
                -option_type,
            );
    }
    intrinsic.max(
        (forward_price.sqrt() * strike_price.sqrt())
            * normalised_black(
                (forward_price / strike_price).ln(),
                sigma * time_to_maturity.sqrt(),
                q,
            ),
    )
}

pub fn implied_volatility_from_a_transformed_rational_guess(
    market_price: f64,
    forward_price: f64,
    strike_price: f64,
    time_to_maturity: f64,
    option_type: OptionType,
) -> f64 {
    so_rational::implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(
        market_price,
        forward_price,
        strike_price,
        time_to_maturity,
        option_type,
        IMPLIED_VOLATILITY_MAXIMUM_ITERATIONS,
    )
}
