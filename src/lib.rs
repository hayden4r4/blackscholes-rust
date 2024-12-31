//! This library provides an simple, lightweight, and efficient (though not heavily optimized) implementation of various option valuation models for pricing options.
//!
//! Provides methods for pricing options, calculating implied volatility, and calculating the first, second, and third order Greeks.
//!
//! ### Example:
//! ```
//! use option_valuators::{OptionType, Pricing};
//! use option_valuators::valuators::black_scholes::Inputs;
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

pub(crate) mod constants;
pub mod lets_be_rational;
mod option_type;
pub use option_type::OptionType;
mod implied_volatility;
pub mod valuators;
pub use implied_volatility::ImpliedVolatility;
mod pricing;
pub use pricing::Pricing;
mod greeks;
pub use greeks::Greeks;
mod shift;
pub use shift::Shift;
