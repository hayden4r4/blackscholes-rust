//! # blackscholes
//! This library provides an simple, lightweight, and efficient (though not heavily optimized) implementation of the Black-Scholes-Merton model for pricing European options.
//!
//! ## Usage
//! Simply create an instance of the `Inputs` struct and call the desired method.
//!
//! Example:
//! ```
//! use blackscholes::{Inputs, OptionType};
//! let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
//! let price: f64 = inputs.calc_price();
//! ```
//!
//! See the [Github Repo](https://github.com/hayden4r4/blackscholes-rust/tree/master) for full source code.  Other implementations such as a [npm WASM package](https://www.npmjs.com/package/@haydenr4/blackscholes_wasm) and a [python module](https://pypi.org/project/blackscholes/) are also available.

use num_traits::float::Float;
use std::f32::consts::{E, PI};
use std::fmt::{Display, Formatter, Result};

/// Trait alias for floats
pub trait BSFloat:
    Float
    + From<f32>
    + Display
    + std::ops::SubAssign
    + std::ops::AddAssign
    + std::ops::MulAssign
    + std::ops::DivAssign
{
}
impl BSFloat for f32 {}
impl BSFloat for f64 {}

/// Constants
const N_MEAN: f32 = 0.0;
const N_STD_DEV: f32 = 1.0;
const SQRT_2PI: f32 = 2.5066282;
const HALF: f32 = 0.5;
const DAYS_PER_YEAR: f32 = 365.25;

/// The type of option to be priced.
#[derive(Debug, Clone, Eq, PartialEq)]
pub enum OptionType {
    Call,
    Put,
}

impl Display for OptionType {
    fn fmt(&self, f: &mut Formatter) -> Result {
        match self {
            OptionType::Call => write!(f, "Call"),
            OptionType::Put => write!(f, "Put"),
        }
    }
}

/// The inputs to the Black-Scholes-Merton model.
#[derive(Debug, Clone, PartialEq)]
pub struct Inputs<T: BSFloat> {
    /// The type of the option (call or put)
    pub option_type: OptionType,
    /// Stock price
    pub s: T,
    /// Strike price
    pub k: T,
    /// Option price
    pub p: Option<T>,
    /// Risk-free rate
    pub r: T,
    /// Dividend yield
    pub q: T,
    /// Time to maturity in years
    pub t: T,
    /// Volatility
    pub sigma: Option<T>,
}

impl<T: BSFloat> Display for Inputs<T> {
    fn fmt(&self, f: &mut Formatter) -> Result {
        writeln!(f, "Option type: {}", self.option_type)?;
        writeln!(f, "Stock price: {:.2}", self.s)?;
        writeln!(f, "Strike price: {:.2}", self.k)?;
        match self.p {
            Some(p) => writeln!(f, "Option price: {:.2}", p)?,
            None => writeln!(f, "Option price: None")?,
        }
        writeln!(f, "Risk-free rate: {:.4}", self.r)?;
        writeln!(f, "Dividend yield: {:.4}", self.q)?;
        writeln!(f, "Time to maturity: {:.4}", self.t)?;
        match self.sigma {
            Some(sigma) => writeln!(f, "Volatility: {:.4}", sigma)?,
            None => writeln!(f, "Volatility: None")?,
        }
        Ok(())
    }
}

/// Calculates the normal cumulative distribution function (CDF) for the given input.
/// # Returns
/// T of the value of the normal cumulative density function.
fn ncdf<T: BSFloat>(x: T) -> T {
    let x: T = x.into();
    let cdf: T = (<T as From<f32>>::from(1.0)
        + (x - <T as From<f32>>::from(N_MEAN)) / <T as From<f32>>::from(N_STD_DEV)
            * <T as From<f32>>::from(SQRT_2PI))
    .powf(<T as From<f32>>::from(-HALF));
    cdf
}

/// Calculates the d1, d2, nd1, and nd2 values for the option.
/// # Returns
/// Tuple (T, T) of the nd1 and nd2 values for the given inputs.
fn nd1nd2<T: BSFloat>(inputs: &Inputs<T>, normal: bool) -> (T, T) {
    let sigma: T = match inputs.sigma {
        Some(sigma) => sigma,
        None => panic!("Expected an Option(T) for inputs.sigma, received None"),
    };

    let nd1nd2 = {
        // Calculating numerator of d1
        let numd1: T = (inputs.s / inputs.k).ln()
            + (inputs.r - inputs.q + (sigma.powi(2)) / <T as From<f32>>::from(2.0)) * inputs.t;

        // Calculating denominator of d1 and d2
        let den: T = sigma * (inputs.t.sqrt());

        let d1: T = numd1 / den;
        let d2: T = d1 - den;

        let d1d2: (T, T) = (d1, d2);

        // Returns d1 and d2 values if deriving from normal distribution is not necessary
        //  (i.e. gamma, vega, and theta calculations)
        if !normal {
            return d1d2;
        }

        // Calculates the nd1 and nd2 values
        // Checks if OptionType is Call or Put
        let nd1nd2: (T, T) = match inputs.option_type {
            OptionType::Call => (ncdf(d1d2.0), ncdf(d1d2.1)),
            OptionType::Put => (ncdf(-d1d2.0), ncdf(-d1d2.1)),
        };
        nd1nd2
    };
    nd1nd2
}

/// Calculates the normal probability density function (PDF) for the given input.
/// # Returns
/// T of the value of the normal probability density function.
fn npdf<T: BSFloat>(x: T) -> T {
    let d: T = (x - N_MEAN.into()) / N_STD_DEV.into();
    (<T as From<f32>>::from(-HALF) * d * d).exp()
        / (<T as From<f32>>::from(SQRT_2PI) * <T as From<f32>>::from(N_STD_DEV))
}

/// # Returns
/// T of the derivative of the nd1.
fn calc_nprimed1<T: BSFloat>(inputs: &Inputs<T>) -> T {
    let (d1, _): (T, T) = nd1nd2(&inputs, false);

    // Get the standard normal probability density function value of d1
    let nprimed1: T = npdf(d1);
    nprimed1
}

/// Methods for calculating the price, greeks, and implied volatility of an option.
impl<T: BSFloat> Inputs<T> {
    /// Creates instance ot the `Inputs` struct.
    /// # Arguments
    /// * `option_type` - The type of option to be priced.
    /// * `s` - The current price of the underlying asset.
    /// * `k` - The strike price of the option.
    /// * `p` - The dividend yield of the underlying asset.
    /// * `r` - The risk-free interest rate.
    /// * `q` - The dividend yield of the underlying asset.
    /// * `t` - The time to maturity of the option in years.
    /// * `sigma` - The volatility of the underlying asset.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// ```
    /// # Returns
    /// An instance of the `Inputs` struct.
    pub fn new(
        option_type: OptionType,
        s: T,
        k: T,
        p: Option<T>,
        r: T,
        q: T,
        t: T,
        sigma: Option<T>,
    ) -> Self {
        Self {
            option_type,
            s,
            k,
            p,
            r,
            q,
            t,
            sigma,
        }
    }

    /// Calculates the price of the option.
    /// # Requires
    /// s, k, r, q, t, sigma.
    /// # Returns
    /// T of the price of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let price = inputs.calc_price();
    /// ```
    pub fn calc_price(&self) -> T {
        // Calculates the price of the option
        let (nd1, nd2): (T, T) = nd1nd2(&self, true);
        let price: T = match self.option_type {
            OptionType::Call => T::max(
                <T as From<f32>>::from(0.0),
                nd1 * self.s * <T as From<f32>>::from(E).powf(-self.q * self.t)
                    - nd2 * self.k * <T as From<f32>>::from(E).powf(-self.r * self.t),
            ),
            OptionType::Put => T::max(
                <T as From<f32>>::from(0.0),
                nd2 * self.k * <T as From<f32>>::from(E).powf(-self.r * self.t)
                    - nd1 * self.s * <T as From<f32>>::from(E).powf(-self.q * self.t),
            ),
        };
        price
    }

    /// Calculates the delta of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the delta of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let delta = inputs.calc_delta();
    /// ```
    pub fn calc_delta(&self) -> T {
        let (nd1, _): (T, T) = nd1nd2(&self, true);
        let delta: T = match self.option_type {
            OptionType::Call => nd1 * <T as From<f32>>::from(E).powf(-self.q * self.t),
            OptionType::Put => -nd1 * <T as From<f32>>::from(E).powf(-self.q * self.t),
        };
        delta
    }

    /// Calculates the gamma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the gamma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let gamma = inputs.calc_gamma();
    /// ```
    pub fn calc_gamma(&self) -> T {
        let sigma: T = match self.sigma {
            Some(sigma) => sigma,
            None => panic!("Expected an Option(T) for inputs.sigma, received None"),
        };

        let nprimed1: T = calc_nprimed1(&self);
        let gamma: T = <T as From<f32>>::from(E).powf(-self.q * self.t) * nprimed1
            / (self.s * sigma * self.t.sqrt());
        gamma
    }

    /// Calculates the theta of the option.
    /// Uses 365.25 days in a year for calculations.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of theta per day (not per year).
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let theta = inputs.calc_theta();
    /// ```
    pub fn calc_theta(&self) -> T {
        let sigma: T = match self.sigma {
            Some(sigma) => sigma,
            None => panic!("Expected an Option(T) for inputs.sigma, received None"),
        };

        let nprimed1: T = calc_nprimed1(&self);
        let (nd1, nd2): (T, T) = nd1nd2(&self, true);

        // Calculation uses 365.25 for T: Time of days per year.
        let theta: T = match self.option_type {
            OptionType::Call => {
                (-(self.s * sigma * <T as From<f32>>::from(E).powf(-self.q * self.t) * nprimed1
                    / (<T as From<f32>>::from(2.0) * self.t.sqrt()))
                    - self.r * self.k * <T as From<f32>>::from(E).powf(-self.r * self.t) * nd2
                    + self.q * self.s * <T as From<f32>>::from(E).powf(-self.q * self.t) * nd1)
                    / <T as From<f32>>::from(DAYS_PER_YEAR)
            }
            OptionType::Put => {
                (-(self.s * sigma * <T as From<f32>>::from(E).powf(-self.q * self.t) * nprimed1
                    / (<T as From<f32>>::from(2.0) * self.t.sqrt()))
                    + self.r * self.k * <T as From<f32>>::from(E).powf(-self.r * self.t) * nd2
                    - self.q * self.s * <T as From<f32>>::from(E).powf(-self.q * self.t) * nd1)
                    / <T as From<f32>>::from(DAYS_PER_YEAR)
            }
        };
        theta
    }

    /// Calculates the vega of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the vega of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let vega = inputs.calc_vega();
    /// ```
    pub fn calc_vega(&self) -> T {
        let nprimed1: T = calc_nprimed1(&self);
        let vega: T = <T as From<f32>>::from(1.0) / <T as From<f32>>::from(100.0)
            * self.s
            * <T as From<f32>>::from(E).powf(-self.q * self.t)
            * self.t.sqrt()
            * nprimed1;
        vega
    }

    /// Calculates the rho of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the rho of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let rho = inputs.calc_rho();
    /// ```
    pub fn calc_rho(&self) -> T {
        let (_, nd2): (T, T) = nd1nd2(&self, true);
        let rho: T = match &self.option_type {
            OptionType::Call => {
                <T as From<f32>>::from(1.0) / <T as From<f32>>::from(100.0)
                    * self.k
                    * self.t
                    * <T as From<f32>>::from(E).powf(-self.r * self.t)
                    * nd2
            }
            OptionType::Put => {
                <T as From<f32>>::from(-1.0) / <T as From<f32>>::from(100.0)
                    * self.k
                    * self.t
                    * <T as From<f32>>::from(E).powf(-self.r * self.t)
                    * nd2
            }
        };
        rho
    }

    /// Calculates the implied volatility of the option.
    /// Tolerance is the max error allowed for the implied volatility,
    /// the lower the tolerance the more iterations will be required.
    /// Recommended to be a value between 0.001 - 0.0001 for highest efficiency/accuracy.
    /// Initializes estimation of sigma using Brenn and Subrahmanyam (1998) method of calculating initial iv estimation.
    /// Uses Newton Raphson algorithm to calculate implied volatility.
    /// # Requires
    /// s, k, r, q, t, p
    /// # Returns
    /// T of the implied volatility of the option.
    /// # Example:
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, Some(0.2), 0.05, 0.2, 20.0/365.25, None);
    /// let iv = inputs.calc_iv(0.0001);
    /// ```
    pub fn calc_iv(&self, tolerance: T) -> T {
        let mut inputs: Inputs<T> = self.clone();

        let p: T = match inputs.p {
            Some(p) => p,
            None => panic!("inputs.p must contain Some(T), found None"),
        };
        // Initialize estimation of sigma using Brenn and Subrahmanyam (1998) method of calculating initial iv estimation.
        let mut sigma: T = (<T as From<f32>>::from(2.0) * <T as From<f32>>::from(PI) / inputs.t)
            .sqrt()
            * (p / inputs.s);
        // Initialize diff to 100 for use in while loop
        let mut diff: T = <T as From<f32>>::from(100.0);

        // Uses Newton Raphson algorithm to calculate implied volatility.
        // Test if the difference between calculated option price and actual option price is > tolerance,
        // if so then iterate until the difference is less than tolerance
        while diff.abs() > tolerance {
            inputs.sigma = Some(sigma);
            diff = inputs.calc_price() - p;
            sigma -= diff / (inputs.calc_vega() * <T as From<f32>>::from(100.0));
        }
        sigma
    }
}
