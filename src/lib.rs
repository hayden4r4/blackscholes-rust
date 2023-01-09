//! # blackscholes
//! This library provides an simple, lightweight, and efficient (though not heavily optimized) implementation of the Black-Scholes-Merton model for pricing European options.
//!
//! ## Usage
//! Simply create an instance of the `Inputs` struct and call the desired method.
//! This crate is generic over all floating point types (f32, f64)
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
use num_traits::NumCast;
use statrs::distribution::{ContinuousCDF, Normal};
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

/// Calculates the d1, d2, nd1, and nd2 values for the option.
/// # Requires
/// * inputs: &Inputs<T> - The inputs to the Black-Scholes-Merton model
/// * n: bool - Whether or not to calculate the nd1 and nd2 values, if false, only d1 and d2 are calculated and returned
/// # Returns
/// Tuple (T, T) of either (nd1, nd2) if n==True or (d1, d2) if n==false
fn calc_nd1nd2<T: BSFloat>(inputs: &Inputs<T>, n: bool) -> (T, T) {
    let sigma: T = match inputs.sigma {
        Some(sigma) => sigma,
        None => panic!("Expected Some(T) for inputs.sigma, received None"),
    };

    let nd1nd2 = {
        // Calculating numerator of d1
        let numd1 = (inputs.s / inputs.k).ln()
            + (inputs.r - inputs.q + (sigma.powi(2)) / <T as From<f32>>::from(2.0)) * inputs.t;

        // Calculating denominator of d1 and d2
        let den = sigma * (inputs.t.sqrt());

        let d1 = numd1 / den;
        let d2 = d1 - den;

        let d1d2 = (d1, d2);

        // Returns d1 and d2 values if deriving from n distribution is not necessary
        //  (i.e. gamma, vega, and theta calculations)
        if !n {
            return d1d2;
        }

        let n: Normal = Normal::new(N_MEAN as f64, N_STD_DEV as f64).unwrap();

        // Calculates the nd1 and nd2 values
        // Checks if OptionType is Call or Put
        let nd1nd2 = match inputs.option_type {
            OptionType::Call => (
                NumCast::from(n.cdf(NumCast::from(d1d2.0).unwrap())).unwrap(),
                NumCast::from(n.cdf(NumCast::from(d1d2.1).unwrap())).unwrap(),
            ),
            OptionType::Put => (
                NumCast::from(n.cdf(NumCast::from(-d1d2.0).unwrap())).unwrap(),
                NumCast::from(n.cdf(NumCast::from(-d1d2.1).unwrap())).unwrap(),
            ),
        };
        nd1nd2
    };
    nd1nd2
}

/// Calculates the n probability density function (PDF) for the given input.
/// # Returns
/// T of the value of the n probability density function.
fn calc_npdf<T: BSFloat>(x: T) -> T {
    let d: T = (x - N_MEAN.into()) / N_STD_DEV.into();
    (<T as From<f32>>::from(-HALF) * d * d).exp()
        / (<T as From<f32>>::from(SQRT_2PI) * <T as From<f32>>::from(N_STD_DEV))
}

/// # Returns
/// T of the derivative of the nd1.
fn calc_nprimed1<T: BSFloat>(inputs: &Inputs<T>) -> T {
    let (d1, _) = calc_nd1nd2(&inputs, false);

    // Get the standard n probability density function value of d1
    let nprimed1 = calc_npdf(d1);
    nprimed1
}

/// # Returns
/// T of the derivative of the nd2.
fn calc_nprimed2<T: BSFloat>(inputs: &Inputs<T>) -> T {
    let (_, d2) = calc_nd1nd2(&inputs, false);

    // Get the standard n probability density function value of d1
    let nprimed2 = calc_npdf(d2);
    nprimed2
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
        let (nd1, nd2): (T, T) = calc_nd1nd2(&self, true);
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
        let (nd1, _): (T, T) = calc_nd1nd2(&self, true);
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
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            panic!("Expected Some(T) for inputs.sigma, received None")
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
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            panic!("Expected Some(T) for inputs.sigma, received None")
        };

        let nprimed1: T = calc_nprimed1(&self);
        let (nd1, nd2): (T, T) = calc_nd1nd2(&self, true);

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
        let vega: T = <T as From<f32>>::from(0.01)
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
        let (_, nd2): (T, T) = calc_nd1nd2(&self, true);
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

    // The formulas for the greeks below are from the wikipedia page for the Black-Scholes greeks
    // https://en.wikipedia.org/wiki/Greeks_(finance)#Black.E2.80.93Scholes_Greeks
    // Some sources I reviewed contain variations of these formulas and/or varying values, therefore the
    // values returned by this library may not match other libraries or sources.
    // These functions have not been throughouly tested.

    /// Calculates the epsilon of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the epsilon of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let epsilon = inputs.calc_epsilon();
    /// ```
    pub fn calc_epsilon(&self) -> T {
        let (nd1, _) = calc_nd1nd2(&self, true);
        let e_negqt = <T as From<f32>>::from(E).powf(-self.q * self.t);
        let epsilon: T = match &self.option_type {
            OptionType::Call => -self.s * self.t * e_negqt * nd1,
            OptionType::Put => self.s * self.t * e_negqt * nd1,
        };
        epsilon
    }

    /// Calculates the vanna of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the vanna of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let vanna = inputs.calc_vanna();
    /// ```
    pub fn calc_vanna(&self) -> T {
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            panic!("Expected Some(T) for inputs.sigma, received None")
        };

        let nprimed1 = calc_nprimed1(&self);
        let (_, d2) = calc_nd1nd2(&self, false);
        let vanna: T = d2
            * <T as From<f32>>::from(E).powf(-self.q * self.t)
            * nprimed1
            * <T as From<f32>>::from(-0.01)
            / sigma;
        vanna
    }

    // /// Calculates the charm of the option.
    // /// # Requires
    // /// s, k, r, q, t, sigma
    // /// # Returns
    // /// T of the charm of the option.
    // /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let charm = inputs.calc_charm();
    /// ```
    pub fn calc_charm(&self) -> T {
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            panic!("Expected Some(T) for inputs.sigma, received None")
        };

        let nprimed1 = calc_nprimed1(&self);
        let (nd1, _) = calc_nd1nd2(&self, true);
        let (_, d2) = calc_nd1nd2(&self, false);
        let e_negqt = <T as From<f32>>::from(E).powf(-self.q * self.t);

        let charm = match &self.option_type {
            OptionType::Call => {
                self.q * e_negqt * nd1
                    - e_negqt
                        * nprimed1
                        * (<T as From<f32>>::from(2.0) * (self.r - self.q) * self.t
                            - d2 * sigma * self.t.sqrt())
                        / (<T as From<f32>>::from(2.0) * self.t * sigma * self.t.sqrt())
            }
            OptionType::Put => {
                -self.q * e_negqt * nd1
                    - e_negqt
                        * nprimed1
                        * (<T as From<f32>>::from(2.0) * (self.r - self.q) * self.t
                            - d2 * sigma * self.t.sqrt())
                        / (<T as From<f32>>::from(2.0) * self.t * sigma * self.t.sqrt())
            }
        };
        charm
    }

    /// Calculates the veta of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the veta of the option.
    pub fn calc_veta(&self) -> T {
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            panic!("Expected Some(T) for inputs.sigma, received None")
        };

        let nprimed1 = calc_nprimed1(&self);
        let (d1, d2) = calc_nd1nd2(&self, false);
        let e_negqt = <T as From<f32>>::from(E).powf(-self.q * self.t);

        let veta = -self.s
            * e_negqt
            * nprimed1
            * self.t.sqrt()
            * (self.q + ((self.r - self.q) * d1) / (sigma * self.t.sqrt())
                - ((<T as From<f32>>::from(1.0) + d1 * d2)
                    / (<T as From<f32>>::from(2.0) * self.t)));
        veta
    }

    /// Calculates the vomma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the vomma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let vomma = inputs.calc_vomma();
    /// ```
    pub fn calc_vomma(&self) -> T {
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            panic!("Expected Some(t) for inputs.sigma, received None")
        };

        let (d1, d2) = calc_nd1nd2(&self, false);

        let vomma = self.calc_vega() * ((d1 * d2) / sigma);
        vomma
    }

    /// Calculates the speed of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the speed of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let speed = inputs.calc_speed();
    /// ```
    pub fn calc_speed(&self) -> T {
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            panic!("Expected Some(t) for inputs.sigma, received None")
        };

        let (d1, _) = calc_nd1nd2(&self, false);
        let gamma = self.calc_gamma();

        let speed = -gamma / self.s * (d1 / (sigma * self.t.sqrt()) + <T as From<f32>>::from(1.0));
        speed
    }

    /// Calculates the zomma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the zomma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let zomma = inputs.calc_zomma();
    /// ```
    pub fn calc_zomma(&self) -> T {
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            panic!("Expected Some(t) for inputs.sigma, received None")
        };

        let (d1, d2) = calc_nd1nd2(&self, false);
        let gamma = self.calc_gamma();

        let zomma = gamma * ((d1 * d2 - <T as From<f32>>::from(1.0)) / sigma);
        zomma
    }

    /// Calculates the color of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the color of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let color = inputs.calc_color();
    /// ```
    pub fn calc_color(&self) -> T {
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            panic!("Expected Some(t) for inputs.sigma, received None")
        };

        let (d1, d2) = calc_nd1nd2(&self, false);
        let nprimed1 = calc_nprimed1(&self);
        let e_negqt = <T as From<f32>>::from(E).powf(-self.q * self.t);

        let color = -e_negqt
            * (nprimed1 / (<T as From<f32>>::from(2.0) * self.s * self.t * sigma * self.t.sqrt()))
            * (<T as From<f32>>::from(2.0) * self.q * self.t
                + <T as From<f32>>::from(1.0)
                + (<T as From<f32>>::from(2.0) * (self.r - self.q) * self.t
                    - d2 * sigma * self.t.sqrt())
                    / (sigma * self.t.sqrt())
                    * d1);
        color
    }

    /// Calculates the ultima of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the ultima of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let ultima = inputs.calc_ultima();
    /// ```
    pub fn calc_ultima(&self) -> T {
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            panic!("Expected Some(t) for inputs.sigma, received None")
        };

        let (d1, d2) = calc_nd1nd2(&self, false);
        let vega = self.calc_vega();

        let ultima = -vega / sigma.powf(<T as From<f32>>::from(2.0))
            * (d1 * d2 * (<T as From<f32>>::from(1.0) - d1 * d2)
                + d1.powf(<T as From<f32>>::from(2.0))
                + d2.powf(<T as From<f32>>::from(2.0)));
        ultima
    }

    /// Calculates the dual delta of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the dual delta of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let dual_delta = inputs.calc_dual_delta();
    /// ```
    pub fn calc_dual_delta(&self) -> T {
        let (_, nd2) = calc_nd1nd2(&self, true);
        let e_negqt = <T as From<f32>>::from(E).powf(-self.q * self.t);

        let dual_delta = match self.option_type {
            OptionType::Call => -e_negqt * nd2,
            OptionType::Put => e_negqt * nd2,
        };
        dual_delta
    }

    /// Calculates the dual gamma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the dual gamma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let dual_gamma = inputs.calc_dual_gamma();
    /// ```
    pub fn calc_dual_gamma(&self) -> T {
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            panic!("Expected Some(t) for inputs.sigma, received None")
        };

        let nprimed2 = calc_nprimed2(&self);
        let e_negqt = <T as From<f32>>::from(E).powf(-self.q * self.t);

        let dual_gamma = e_negqt * (nprimed2 / (self.k * sigma * self.t.sqrt()));
        dual_gamma
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

        let p = if let Some(p) = inputs.p {
            p
        } else {
            panic!("inputs.p must contain Some(T), found None")
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
