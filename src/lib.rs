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
//! let price: f32 = inputs.calc_price().unwrap();
//! ```
//!
//! See the [Github Repo](https://github.com/hayden4r4/blackscholes-rust/tree/master) for full source code.  Other implementations such as a [npm WASM package](https://www.npmjs.com/package/@haydenr4/blackscholes_wasm) and a [python module](https://pypi.org/project/blackscholes/) are also available.

use num_traits::NumCast;
use statrs::distribution::{ContinuousCDF, Normal};
use std::f32::consts::{E, PI};
use std::fmt::{Display, Formatter, Result as fmtResult};

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
    fn fmt(&self, f: &mut Formatter) -> fmtResult {
        match self {
            OptionType::Call => write!(f, "Call"),
            OptionType::Put => write!(f, "Put"),
        }
    }
}

/// The inputs to the Black-Scholes-Merton model.
#[derive(Debug, Clone, PartialEq)]
pub struct Inputs {
    /// The type of the option (call or put)
    pub option_type: OptionType,
    /// Stock price
    pub s: f32,
    /// Strike price
    pub k: f32,
    /// Option price
    pub p: Option<f32>,
    /// Risk-free rate
    pub r: f32,
    /// Dividend yield
    pub q: f32,
    /// Time to maturity in years
    pub t: f32,
    /// Volatility
    pub sigma: Option<f32>,
}

impl Display for Inputs {
    fn fmt(&self, f: &mut Formatter) -> fmtResult {
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

/// Calculates the d1 and d2 values for the option.
/// # Requires
/// s, k, r, q, t, sigma.
/// # Returns
/// Tuple (f32, f32) of (d1, d2)
fn calc_d1d2(inputs: &Inputs) -> Result<(f32, f32), String> {
    let sigma = if let Some(sigma) = inputs.sigma {
        sigma
    } else {
        return Err("Expected Some(f32) for inputs.sigma, received None".into());
    };
    // Calculating numerator of d1
    let numd1 =
        (inputs.s / inputs.k).ln() + (inputs.r - inputs.q + (sigma.powi(2)) / 2.0) * inputs.t;

    // Calculating denominator of d1 and d2
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
fn calc_nd1nd2(inputs: &Inputs) -> Result<(f32, f32), String> {
    let nd1nd2 = {
        let d1d2 = calc_d1d2(inputs)?;

        let n: Normal = Normal::new(N_MEAN as f64, N_STD_DEV as f64).unwrap();

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
fn calc_npdf(x: f32) -> f32 {
    let d: f32 = (x - N_MEAN) / N_STD_DEV;
    (-HALF * d * d).exp() / (SQRT_2PI * N_STD_DEV)
}

/// # Returns
/// f32 of the derivative of the nd1.
fn calc_nprimed1(inputs: &Inputs) -> Result<f32, String> {
    let (d1, _) = calc_d1d2(&inputs)?;

    // Get the standard n probability density function value of d1
    let nprimed1 = calc_npdf(d1);
    Ok(nprimed1)
}

/// # Returns
/// f32 of the derivative of the nd2.
fn calc_nprimed2(inputs: &Inputs) -> Result<f32, String> {
    let (_, d2) = calc_d1d2(&inputs)?;

    // Get the standard n probability density function value of d1
    let nprimed2 = calc_npdf(d2);
    Ok(nprimed2)
}

/// Methods for calculating the price, greeks, and implied volatility of an option.
impl Inputs {
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
        s: f32,
        k: f32,
        p: Option<f32>,
        r: f32,
        q: f32,
        t: f32,
        sigma: Option<f32>,
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
    /// f32 of the price of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let price = inputs.calc_price().unwrap();
    /// ```
    pub fn calc_price(&self) -> Result<f32, String> {
        // Calculates the price of the option
        let (nd1, nd2): (f32, f32) = calc_nd1nd2(&self)?;
        let price: f32 = match self.option_type {
            OptionType::Call => f32::max(
                0.0,
                nd1 * self.s * E.powf(-self.q * self.t) - nd2 * self.k * E.powf(-self.r * self.t),
            ),
            OptionType::Put => f32::max(
                0.0,
                nd2 * self.k * E.powf(-self.r * self.t) - nd1 * self.s * E.powf(-self.q * self.t),
            ),
        };
        Ok(price)
    }

    /// Calculates the delta of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the delta of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let delta = inputs.calc_delta().unwrap();
    /// ```
    pub fn calc_delta(&self) -> Result<f32, String> {
        let (nd1, _): (f32, f32) = calc_nd1nd2(&self)?;
        let delta: f32 = match self.option_type {
            OptionType::Call => nd1 * E.powf(-self.q * self.t),
            OptionType::Put => -nd1 * E.powf(-self.q * self.t),
        };
        Ok(delta)
    }

    /// Calculates the gamma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the gamma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let gamma = inputs.calc_gamma().unwrap();
    /// ```
    pub fn calc_gamma(&self) -> Result<f32, String> {
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            panic!("Expected Some(f32) for inputs.sigma, received None")
        };

        let nprimed1: f32 = calc_nprimed1(&self)?;
        let gamma: f32 = E.powf(-self.q * self.t) * nprimed1 / (self.s * sigma * self.t.sqrt());
        Ok(gamma)
    }

    /// Calculates the theta of the option.
    /// Uses 365.25 days in a year for calculations.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of theta per day (not per year).
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let theta = inputs.calc_theta().unwrap();
    /// ```
    pub fn calc_theta(&self) -> Result<f32, String> {
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            panic!("Expected Some(f32) for inputs.sigma, received None")
        };

        let nprimed1: f32 = calc_nprimed1(&self)?;
        let (nd1, nd2): (f32, f32) = calc_nd1nd2(&self)?;

        // Calculation uses 365.25 for f32: Time of days per year.
        let theta: f32 = match self.option_type {
            OptionType::Call => {
                (-(self.s * sigma * E.powf(-self.q * self.t) * nprimed1 / (2.0 * self.t.sqrt()))
                    - self.r * self.k * E.powf(-self.r * self.t) * nd2
                    + self.q * self.s * E.powf(-self.q * self.t) * nd1)
                    / DAYS_PER_YEAR
            }
            OptionType::Put => {
                (-(self.s * sigma * E.powf(-self.q * self.t) * nprimed1 / (2.0 * self.t.sqrt()))
                    + self.r * self.k * E.powf(-self.r * self.t) * nd2
                    - self.q * self.s * E.powf(-self.q * self.t) * nd1)
                    / DAYS_PER_YEAR
            }
        };
        Ok(theta)
    }

    /// Calculates the vega of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the vega of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let vega = inputs.calc_vega().unwrap();
    /// ```
    pub fn calc_vega(&self) -> Result<f32, String> {
        let nprimed1: f32 = calc_nprimed1(&self)?;
        let vega: f32 = 0.01 * self.s * E.powf(-self.q * self.t) * self.t.sqrt() * nprimed1;
        Ok(vega)
    }

    /// Calculates the rho of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the rho of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let rho = inputs.calc_rho().unwrap();
    /// ```
    pub fn calc_rho(&self) -> Result<f32, String> {
        let (_, nd2): (f32, f32) = calc_nd1nd2(&self)?;
        let rho: f32 = match &self.option_type {
            OptionType::Call => 1.0 / 100.0 * self.k * self.t * E.powf(-self.r * self.t) * nd2,
            OptionType::Put => -1.0 / 100.0 * self.k * self.t * E.powf(-self.r * self.t) * nd2,
        };
        Ok(rho)
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
    /// f32 of the epsilon of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let epsilon = inputs.calc_epsilon().unwrap();
    /// ```
    pub fn calc_epsilon(&self) -> Result<f32, String> {
        let (nd1, _) = calc_nd1nd2(&self)?;
        let e_negqt = E.powf(-self.q * self.t);
        let epsilon: f32 = match &self.option_type {
            OptionType::Call => -self.s * self.t * e_negqt * nd1,
            OptionType::Put => self.s * self.t * e_negqt * nd1,
        };
        Ok(epsilon)
    }

    /// Calculates the vanna of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the vanna of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let vanna = inputs.calc_vanna().unwrap();
    /// ```
    pub fn calc_vanna(&self) -> Result<f32, String> {
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            return Err("Expected Some(f32) for inputs.sigma, received None".into());
        };

        let nprimed1 = calc_nprimed1(&self)?;
        let (_, d2) = calc_d1d2(&self)?;
        let vanna: f32 = d2 * E.powf(-self.q * self.t) * nprimed1 * -0.01 / sigma;
        Ok(vanna)
    }

    // /// Calculates the charm of the option.
    // /// # Requires
    // /// s, k, r, q, t, sigma
    // /// # Returns
    // /// f32 of the charm of the option.
    // /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let charm = inputs.calc_charm().unwrap();
    /// ```
    pub fn calc_charm(&self) -> Result<f32, String> {
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            panic!("Expected Some(f32) for inputs.sigma, received None")
        };

        let nprimed1 = calc_nprimed1(&self)?;
        let (nd1, _) = calc_nd1nd2(&self)?;
        let (_, d2) = calc_d1d2(&self)?;
        let e_negqt = E.powf(-self.q * self.t);

        let charm = match &self.option_type {
            OptionType::Call => {
                self.q * e_negqt * nd1
                    - e_negqt
                        * nprimed1
                        * (2.0 * (self.r - self.q) * self.t - d2 * sigma * self.t.sqrt())
                        / (2.0 * self.t * sigma * self.t.sqrt())
            }
            OptionType::Put => {
                -self.q * e_negqt * nd1
                    - e_negqt
                        * nprimed1
                        * (2.0 * (self.r - self.q) * self.t - d2 * sigma * self.t.sqrt())
                        / (2.0 * self.t * sigma * self.t.sqrt())
            }
        };
        Ok(charm)
    }

    /// Calculates the veta of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the veta of the option.
    pub fn calc_veta(&self) -> Result<f32, String> {
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            panic!("Expected Some(f32) for inputs.sigma, received None")
        };

        let nprimed1 = calc_nprimed1(&self)?;
        let (d1, d2) = calc_d1d2(&self)?;
        let e_negqt = E.powf(-self.q * self.t);

        let veta = -self.s
            * e_negqt
            * nprimed1
            * self.t.sqrt()
            * (self.q + ((self.r - self.q) * d1) / (sigma * self.t.sqrt())
                - ((1.0 + d1 * d2) / (2.0 * self.t)));
        Ok(veta)
    }

    /// Calculates the vomma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the vomma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let vomma = inputs.calc_vomma().unwrap();
    /// ```
    pub fn calc_vomma(&self) -> Result<f32, String> {
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            return Err("Expected Some(f32) for inputs.sigma, received None".into());
        };

        let (d1, d2) = calc_d1d2(&self)?;

        let vomma = self.calc_vega()? * ((d1 * d2) / sigma);
        Ok(vomma)
    }

    /// Calculates the speed of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the speed of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let speed = inputs.calc_speed().unwrap();
    /// ```
    pub fn calc_speed(&self) -> Result<f32, String> {
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            return Err("Expected Some(f32) for inputs.sigma, received None".into());
        };

        let (d1, _) = calc_d1d2(&self)?;
        let gamma = self.calc_gamma()?;

        let speed = -gamma / self.s * (d1 / (sigma * self.t.sqrt()) + 1.0);
        Ok(speed)
    }

    /// Calculates the zomma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the zomma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let zomma = inputs.calc_zomma().unwrap();
    /// ```
    pub fn calc_zomma(&self) -> Result<f32, String> {
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            return Err("Expected Some(f32) for inputs.sigma, received None".into());
        };

        let (d1, d2) = calc_d1d2(&self)?;
        let gamma = self.calc_gamma()?;

        let zomma = gamma * ((d1 * d2 - 1.0) / sigma);
        Ok(zomma)
    }

    /// Calculates the color of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the color of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let color = inputs.calc_color().unwrap();
    /// ```
    pub fn calc_color(&self) -> Result<f32, String> {
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            return Err("Expected Some(f32) for inputs.sigma, received None".into());
        };

        let (d1, d2) = calc_d1d2(&self)?;
        let nprimed1 = calc_nprimed1(&self)?;
        let e_negqt = E.powf(-self.q * self.t);

        let color = -e_negqt
            * (nprimed1 / (2.0 * self.s * self.t * sigma * self.t.sqrt()))
            * (2.0 * self.q * self.t
                + 1.0
                + (2.0 * (self.r - self.q) * self.t - d2 * sigma * self.t.sqrt())
                    / (sigma * self.t.sqrt())
                    * d1);
        Ok(color)
    }

    /// Calculates the ultima of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the ultima of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let ultima = inputs.calc_ultima().unwrap();
    /// ```
    pub fn calc_ultima(&self) -> Result<f32, String> {
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            return Err("Expected Some(f32) for inputs.sigma, received None".into());
        };

        let (d1, d2) = calc_d1d2(&self)?;
        let vega = self.calc_vega()?;

        let ultima =
            -vega / sigma.powf(2.0) * (d1 * d2 * (1.0 - d1 * d2) + d1.powf(2.0) + d2.powf(2.0));
        Ok(ultima)
    }

    /// Calculates the dual delta of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the dual delta of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let dual_delta = inputs.calc_dual_delta().unwrap();
    /// ```
    pub fn calc_dual_delta(&self) -> Result<f32, String> {
        let (_, nd2) = calc_nd1nd2(&self)?;
        let e_negqt = E.powf(-self.q * self.t);

        let dual_delta = match self.option_type {
            OptionType::Call => -e_negqt * nd2,
            OptionType::Put => e_negqt * nd2,
        };
        Ok(dual_delta)
    }

    /// Calculates the dual gamma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the dual gamma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let dual_gamma = inputs.calc_dual_gamma().unwrap();
    /// ```
    pub fn calc_dual_gamma(&self) -> Result<f32, String> {
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else {
            return Err("Expected Some(f32) for inputs.sigma, received None".into());
        };

        let nprimed2 = calc_nprimed2(&self)?;
        let e_negqt = E.powf(-self.q * self.t);

        let dual_gamma = e_negqt * (nprimed2 / (self.k * sigma * self.t.sqrt()));
        Ok(dual_gamma)
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
    /// f32 of the implied volatility of the option.
    /// # Example:
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, Some(0.2), 0.05, 0.2, 20.0/365.25, None);
    /// let iv = inputs.calc_iv(0.0001).unwrap();
    /// ```
    pub fn calc_iv(&self, tolerance: f32) -> Result<f32, String> {
        let mut inputs: Inputs = self.clone();

        let p = if let Some(p) = inputs.p {
            p
        } else {
            panic!("inputs.p must contain Some(f32), found None")
        };
        // Initialize estimation of sigma using Brenn and Subrahmanyam (1998) method of calculating initial iv estimation.
        let mut sigma: f32 = (2.0 * PI / inputs.t).sqrt() * (p / inputs.s);
        // Initialize diff to 100 for use in while loop
        let mut diff: f32 = 100.0;

        // Uses Newton Raphson algorithm to calculate implied volatility.
        // Test if the difference between calculated option price and actual option price is > tolerance,
        // if so then iterate until the difference is less than tolerance
        while diff.abs() > tolerance {
            inputs.sigma = Some(sigma);
            diff = inputs.calc_price()? - p;
            sigma -= diff / (inputs.calc_vega()? * 100.0);
        }
        Ok(sigma)
    }
}
