//! # blackscholes_python
//! This library provides an simple, lightweight, and efficient (though not heavily optimized) implementation of the Black-Scholes-Merton model for pricing European options.
//!
//! ## Usage
//! Simply create an instance of the `Inputs` struct and call the desired method.
//!
//! Example:
//! ```
//! let inputs: blackscholes::Inputs = blackscholes::Inputs.new(blackscholes::OptionType::Call, 100.0, 100.0, None, 0.05, 0.02, 20.0 / 365.25, Some(0.2));
//! let price: f64 = inputs.calc_price();
//! ```
//!
//! See the [Github Repo](https://github.com/hayden4r4/blackscholes-rust/tree/python_package) for full source code.  Other implementations such as a WASM crate and a [pure Rust Crate](https://crates.io/crates/blackscholes) are also available.

use statrs::distribution::{Continuous, ContinuousCDF, Normal};
use std::f64::consts::{E, PI};
use std::fmt::{Display, Formatter, Result};

use pyo3::prelude::*;

#[pymodule]
fn blackscholes_python(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<OptionType>()?;
    m.add_class::<Inputs>()?;
    Ok(())
}

/// The type of option to be priced. Call or Put.
#[derive(Debug, Clone, Eq, PartialEq)]
#[pyclass(text_signature = "(Call, Put, /)")]
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

#[pymethods]
impl OptionType {
    /// Creates an instance of the `OptionType` enum.
    /// # Arguments
    /// * `option_type` - The type of option to be priced. Call or Put.
    /// # Example
    /// ```
    /// use blackscholes::OptionType;
    /// let option_type: OptionType = OptionType::Call;
    /// ```
    /// # Returns
    /// An instance of the `OptionType` enum.
    #[new]
    pub fn new(option_type: &str) -> Self {
        match option_type {
            "Call" => OptionType::Call,
            "Put" => OptionType::Put,
            _ => panic!("Option type must be either Call or Put"),
        }
    }
    /// # Returns
    /// A string representation of the `OptionType` enum.
    pub fn __str__(&self) -> String {
        format!("{}", self)
    }
}

/// The inputs to the Black-Scholes-Merton model.
#[derive(Debug, Clone)]
#[pyclass(text_signature = "(option_type, s, k, p, r, q, t, sigma, /)")]
pub struct Inputs {
    /// The type of the option (call or put)
    pub option_type: OptionType,
    /// Stock price
    pub s: f64,
    /// Strike price
    pub k: f64,
    /// Option price
    pub p: Option<f64>,
    /// Risk-free rate
    pub r: f64,
    /// Dividend yield
    pub q: f64,
    /// Time to maturity in years
    pub t: f64,
    /// Volatility
    pub sigma: Option<f64>,
}

impl Display for Inputs {
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
/// # Returns
/// Tuple (f64, f64) of the nd1 and nd2 values for the given inputs.
fn nd1nd2(inputs: &Inputs, normal: bool) -> (f64, f64) {
    let sigma: f64 = match inputs.sigma {
        Some(sigma) => sigma,
        None => panic!("Expected an Option(f64) for inputs.sigma, received None"),
    };

    let nd1nd2 = {
        // Calculating numerator of d1
        let numd1: f64 =
            (inputs.s / inputs.k).ln() + (inputs.r - inputs.q + (sigma.powi(2)) / 2.0) * inputs.t;

        // Calculating denominator of d1 and d2
        let den: f64 = sigma * (inputs.t.sqrt());

        let d1: f64 = numd1 / den;
        let d2: f64 = d1 - den;

        let d1d2: (f64, f64) = (d1, d2);

        // Returns d1 and d2 values if deriving from normal distribution is not necessary
        //  (i.e. gamma, vega, and theta calculations)
        if !normal {
            return d1d2;
        }

        // Creating normal distribution
        let n: Normal = Normal::new(0.0, 1.0).unwrap();

        // Calculates the nd1 and nd2 values
        // Checks if OptionType is Call or Put
        let nd1nd2: (f64, f64) = match inputs.option_type {
            OptionType::Call => (n.cdf(d1d2.0), n.cdf(d1d2.1)),
            OptionType::Put => (n.cdf(-d1d2.0), n.cdf(-d1d2.1)),
        };
        nd1nd2
    };
    nd1nd2
}

/// # Returns
/// f64 of the derivative of the nd1.
fn calc_nprimed1(inputs: &Inputs) -> f64 {
    let (d1, _): (f64, f64) = nd1nd2(&inputs, false);

    // Generate normal probability distribution
    let n: Normal = Normal::new(0.0, 1.0).unwrap();

    // Get the standard normal probability density function value of d1
    let nprimed1: f64 = n.pdf(d1);
    nprimed1
}

/// Methods for calculating the price, greeks, and implied volatility of an option.
#[pymethods]
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
    /// use blackscholes::Inputs;
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20/365.25, Some(0.2));
    /// ```
    /// # Returns
    /// An instance of the `Inputs` struct.
    #[new]
    pub fn new(
        option_type: OptionType,
        s: f64,
        k: f64,
        p: Option<f64>,
        r: f64,
        q: f64,
        t: f64,
        sigma: Option<f64>,
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

    /// # Returns
    /// string representation of the `Inputs` struct.    
    pub fn __str__(&self) -> String {
        format!(
            "OptionType: {}, S: {}, K: {}, P: {}, R: {}, Q: {}, T: {}, Sigma: {}",
            self.option_type,
            self.s,
            self.k,
            match self.p {
                Some(p) => format!("{}", p),
                None => "None".to_string(),
            },
            self.r,
            self.q,
            self.t,
            match self.sigma {
                Some(sigma) => format!("{}", sigma),
                None => "None".to_string(),
            },
        )
    }
    /// Calculates the price of the option.
    /// # Requires
    /// s, k, r, q, t, sigma.
    /// # Returns
    /// f64 of the price of the option.
    /// # Example
    /// ```
    /// use blackscholes::Inputs;
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20/365.25, Some(0.2));
    /// let price = inputs.calc_price();
    /// ```
    pub fn calc_price(&self) -> PyResult<f64> {
        let (nd1, nd2): (f64, f64) = nd1nd2(self, true);
        let price: f64 = match self.option_type {
            OptionType::Call => f64::max(
                0.0,
                nd1 * self.s * E.powf(-self.q * self.t) - nd2 * self.k * E.powf(-self.r * self.t),
            ),
            OptionType::Put => f64::max(
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
    /// f64 of the delta of the option.
    /// # Example
    /// ```
    /// use blackscholes::Inputs;
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20/365.25, Some(0.2));
    /// let delta = inputs.calc_delta();
    /// ```
    pub fn calc_delta(&self) -> PyResult<f64> {
        let (nd1, _): (f64, f64) = nd1nd2(self, true);
        let delta: f64 = match self.option_type {
            OptionType::Call => nd1 * E.powf(-self.q * self.t),
            OptionType::Put => -nd1 * E.powf(-self.q * self.t),
        };
        Ok(delta)
    }

    /// Calculates the gamma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f64 of the gamma of the option.
    /// # Example
    /// ```
    /// use blackscholes::Inputs;
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20/365.25, Some(0.2));
    /// let gamma = inputs.calc_gamma();
    /// ```
    pub fn calc_gamma(&self) -> PyResult<f64> {
        let sigma: f64 = match self.sigma {
            Some(sigma) => sigma,
            None => panic!("Expected an Option(f64) for self.sigma, received None"),
        };

        let nprimed1: f64 = calc_nprimed1(self);
        let gamma: f64 = E.powf(-self.q * self.t) * nprimed1 / (self.s * sigma * self.t.sqrt());
        Ok(gamma)
    }

    /// Calculates the theta of the option.
    /// Uses 365.25 days in a year for calculations.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f64 of theta per day (not per year).
    /// # Example
    /// ```
    /// use blackscholes::Inputs;
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20/365.25, Some(0.2));
    /// let theta = inputs.calc_theta();
    /// ```
    pub fn calc_theta(&self) -> PyResult<f64> {
        let sigma: f64 = match self.sigma {
            Some(sigma) => sigma,
            None => panic!("Expected an Option(f64) for self.sigma, received None"),
        };

        let nprimed1: f64 = calc_nprimed1(self);
        let (nd1, nd2): (f64, f64) = nd1nd2(self, true);

        // Calculation uses 360 for T: Time of days per year.
        let theta: f64 = match self.option_type {
            OptionType::Call => {
                (-(self.s * sigma * E.powf(-self.q * self.t) * nprimed1 / (2.0 * self.t.sqrt()))
                    - self.r * self.k * E.powf(-self.r * self.t) * nd2
                    + self.q * self.s * E.powf(-self.q * self.t) * nd1)
                    / 365.25
            }
            OptionType::Put => {
                (-(self.s * sigma * E.powf(-self.q * self.t) * nprimed1 / (2.0 * self.t.sqrt()))
                    + self.r * self.k * E.powf(-self.r * self.t) * nd2
                    - self.q * self.s * E.powf(-self.q * self.t) * nd1)
                    / 365.25
            }
        };
        Ok(theta)
    }

    /// Calculates the vega of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f64 of the vega of the option.
    /// # Example
    /// ```
    /// use blackscholes::Inputs;
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20/365.25, Some(0.2));
    /// let vega = inputs.calc_vega();
    /// ```
    pub fn calc_vega(&self) -> PyResult<f64> {
        let nprimed1: f64 = calc_nprimed1(self);
        let vega: f64 = 1.0 / 100.0 * self.s * E.powf(-self.q * self.t) * self.t.sqrt() * nprimed1;
        Ok(vega)
    }

    /// Calculates the rho of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f64 of the rho of the option.
    /// # Example
    /// ```
    /// use blackscholes::Inputs;
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20/365.25, Some(0.2));
    /// let rho = inputs.calc_rho();
    /// ```
    pub fn calc_rho(&self) -> PyResult<f64> {
        let (_, nd2): (f64, f64) = nd1nd2(self, true);
        let rho: f64 = match self.option_type {
            OptionType::Call => 1.0 / 100.0 * self.k * self.t * E.powf(-self.r * self.t) * nd2,
            OptionType::Put => -1.0 / 100.0 * self.k * self.t * E.powf(-self.r * self.t) * nd2,
        };
        Ok(rho)
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
    /// f64 of the implied volatility of the option.
    /// # Example:
    /// ```
    /// use blackscholes::Inputs;
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, Some(10), 0.05, 0.02, 20.0 / 365.25, None);
    /// let iv = inputs.calc_iv(0.0001);
    /// ```
    pub fn calc_iv(&self, tolerance: f64) -> PyResult<f64> {
        let mut inputs: Inputs = self.clone();

        let p: f64 = match inputs.p {
            Some(p) => p,
            None => panic!("inputs.p must contain Some(f64), found None"),
        };
        // Initialize estimation of sigma using Brenn and Subrahmanyam (1998) method of calculating initial iv estimation
        let mut sigma: f64 = (2.0 * PI / inputs.t).sqrt() * (p / inputs.s);
        // Initialize diff to 100 for use in while loop
        let mut diff: f64 = 100.0;

        // Uses Newton Raphson algorithm to calculate implied volatility
        // Test if the difference between calculated option price and actual option price is > tolerance
        // If so then iterate until the difference is less than tolerance
        while diff.abs() > tolerance {
            inputs.sigma = Some(sigma);
            diff = inputs.calc_price().unwrap() - p;
            sigma -= diff / (inputs.calc_vega().unwrap() * 100.0);
        }
        Ok(sigma)
    }
}
