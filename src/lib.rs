use statrs::distribution::{Continuous, ContinuousCDF, Normal};
use std::f64::consts::{E, PI};
use std::fmt::{Display, Formatter, Result};

use pyo3::prelude::*;

#[pymodule]
fn blackscholes_python(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<OptionType>()?;
    m.add_class::<Inputs>()?;
    m.add_class::<Price>()?;
    m.add_class::<Greeks>()?;
    m.add_class::<Volatility>()?;
    Ok(())
}

#[derive(Debug, Clone, Eq, PartialEq)]
#[pyclass(text_signature= "(Call, Put, /)")]
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
    #[new]
    pub fn new(option_type: &str) -> Self {
        match option_type {
            "Call" => OptionType::Call,
            "Put" => OptionType::Put,
            _ => panic!("Option type must be either Call or Put"),
        }
    }
    pub fn __str__(&self) -> String {
        format!("{}", self)
    }
}

#[derive(Debug, Clone)]
#[pyclass(text_signature = "(option_type, s, k, p, r, q, t, sigma, /)")]
pub struct Inputs {
    // The type of the option (call or put)
    pub option_type: OptionType,
    // Stock price
    pub s: f64,
    // Strike price
    pub k: f64,
    // Option price
    pub p: Option<f64>,
    // Risk-free rate
    pub r: f64,
    // Dividend yield
    pub q: f64,
    // Time to maturity as fraction of year
    pub t: f64,
    // Volatility
    pub sigma: Option<f64>,
}

#[pymethods]
impl Inputs {
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

    pub fn __str__ (&self) -> String {
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

fn nd1nd2(inputs: &Inputs, normal: bool) -> (f64, f64) {
    // Returns the nd1 and nd2 values for the given inputs

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

fn calc_nprimed1(inputs: &Inputs) -> f64 {
    // Returns the derivative of the nd1

    let (d1, _): (f64, f64) = nd1nd2(&inputs, false);

    // Generate normal probability distribution
    let n: Normal = Normal::new(0.0, 1.0).unwrap();

    // Get the standard normal probability density function value of d1
    let nprimed1: f64 = n.pdf(d1);
    nprimed1
}

#[pyclass]
pub struct Price {}

#[pymethods]
impl Price {
    #[staticmethod]
    pub fn calc_price(inputs: &Inputs) -> PyResult<f64> {
        // Returns the price of the option
        // Requires s k r q t sigma

        // Calculates the price of the option
        let (nd1, nd2): (f64, f64) = nd1nd2(&inputs, true);
        let price: f64 = match &inputs.option_type {
            OptionType::Call => f64::max(
                0.0,
                nd1 * &inputs.s * E.powf(-&inputs.q * &inputs.t)
                    - nd2 * &inputs.k * E.powf(-&inputs.r * &inputs.t),
            ),
            OptionType::Put => f64::max(
                0.0,
                nd2 * inputs.k * E.powf(-inputs.r * inputs.t)
                    - nd1 * inputs.s * E.powf(-inputs.q * inputs.t),
            ),
        };
        Ok(price)
    }
}

#[pyclass]
pub struct Greeks {}

#[pymethods]
impl Greeks {
    #[staticmethod]
    // Requires s k r q t sigma
    pub fn calc_delta(inputs: &Inputs) -> PyResult<f64> {
        // Calculates the delta of the option

        let (nd1, _): (f64, f64) = nd1nd2(&inputs, true);
        let delta: f64 = match &inputs.option_type {
            OptionType::Call => nd1 * E.powf(-&inputs.q * &inputs.t),
            OptionType::Put => -nd1 * E.powf(-&inputs.q * &inputs.t),
        };
        Ok(delta)
    }

    #[staticmethod]
    pub fn calc_gamma(inputs: &Inputs) -> PyResult<f64> {
        // Calculates the gamma of the option

        let sigma: f64 = match &inputs.sigma {
            Some(sigma) => *sigma,
            None => panic!("Expected an Option(f64) for inputs.sigma, received None"),
        };

        let nprimed1: f64 = calc_nprimed1(&inputs);
        let gamma: f64 =
            E.powf(-&inputs.q * &inputs.t) * nprimed1 / (&inputs.s * sigma * &inputs.t.sqrt());
        Ok(gamma)
    }

    #[staticmethod]
    pub fn calc_theta(inputs: &Inputs) -> PyResult<f64> {
        // Calculates the theta of the option

        let sigma: f64 = match &inputs.sigma {
            Some(sigma) => *sigma,
            None => panic!("Expected an Option(f64) for inputs.sigma, received None"),
        };

        let nprimed1: f64 = calc_nprimed1(&inputs);
        let (nd1, nd2): (f64, f64) = nd1nd2(&inputs, true);

        // Calculation uses 360 for T: Time of days per year.
        let theta: f64 = match &inputs.option_type {
            OptionType::Call => {
                (-(&inputs.s * sigma * E.powf(-&inputs.q * &inputs.t) * nprimed1
                    / (2.0 * &inputs.t.sqrt()))
                    - &inputs.r * &inputs.k * E.powf(-&inputs.r * &inputs.t) * nd2
                    + &inputs.q * &inputs.s * E.powf(-&inputs.q * &inputs.t) * nd1)
                    / 365.25
            }
            OptionType::Put => {
                (-(&inputs.s * sigma * E.powf(-&inputs.q * &inputs.t) * nprimed1
                    / (2.0 * &inputs.t.sqrt()))
                    + &inputs.r * &inputs.k * E.powf(-&inputs.r * &inputs.t) * nd2
                    - &inputs.q * &inputs.s * E.powf(-&inputs.q * &inputs.t) * nd1)
                    / 365.25
            }
        };
        Ok(theta)
    }

    #[staticmethod]
    pub fn calc_vega(inputs: &Inputs) -> PyResult<f64> {
        // Calculates the vega of the option

        let nprimed1: f64 = calc_nprimed1(&inputs);
        let vega: f64 =
            1.0 / 100.0 * &inputs.s * E.powf(-&inputs.q * &inputs.t) * &inputs.t.sqrt() * nprimed1;
        Ok(vega)
    }

    #[staticmethod]
    pub fn calc_rho(inputs: &Inputs) -> PyResult<f64> {
        // Calculates the rho of the option

        let (_, nd2): (f64, f64) = nd1nd2(&inputs, true);
        let rho: f64 = match &inputs.option_type {
            OptionType::Call => {
                1.0 / 100.0 * &inputs.k * &inputs.t * E.powf(-&inputs.r * &inputs.t) * nd2
            }
            OptionType::Put => {
                -1.0 / 100.0 * &inputs.k * &inputs.t * E.powf(-&inputs.r * &inputs.t) * nd2
            }
        };
        Ok(rho)
    }
}

#[pyclass]
pub struct Volatility {}

#[pymethods]
impl Volatility {
    #[staticmethod]
    pub fn calc_iv(inputs: &mut Inputs, tolerance: f64) -> PyResult<f64> {
        // Calculates the implied volatility of the option
        // Tolerance is the max error allowed for the implied volatility,
        // the lower the tolerance the more iterations will be required.
        // Recommended to be a value between 0.001 - 0.0001 for highest efficiency/accuracy
        // Requires s k r q t price

        let p: f64 = match &inputs.p {
            Some(p) => *p,
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
            diff = Price::calc_price(&inputs).unwrap() - p;
            sigma -= diff / (Greeks::calc_vega(&inputs).unwrap() * 100.0);
        }
        Ok(sigma)
    }
}
