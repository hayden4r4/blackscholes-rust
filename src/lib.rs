use statrs::distribution::{Continuous, ContinuousCDF, Normal};
use std::f64::consts::{E, PI};

#[derive(Clone)]
pub enum OptionType {
    Call,
    Put,
}

pub struct Inputs {
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

fn nd1nd2(inputs: &Inputs, normal: bool, sigma: Option<f64>) -> (f64, f64) {
    // Returns the first and second order moments of the normal distribution
    let sigma: f64 = match sigma {
        Some(sigma) => sigma,
        None => match inputs.sigma {
            Some(sigma) => sigma,
            None => panic!("Expected an Option(f64) for inputs.sigma, received None"),
        },
    };

    let nd1nd2 = {
        // Calculating numerator of the first moment of the normal distribution
        let numd1: f64 =
            (inputs.s / inputs.k).ln() + (inputs.r - inputs.q + (sigma.powi(2)) / 2.0) * inputs.t;

        // Calculating denominator of the first and second moment of the normal distribution
        let den: f64 = sigma * (inputs.t.sqrt());

        let d1: f64 = numd1 / den;
        let d2: f64 = d1 - den;

        let d1d2: (f64, f64) = (d1, d2);

        // Returns first and second moments if deriving from normal distribution is not necessary
        if !normal {
            return d1d2;
        }

        // Creating normal distribution
        let n: Normal = Normal::new(0.0, 1.0).unwrap();

        // Calculates the first and second order moments of the normal distribution
        // Checks if OptionType is Call or Put
        let nd1nd2: (f64, f64) = match inputs.option_type {
            OptionType::Call => (n.cdf(d1d2.0), n.cdf(d1d2.1)),
            OptionType::Put => (n.cdf(-d1d2.0), n.cdf(-d1d2.1)),
        };
        nd1nd2
    };
    nd1nd2
}

fn _calc_price(inputs: &Inputs, custom_sigma: Option<f64>) -> f64 {
    // Returns the price of the option

    // Calculates the price of the option
    let (nd1, nd2): (f64, f64) = nd1nd2(&inputs, true, custom_sigma);
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
    price
}

pub fn calc_price(inputs: &Inputs) -> f64 {
    // Returns the price of the option (abstraction for _calc_price)

    _calc_price(inputs, None)
}

pub fn calc_delta(inputs: &Inputs) -> f64 {
    let (nd1, _): (f64, f64) = nd1nd2(&inputs, true, None);
    let delta: f64 = match &inputs.option_type {
        OptionType::Call => nd1 * E.powf(-&inputs.q * &inputs.t),
        OptionType::Put => -nd1 * E.powf(-&inputs.q * &inputs.t),
    };
    delta
}

fn calc_nprimed1(inputs: &Inputs, sigma: Option<f64>) -> f64 {
    // Returns the derivative of the first order moment of the normal distribution

    let (d1, _): (f64, f64) = nd1nd2(&inputs, false, sigma);

    // Generate normal probability distribution
    let n: Normal = Normal::new(0.0, 1.0).unwrap();

    // Get the standard normal probability density function value of d1
    let nprimed1: f64 = n.pdf(d1);
    nprimed1
}

pub fn calc_gamma(inputs: &Inputs) -> f64 {
    // Calculates the gamma of the option

    let sigma: f64 = match &inputs.sigma {
        Some(sigma) => *sigma,
        None => panic!("Expected an Option(f64) for inputs.sigma, received None"),
    };

    let nprimed1: f64 = calc_nprimed1(&inputs, None);
    let gamma: f64 =
        E.powf(-&inputs.q * &inputs.t) * nprimed1 / (&inputs.s * sigma * &inputs.t.sqrt());
    gamma
}

pub fn calc_theta(inputs: &Inputs) -> f64 {
    // Calculates the theta of the option

    let sigma: f64 = match &inputs.sigma {
        Some(sigma) => *sigma,
        None => panic!("Expected an Option(f64) for inputs.sigma, received None"),
    };

    let nprimed1: f64 = calc_nprimed1(&inputs, None);
    let (nd1, nd2): (f64, f64) = nd1nd2(&inputs, true, None);

    // Calculation uses 360 for T: Time of days per year.
    let theta: f64 = match &inputs.option_type {
        OptionType::Call => {
            (-(&inputs.s * sigma * E.powf(-&inputs.q * &inputs.t) * nprimed1
                / (2.0 * &inputs.t.sqrt()))
                - &inputs.r * &inputs.k * E.powf(-&inputs.r * &inputs.t) * nd2
                + &inputs.q * &inputs.s * E.powf(-&inputs.q * &inputs.t) * nd1)
                / 360.0
        }
        OptionType::Put => {
            (-(&inputs.s * sigma * E.powf(-&inputs.q * &inputs.t) * nprimed1
                / (2.0 * &inputs.t.sqrt()))
                + &inputs.r * &inputs.k * E.powf(-&inputs.r * &inputs.t) * nd2
                - &inputs.q * &inputs.s * E.powf(-&inputs.q * &inputs.t) * nd1)
                / 360.0
        }
    };
    theta
}

fn _calc_vega(inputs: &Inputs, sigma: Option<f64>) -> f64 {
    // Calculates the vega of the option

    let sigma: f64 = match sigma {
        Some(sigma) => sigma,
        None => match inputs.sigma {
            Some(sigma) => sigma,
            None => panic!("Expected an Option(f64) for inputs.sigma, received None"),
        },
    };

    let nprimed1: f64 = calc_nprimed1(&inputs, Some(sigma));
    let vega: f64 =
        1.0 / 100.0 * &inputs.s * E.powf(-&inputs.q * &inputs.t) * &inputs.t.sqrt() * nprimed1;
    vega
}

pub fn calc_vega(inputs: &Inputs) -> f64 {
    // Calculates the vega of the option (abstraction for _calc_vega)

    _calc_vega(inputs, None)
}

pub fn calc_rho(inputs: &Inputs) -> f64 {
    // Calculates the rho of the option

    let (_, nd2): (f64, f64) = nd1nd2(&inputs, true, None);
    let rho: f64 = match &inputs.option_type {
        OptionType::Call => {
            1.0 / 100.0 * &inputs.k * &inputs.t * E.powf(-&inputs.r * &inputs.t) * nd2
        }
        OptionType::Put => {
            -1.0 / 100.0 * &inputs.k * &inputs.t * E.powf(-&inputs.r * &inputs.t) * nd2
        }
    };
    rho
}

pub fn calc_iv(inputs: &Inputs, tolerance: f64) -> f64 {
    // Calculates the implied volatility of the option
    // Tolerance is the max error allowed for the implied volatility,
    // the lower the tolerance the more iterations will be required.
    // Recommended to be a value between 0.001 - 0.0001 for highest efficiency/accuracy

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
        diff = _calc_price(&inputs, Some(sigma)) - p;
        sigma -= diff / (_calc_vega(&inputs, Some(sigma)) * 100.0);
    }
    sigma
}
