use polars::prelude::*;
use statrs::distribution::{ContinuousCDF, Normal};
use std::f64::consts::E;
use std::fs::File;

#[derive(Clone)]
enum OptionType {
    Call,
    Put,
}

struct Inputs {
    option_type: OptionType,
    // Stock price
    s: f64,
    // Strike price
    k: f64,
    // risk-free rate
    r: f64,
    // time to maturity as fraction of year
    t: f64,
    // volatility
    sigma: f64,
}

fn load_returns() -> Result<DataFrame> {
    // Loads return data from csv file
    let file = File::open("SP500_returns.csv")?;

    CsvReader::new(file)
        .infer_schema(None)
        .has_header(true)
        .finish()
}

fn std_dev_log_return(returns: Series) -> f64 {
    // Returns the std dev of log returns
    // Calculates the standard deviation of log returns
    let mut returns_vec: Vec<f64> = Vec::new();
    let mut mean: f64 = 0.0;
    for i in 0..returns.len() {
        if i <= returns.len() - 2 {
            returns_vec.push(
                (returns.str_value(i + 1).to_string().parse::<f64>().unwrap()
                    / returns.str_value(i).to_string().parse::<f64>().unwrap())
                .ln(),
            );
            mean += returns_vec[i];
        }
    }

    mean /= returns_vec.len() as f64;
    // Calculating variance of log returns
    let mut variance: f64 = 0.0;
    for i in 0..returns_vec.len() {
        variance += (returns_vec[i] - mean).powi(2);
    }
    variance /= returns_vec.len() as f64;
    // Returns square root of variance = std dev.
    variance.sqrt()
}

fn nd1nd2(inputs: &Inputs) -> (f64, f64) {
    // Returns the first and second order moments of the normal distribution

    // Calculating numerator of the first moment of the normal distribution
    let numd1: f64 =
        (inputs.s / inputs.k).ln() + (inputs.r + (inputs.sigma.powi(2)) / 2.0) * inputs.t;
    // Calculating numerator of the second moment of the normal distribution
    let numd2: f64 =
        (inputs.s / inputs.k).ln() + (inputs.r - (inputs.sigma.powi(2)) / 2.0) * inputs.t;
    // Calculating denominator of the first and second moment of the normal distribution
    let den: f64 = inputs.sigma * (inputs.t.sqrt());
    let d1d2: (f64, f64) = (numd1 / den, numd2 / den);
    // Creating normal distribution
    let n: Normal = Normal::new(0.0, 1.0).unwrap();
    // Calculates the first and second order moments of the normal distribution
    let nd1nd2: (f64, f64) = match inputs.option_type {
        OptionType::Call => (n.cdf(d1d2.0), n.cdf(d1d2.1)),
        OptionType::Put => (n.cdf(-d1d2.0), n.cdf(-d1d2.1)),
    };
    nd1nd2
}

fn calc_price(inputs: &Inputs) -> f64 {
    // Returns the price of the option

    // Calculates the price of the option
    let (nd1, nd2): (f64, f64) = nd1nd2(inputs);
    let price: f64 = match inputs.option_type {
        OptionType::Call => f64::max(
            0.0,
            nd1 * inputs.s - nd2 * inputs.k * E.powf(-inputs.r * inputs.t),
        ),
        OptionType::Put => f64::max(
            0.0,
            nd2 * inputs.k * E.powf(-inputs.r * inputs.t) - nd1 * inputs.s,
        ),
    };
    price
}

fn main() {
    let returns: Series = load_returns().unwrap()["Close"].clone();
    let std_dev: f64 = std_dev_log_return(returns);
    let inputs: Inputs = Inputs {
        option_type: OptionType::Call,
        s: 4123.34,
        k: 4120.0,
        r: 0.0312,
        t: 176.0 / 360.0,
        sigma: std_dev,
    };
    let price: f64 = calc_price(&inputs);
    println!("{}", price);
}
