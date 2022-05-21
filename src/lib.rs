use crossbeam_channel::bounded;
use crossbeam_utils::thread;
use statrs::distribution::{ContinuousCDF, Normal};
use std::f64::consts::{E, PI};
use std::sync::RwLock;
// use std::thread;

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
    // risk-free rate
    pub r: f64,
    // dividend yield
    pub q: f64,
    // time to maturity as fraction of year
    pub t: f64,
    // volatility
    pub sigma: f64,
}

fn nd1nd2(inputs: &Inputs, normal: bool) -> (f64, f64) {
    // Returns the first and second order moments of the normal distribution

    // Creates Read-write lock on inputs stuct
    let inputs = RwLock::new(inputs);

    // Creating channels with max capacity of one message to communicate between threads
    let (numd1_tx, numd1_rx) = bounded(1);
    let (den_tx, den_rx) = bounded(1);
    let (d1d2_tx, d1d2_rx) = bounded(1);
    let (n_tx, n_rx) = bounded(1);

    // Spawning a scoped thread to avoid the need for arcs (adds overhead as they are stored on heap)
    let nd1nd2 = thread::scope(|s| {
        // Creating thread for first order moment
        s.spawn(|_| {
            // Calculating numerator of the first moment of the normal distribution
            let inputs = inputs.read().unwrap();

            let numd1: f64 = (inputs.s / inputs.k).ln()
                + (inputs.r - inputs.q + (inputs.sigma.powi(2)) / 2.0) * inputs.t;
            // Send the result to the channel
            numd1_tx.send(numd1).unwrap();
        });

        s.spawn(|_| {
            // Calculating denominator of the first and second moment of the normal distribution
            let inputs = inputs.read().unwrap();

            let den: f64 = inputs.sigma * (inputs.t.sqrt());
            den_tx.send(den).unwrap();
        });

        s.spawn(|_| {
            // Calculating the first and second moment of the normal distribution
            let den: f64 = den_rx.recv().unwrap();

            let d1: f64 = numd1_rx.recv().unwrap() / den;
            let d2: f64 = d1 - den;

            let d1d2: (f64, f64) = (d1, d2);
            d1d2_tx.send(d1d2).unwrap();
        });

        s.spawn(|_| {
            // Creating normal distribution
            let n: Normal = Normal::new(0.0, 1.0).unwrap();
            n_tx.send(n).unwrap();
        });

        let n: Normal = n_rx.recv().unwrap();

        let d1d2: (f64, f64) = d1d2_rx.recv().unwrap();

        // Returns first and second moments if deriving from normal distribution is not necessary
        if !normal {
            return d1d2;
        }

        // Calculates the first and second order moments of the normal distribution
        // Checks if OptionType is Call or Put
        let nd1nd2: (f64, f64) = match inputs.read().unwrap().option_type {
            OptionType::Call => (n.cdf(d1d2.0), n.cdf(d1d2.1)),
            OptionType::Put => (n.cdf(-d1d2.0), n.cdf(-d1d2.1)),
        };
        nd1nd2
    })
    .unwrap();

    nd1nd2
}

pub fn calc_price(inputs: Inputs) -> f64 {
    // Returns the price of the option

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
    price
}

pub fn calc_delta(inputs: Inputs) -> f64 {
    let (nd1, _): (f64, f64) = nd1nd2(&inputs, true);
    let delta: f64 = match &inputs.option_type {
        OptionType::Call => nd1 * E.powf(-&inputs.q * &inputs.t),
        OptionType::Put => -nd1 * E.powf(-&inputs.q * &inputs.t),
    };
    delta
}

pub fn calc_nprimed1(inputs: &Inputs) -> f64 {
    // Returns the derivative of the first order moment of the normal distribution
    let (d1, _): (f64, f64) = nd1nd2(&inputs, false);
    let nprimed1: f64 = (1.0 / (2.0 * PI).sqrt()) * E.powf((-d1).powf(2.0) / 2.0);
    nprimed1
}

pub fn calc_gamma(inputs: Inputs) -> f64 {
    // Calculates the gamma of an option
    let nprimed1: f64 = calc_nprimed1(&inputs);
    let gamma: f64 =
        E.powf(-&inputs.q * &inputs.t) * nprimed1 / (&inputs.s * &inputs.sigma * &inputs.t.sqrt());
    gamma
}

pub fn calc_theta(inputs: Inputs) -> f64 {
    // Calculates the theta of the option

    let nprimed1: f64 = calc_nprimed1(&inputs);
    let (nd1, nd2): (f64, f64) = nd1nd2(&inputs, true);

    // Calculation uses 360 for T: Time of days per year.
    let theta: f64 = match &inputs.option_type {
        OptionType::Call => {
            (-(&inputs.s * &inputs.sigma * E.powf(-&inputs.q * &inputs.t) * nprimed1
                / (2.0 * &inputs.t.sqrt()))
                - &inputs.r * &inputs.k * E.powf(-&inputs.r * &inputs.t) * nd2
                + &inputs.q * &inputs.s * E.powf(-&inputs.q * &inputs.t) * nd1)
                / 360.0
        }
        OptionType::Put => {
            (-(&inputs.s * &inputs.sigma * E.powf(-&inputs.q * &inputs.t) * nprimed1
                / (2.0 * &inputs.t.sqrt()))
                + &inputs.r * &inputs.k * E.powf(-&inputs.r * &inputs.t) * nd2
                - &inputs.q * &inputs.s * E.powf(-&inputs.q * &inputs.t) * nd1)
                / 360.0
        }
    };
    theta
}

pub fn calc_vega(inputs: Inputs) -> f64 {
    // Calculates the vega of the option

    let nprimed1: f64 = calc_nprimed1(&inputs);
    let vega: f64 =
        1.0 / 100.0 * &inputs.s * E.powf(-&inputs.q * &inputs.t) * &inputs.t.sqrt() * nprimed1;
    vega
}

pub fn calc_rho(inputs: Inputs) -> f64 {
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
    rho
}
