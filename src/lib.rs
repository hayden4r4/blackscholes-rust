use crossbeam_channel::bounded;
use statrs::distribution::{ContinuousCDF, Normal};
use std::f64::consts::E;
use std::sync::{Arc, RwLock};
use std::thread;

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

fn nd1nd2(inputs: Arc<RwLock<Inputs>>) -> (f64, f64) {
    // Returns the first and second order moments of the normal distribution

    // Creating a vector to hold threads before joining
    let mut handles = Vec::new();

    // Creating a channel with max capacity of one message to communicate between threads
    let (numd1_tx, numd1_rx) = bounded(1);
    // Creating mutex for inputs
    let inputs1: Arc<RwLock<Inputs>> = Arc::clone(&inputs);
    // Creating thread for first order moment
    let handle_numd1 = thread::spawn(move || {
        // Calculating numerator of the first moment of the normal distribution
        // Obtains a read lock on inputs1
        let inputs = inputs1.read().unwrap();

        let numd1: f64 = (inputs.s / inputs.k).ln()
            + (inputs.r - inputs.q + (inputs.sigma.powi(2)) / 2.0) * inputs.t;
        // Send the result to the channel
        numd1_tx.send(numd1).unwrap();
    });
    // Pushing the handle to the handles vector
    handles.push(handle_numd1);

    let (den_tx, den_rx) = bounded(1);
    let inputs2: Arc<RwLock<Inputs>> = Arc::clone(&inputs);
    let handle_den = thread::spawn(move || {
        // Calculating denominator of the first and second moment of the normal distribution
        let inputs = inputs2.read().unwrap();

        let den: f64 = inputs.sigma * (inputs.t.sqrt());
        den_tx.send(den).unwrap();
    });
    handles.push(handle_den);

    let (d1d2_tx, d1d2_rx) = bounded(1);
    let handle_d1d2 = thread::spawn(move || {
        // Calculating the first and second moment of the normal distribution
        let den: f64 = den_rx.recv().unwrap();

        let d1: f64 = numd1_rx.recv().unwrap() / den;
        let d2: f64 = d1 - den;

        let d1d2: (f64, f64) = (d1, d2);
        d1d2_tx.send(d1d2).unwrap();
    });
    handles.push(handle_d1d2);

    let (n_tx, n_rx) = bounded(1);
    let handle_n = thread::spawn(move || {
        // Creating normal distribution

        let n: Normal = Normal::new(0.0, 1.0).unwrap();
        n_tx.send(n).unwrap();
    });
    handles.push(handle_n);

    // Calculates the first and second order moments of the normal distribution
    let n: Normal = n_rx.recv().unwrap();

    let d1d2: (f64, f64) = d1d2_rx.recv().unwrap();

    // Checks if OptionType is Call or Put
    let nd1nd2: (f64, f64) = match inputs.read().unwrap().option_type {
        OptionType::Call => (n.cdf(d1d2.0), n.cdf(d1d2.1)),
        OptionType::Put => (n.cdf(-d1d2.0), n.cdf(-d1d2.1)),
    };

    for handle in handles {
        handle.join().unwrap();
    }

    nd1nd2
}

pub fn calc_price(inputs: Inputs) -> f64 {
    // Returns the price of the option

    // Calculates the price of the option
    let inputs_arc = Arc::new(RwLock::new(inputs));
    let (nd1, nd2): (f64, f64) = nd1nd2(Arc::clone(&inputs_arc));
    let inputs = inputs_arc.read().unwrap();
    let price: f64 = match inputs.option_type {
        OptionType::Call => f64::max(
            0.0,
            nd1 * inputs.s * E.powf(-inputs.q * inputs.t)
                - nd2 * inputs.k * E.powf(-inputs.r * inputs.t),
        ),
        OptionType::Put => f64::max(
            0.0,
            nd2 * inputs.k * E.powf(-inputs.r * inputs.t)
                - nd1 * inputs.s * E.powf(-inputs.q * inputs.t),
        ),
    };
    price
}
