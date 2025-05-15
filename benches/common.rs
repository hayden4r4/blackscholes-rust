use std::f64::consts::PI;
use blackscholes::{Inputs, OptionType};
use rand::prelude::*;
use rand::distributions::{Distribution, Standard, Uniform};

/// Defines different option moneyness categories
#[derive(Debug, Clone, Copy)]
pub enum Moneyness {
    DeepInTheMoney,
    InTheMoney,
    AtTheMoney,
    OutOfTheMoney,
    DeepOutOfTheMoney,
}

/// Defines time to maturity ranges
#[derive(Debug, Clone, Copy)]
pub enum TimeToMaturity {
    ShortTerm,   // Less than 30 days
    MediumTerm,  // 30-90 days
    LongTerm,    // More than 90 days
}

/// Defines volatility levels
#[derive(Debug, Clone, Copy)]
pub enum VolatilityLevel {
    Low,     // < 15%
    Medium,  // 15-30%
    High,    // > 30%
}

/// Generates a standard set of inputs for benchmarking
pub fn generate_standard_inputs(
    option_type: OptionType,
    moneyness: Moneyness,
    maturity: TimeToMaturity,
    vol_level: VolatilityLevel,
) -> Inputs {
    let spot = 100.0;
    
    // Set strike based on moneyness
    let strike = match moneyness {
        Moneyness::DeepInTheMoney => if option_type == OptionType::Call { 70.0 } else { 130.0 },
        Moneyness::InTheMoney => if option_type == OptionType::Call { 90.0 } else { 110.0 },
        Moneyness::AtTheMoney => 100.0,
        Moneyness::OutOfTheMoney => if option_type == OptionType::Call { 110.0 } else { 90.0 },
        Moneyness::DeepOutOfTheMoney => if option_type == OptionType::Call { 130.0 } else { 70.0 },
    };
    
    // Set time to maturity in years
    let time = match maturity {
        TimeToMaturity::ShortTerm => 15.0 / 365.25,
        TimeToMaturity::MediumTerm => 60.0 / 365.25,
        TimeToMaturity::LongTerm => 180.0 / 365.25,
    };
    
    // Set volatility level
    let vol = match vol_level {
        VolatilityLevel::Low => 0.10,
        VolatilityLevel::Medium => 0.20,
        VolatilityLevel::High => 0.40,
    };
    
    Inputs::new(
        option_type,
        spot,
        strike,
        None,
        0.05,  // Risk-free rate
        0.01,  // Dividend yield
        time,
        Some(vol),
    )
}

/// Generates a batch of random inputs for benchmarking
pub fn generate_random_inputs(size: usize, rng: &mut ThreadRng) -> Vec<Inputs> {
    let mut inputs = Vec::with_capacity(size);
    
    // Define distributions for parameters
    let spot_dist = Uniform::from(50.0..150.0);
    let strike_dist = Uniform::from(50.0..150.0);
    let rate_dist = Uniform::from(0.0..0.10);
    let div_dist = Uniform::from(0.0..0.05);
    let time_dist = Uniform::from(1.0/365.25..1.0);
    let vol_dist = Uniform::from(0.05..0.50);
    
    for _ in 0..size {
        let option_type = if rng.gen::<bool>() { OptionType::Call } else { OptionType::Put };
        
        inputs.push(Inputs::new(
            option_type,
            spot_dist.sample(rng),
            strike_dist.sample(rng),
            None,
            rate_dist.sample(rng),
            div_dist.sample(rng),
            time_dist.sample(rng),
            Some(vol_dist.sample(rng)),
        ));
    }
    
    inputs
}

/// Struct to hold parameters in struct-of-arrays (SoA) format for SIMD processing
#[derive(Debug, Clone)]
pub struct InputsSoA {
    pub option_types: Vec<OptionType>,
    pub spots: Vec<f64>,
    pub strikes: Vec<f64>,
    pub rates: Vec<f64>,
    pub dividends: Vec<f64>,
    pub times: Vec<f64>,
    pub volatilities: Vec<f64>,
}

impl InputsSoA {
    /// Creates a new InputsSoA with preallocated capacity
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            option_types: Vec::with_capacity(capacity),
            spots: Vec::with_capacity(capacity),
            strikes: Vec::with_capacity(capacity),
            rates: Vec::with_capacity(capacity),
            dividends: Vec::with_capacity(capacity),
            times: Vec::with_capacity(capacity),
            volatilities: Vec::with_capacity(capacity),
        }
    }
    
    /// Convert Vec<Inputs> to InputsSoA format
    pub fn from_inputs(inputs: &[Inputs]) -> Self {
        let mut result = Self::with_capacity(inputs.len());
        
        for input in inputs {
            result.option_types.push(input.option_type);
            result.spots.push(input.s);
            result.strikes.push(input.k);
            result.rates.push(input.r);
            result.dividends.push(input.q);
            result.times.push(input.t);
            result.volatilities.push(input.sigma.unwrap_or(0.0));
        }
        
        result
    }
    
    /// Generate random SoA inputs directly
    pub fn random(size: usize, rng: &mut ThreadRng) -> Self {
        let mut result = Self::with_capacity(size);
        
        // Define distributions for parameters
        let spot_dist = Uniform::from(50.0..150.0);
        let strike_dist = Uniform::from(50.0..150.0);
        let rate_dist = Uniform::from(0.0..0.10);
        let div_dist = Uniform::from(0.0..0.05);
        let time_dist = Uniform::from(1.0/365.25..1.0);
        let vol_dist = Uniform::from(0.05..0.50);
        
        for _ in 0..size {
            let option_type = if rng.gen::<bool>() { OptionType::Call } else { OptionType::Put };
            
            result.option_types.push(option_type);
            result.spots.push(spot_dist.sample(rng));
            result.strikes.push(strike_dist.sample(rng));
            result.rates.push(rate_dist.sample(rng));
            result.dividends.push(div_dist.sample(rng));
            result.times.push(time_dist.sample(rng));
            result.volatilities.push(vol_dist.sample(rng));
        }
        
        result
    }
    
    /// Get the size of the SoA
    pub fn len(&self) -> usize {
        self.spots.len()
    }
    
    /// Check if the SoA is empty
    pub fn is_empty(&self) -> bool {
        self.spots.is_empty()
    }
}

/// Standard benchmark sizes
#[derive(Debug, Clone, Copy)]
pub enum BatchSize {
    Tiny = 10,
    Small = 100,
    Medium = 1_000,
    Large = 10_000,
    Huge = 100_000,
    Massive = 1_000_000,
}

/// Gets appropriate sample sizes and measurement time based on batch size
pub fn get_sample_config(batch_size: usize) -> (usize, std::time::Duration) {
    use std::time::Duration;
    
    match batch_size {
        0..=100 => (100, Duration::from_secs(5)),
        101..=1_000 => (50, Duration::from_secs(5)),
        1_001..=10_000 => (20, Duration::from_secs(10)),
        10_001..=100_000 => (10, Duration::from_secs(10)),
        _ => (3, Duration::from_secs(20)),
    }
} 