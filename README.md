# BlackScholes
A library providing Black-Scholes option pricing, Greek calculations, and implied-volatility solver.

[![Crates.io](https://img.shields.io/crates/v/blackscholes.svg)](https://crates.io/crates/blackscholes)
[![Documentation](https://docs.rs/blackscholes/badge.svg)](https://docs.rs/blackscholes)
[![Benchmarks](https://img.shields.io/badge/Benchmarks-GitHub%20Pages-blue)](https://przemyslawolszewski.github.io/bs-rs/)

A Black-Scholes option pricing, Greek calculation, and implied volatility calculation library.

The library handles both European and American style options for the following option types:
- Vanilla Put/Call
- Binary Put/Call
- Binary OT Range (In/Out)
- Barrier


# Performance Optimization

This library is optimized for both single-option pricing and high-throughput batch processing. We've implemented a comprehensive benchmarking infrastructure to measure and improve performance.

## Key Performance Characteristics

- **Single Option Pricing**: ~35-40 ns per option
- **Rational Pricing Method**: ~55-65 ns per option
- **Delta Calculation**: ~30-35 ns per option
- **Gamma Calculation**: ~14-15 ns per option
- **Batch Processing**: Scales linearly up to large batch sizes
- **All Greeks Calculation**: ~2 ms per 1000 options

## Benchmarking Infrastructure

The library includes a comprehensive benchmarking system for performance tracking:

- **Interactive Charts**: Professional benchmark visualizations on [GitHub Pages](https://przemyslawolszewski.github.io/bs-rs/)
- **Automated Regression Detection**: CI-integrated tests that fail on performance regressions (>10% threshold)
- **Historical Tracking**: Continuous monitoring of performance trends over time
- **Pull Request Comments**: Automatic performance comparison comments on PRs

View live benchmark results at: https://przemyslawolszewski.github.io/bs-rs/

## Usage Examples

```rust
use blackscholes::{Inputs, OptionType, Pricing, Greeks, ImpliedVolatility};

// Basic option pricing
let inputs = Inputs::new(
    OptionType::Call,   // Call option
    100.0,              // Spot price
    100.0,              // Strike price
    None,               // Option price (not needed for pricing)
    0.05,               // Risk-free rate
    0.01,               // Dividend yield
    1.0,                // Time to maturity (in years)
    Some(0.2),          // Volatility
);

// Calculate option price
let price = inputs.calc_price().unwrap();
println!("Option price: {}", price);

// Calculate option Greeks
let delta = inputs.calc_delta().unwrap();
let gamma = inputs.calc_gamma().unwrap();
let theta = inputs.calc_theta().unwrap();
let vega = inputs.calc_vega().unwrap();
let rho = inputs.calc_rho().unwrap();

println!("Delta: {}, Gamma: {}, Vega: {}", delta, gamma, vega);

// Calculate implied volatility from price
let mut iv_inputs = Inputs::new(
    OptionType::Call,
    100.0,
    100.0,
    Some(10.0),  // Option price
    0.05,
    0.01,
    1.0,
    None,        // Volatility is what we're solving for
);

let iv = iv_inputs.calc_iv(0.0001).unwrap();
println!("Implied volatility: {}", iv);
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

# blackscholes
[![Crates.io](https://img.shields.io/crates/v/blackscholes)](https://crates.io/crates/blackscholes)
[![Docs.rs](https://docs.rs/blackscholes/badge.svg)](https://docs.rs/blackscholes)
![License](https://img.shields.io/crates/l/blackscholes)
  
This library provides a simple, lightweight, and efficient (though not heavily optimized) implementation of the Black-Scholes-Merton model for pricing European options.  
  
Includes all first, second, and third order Greeks.  

Implements both:  

- calc_iv() in the ImpliedVolatility trait which uses [Modified Corrado-Miller by Piotr Płuciennik (2007)](https://sin.put.poznan.pl/files/download/37938) for the initial volatility guess and the Newton Raphson algorithm to solve for the implied volatility.
- calc_rational_iv() in the ImpliedVolatility trait which uses "Let's be rational" method from ["Let's be rational" (2016) by Peter Jackel](http://www.jaeckel.org/LetsBeRational.pdf).  Utilizing Jackel's C++ implementation to get convergence within 2 iterations with 64-bit floating point accuracy.
  
## Usage  
  
View the [docs](https://docs.rs/blackscholes) for usage and examples.  
  
**Other packages available:**  
Python: [Pypi](https://pypi.org/project/blackscholes-python/)  
WASM: [npm](https://www.npmjs.com/package/@haydenr4/blackscholes_wasm)  
