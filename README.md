# Black-Scholes Rust

A Black-Scholes pricing model built in Rust. Will calculate the price, greeks, and implied volatility of calls and puts. Simply create an Inputs sturct and call a calc_<>() method on it.  
  
This branch is compilable to a python package using pyo3 and Maturin.  This package features full doc and type annotations. The rust compiled python package is ~1 second slower in pricing an option to 10M iterations than the pure rust crate.
