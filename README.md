# Black-Scholes Rust

A Black-Scholes pricing model built in Rust. Will calculate the price of calls or puts. Simply pass an Inputs struct to the Price::calc_price(), Greek::calc_\<greek>(), or Volatility::calc_iv() function to return the calculated price, desired greek, or implied volatility of an option.  

This branch is compilable to a python package using pyo3 and Maturin.  This package features full doc and type annotations. The rust compiled python package is ~1 second slower in pricing an option to 10M iterations than the pure rust crate.
