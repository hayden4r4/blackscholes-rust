# blackscholes_python  
  
This library provides an simple, lightweight, and efficient (though not heavily optimized) implementation of the Black-Scholes-Merton model for pricing European options.  
  
This crate is compilable to a python package using pyo3 and Maturin.  It features full doc and type annotations. The rust compiled python package is ~1 second slower in pricing an option to 10M iterations than the pure rust crate on a Ryzen 7600x.
  
## Usage  
  
Simply create an instance of the `Inputs` struct and call the desired method.  
  
View the [Rust docs](https://docs.rs/blackscholes_python) and [Python docs](https://pypi.org/project/blackscholes-python/) for usage and examples.  
