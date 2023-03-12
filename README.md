# blackscholes  
  
This library provides an simple, lightweight, and efficient (though not heavily optimized) implementation of the Black-Scholes-Merton model for pricing European options.  
  
Includes all first, second, and third order Greeks.  

- calc_iv() in the ImpliedVolatility trait which uses [Modified Corrado-Miller by Piotr P√luciennik (2007)](https://sin.put.poznan.pl/files/download/37938) for the initial volatility guess and the Newton Raphson algorithm to solve for the implied volatility.
- calc_rational_iv() in the RationalImpliedVolatility trait which uses "Let's be rational" method from ["Let’s be rational" (2016) by Peter Jackel](http://www.jaeckel.org/LetsBeRational.pdf).  Utilizing Jackel's C++ implementation to get convergence within 2 iterations with 64-bit floating point accuracy.
  
## Usage  
  
View the [docs](https://docs.rs/blackscholes) for usage and examples.  
  
**Other packages available:**  
Python: [Pypi](https://pypi.org/project/blackscholes-python/)  
WASM: [npm](https://www.npmjs.com/package/@haydenr4/blackscholes_wasm)  
