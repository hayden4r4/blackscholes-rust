/*
Notice from Jackel's C++ implementation:
Permission to use, copy, modify, and distribute this software is freely granted,
provided that this notice is preserved.

WARRANTY DISCLAIMER
The Software is provided "as is" without warranty of any kind, either express or implied,
including without limitation any implied warranties of condition, uninterrupted use,
merchantability, fitness for a particular purpose, or non-infringement.
 */
use crate::{Inputs, OptionType};
use libc::c_double;

#[link(name = "liblets_be_rational")]
extern "C" {
    fn implied_volatility_from_a_transformed_rational_guess(
        price: c_double,
        F: c_double,
        K: c_double,
        T: c_double,
        q: c_double,
    ) -> c_double;
}

pub trait RationalImpliedVolatility {
    fn calc_rational_iv(&self) -> Result<f64, String>;
}

/// Calculates the implied volatility of the option.
/// Tolerance is the max error allowed for the implied volatility,
/// the lower the tolerance the more iterations will be required.
/// Recommended to be a value between 0.001 - 0.0001 for highest efficiency/accuracy.
/// Initializes estimation of sigma using Brenn and Subrahmanyam (1998) method of calculating initial iv estimation.
/// Uses Newton Raphson algorithm to calculate implied volatility.
/// # Requires
/// s, k, r, t, p
/// # Returns
/// f32 of the implied volatility of the option.
/// # Example:
/// ```
/// use blackscholes::{Inputs, OptionType, RationalImpliedVolatility};
/// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, Some(0.2), 0.05, 20.0/365.25, None);
/// let iv = inputs.calc_rational_iv().unwrap();
/// ```
///
/// Uses the "Let's be rational" method from ["Let’s be rational" (2016) by Peter Jackel](http://www.jaeckel.org/LetsBeRational.pdf)
/// from Jackel's C++ implementation, imported through the C FFI.  The C++ implementation is available at [here](http://www.jaeckel.org/LetsBeRational.7z)
/// Per Jackel's whitepaper, this method can solve for the implied volatility to f64 precision in 2 iterations.
impl RationalImpliedVolatility for Inputs {
    fn calc_rational_iv(&self) -> Result<f64, String> {
        let p: c_double = match self.p {
            Some(p) => p.into(),
            None => return Err("Option price is required".to_string()),
        };
        let s: c_double = self.s.into();
        let k: c_double = self.k.into();
        let t: c_double = self.t.into();
        let q: c_double = match self.option_type {
            OptionType::Call => 1.0,
            OptionType::Put => -1.0,
        };
        let sigma = unsafe { implied_volatility_from_a_transformed_rational_guess(p, s, k, t, q) };

        if sigma.is_nan() || sigma.is_infinite() || sigma < 0.0 {
            Err("Implied volatility failed to converge".to_string())?
        }
        Ok(sigma)
    }
}
