use num_traits::Float;
use once_cell::sync::Lazy;

use crate::lets_be_rational::implied_volatility_from_a_transformed_rational_guess;
use crate::{greeks::Greeks, pricing::Pricing, Inputs, *};

// Static error messages to avoid string allocations
static ERR_MISSING_PRICE: Lazy<String> = Lazy::new(|| "inputs.p must contain Some(f32), found None".to_string());
static ERR_FAILED_CONVERGE: Lazy<String> = Lazy::new(|| "Failed to converge".to_string());
static ERR_IV_FAILED: Lazy<String> = Lazy::new(|| "Implied volatility failed to converge".to_string());
static ERR_OPTION_PRICE_REQUIRED: Lazy<String> = Lazy::new(|| "Option price is required".to_string());

pub trait ImpliedVolatility<T>: Pricing<T> + Greeks<T>
where
    T: Float,
{
    fn calc_iv(&self, tolerance: T) -> Result<T, String>;
    fn calc_rational_iv(&self) -> Result<f64, String>;
}

impl ImpliedVolatility<f32> for Inputs {
    /// Calculates the implied volatility of the option.
    /// Tolerance is the max error allowed for the implied volatility,
    /// the lower the tolerance the more iterations will be required.
    /// Recommended to be a value between 0.001 - 0.0001 for highest efficiency/accuracy.
    /// Initializes estimation of sigma using Brenn and Subrahmanyam (1998) method of calculating initial iv estimation.
    /// Uses Newton Raphson algorithm to calculate implied volatility.
    /// # Requires
    /// s, k, r, q, t, p
    /// # Returns
    /// f32 of the implied volatility of the option.
    /// # Example:
    /// ```
    /// use blackscholes::{Inputs, OptionType, ImpliedVolatility};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, Some(0.5), 0.05, 0.2, 20.0/365.25, None);
    /// let iv = inputs.calc_iv(0.0001).unwrap();
    /// ```
    /// Initial estimation of sigma using Modified Corrado-Miller from ["A MODIFIED CORRADO-MILLER IMPLIED VOLATILITY ESTIMATOR" (2007) by Piotr P√luciennik](https://sin.put.poznan.pl/files/download/37938) method of calculating initial iv estimation.
    /// A more accurate method is the "Let's be rational" method from ["Let's be rational" (2016) by Peter Jackel](http://www.jaeckel.org/LetsBeRational.pdf)
    /// however this method is much more complicated, it is available as calc_rational_iv().
    #[allow(non_snake_case)]
    fn calc_iv(&self, tolerance: f32) -> Result<f32, String> {
        // Avoid cloning the entire input struct by copying only the necessary fields
        let p = self.p.ok_or_else(|| ERR_MISSING_PRICE.clone())?;
        
        // Initialize estimation of sigma using Modified Corrado-Miller method
        let X: f32 = self.k * E.powf(-self.r * self.t);
        let fminusX: f32 = self.s - X;
        let fplusX: f32 = self.s + X;
        let oneoversqrtT: f32 = 1.0 / self.t.sqrt();

        let x: f32 = oneoversqrtT * (SQRT_2PI / (fplusX));
        let y: f32 = p - (self.s - self.k) / 2.0
            + ((p - fminusX / 2.0).powf(2.0) - fminusX.powf(2.0) / PI).sqrt();

        let mut sigma: f32 = oneoversqrtT
            * (SQRT_2PI / fplusX)
            * (p - fminusX / 2.0 + ((p - fminusX / 2.0).powf(2.0) - fminusX.powf(2.0) / PI).sqrt())
            + A
            + B / x
            + C * y
            + D / x.powf(2.0)
            + _E * y.powf(2.0)
            + F * y / x;

        if sigma.is_nan() {
            return Err(ERR_FAILED_CONVERGE.clone());
        }

        // Initialize diff to a value greater than tolerance
        let mut diff: f32 = tolerance * 2.0;
        
        // Create a mutable copy of inputs only when needed
        let mut inputs_copy = Inputs {
            option_type: self.option_type,
            s: self.s,
            k: self.k,
            p: self.p,
            r: self.r,
            q: self.q,
            t: self.t,
            sigma: Some(sigma),
        };

        // Uses Newton Raphson algorithm to calculate implied volatility.
        // Test if the difference between calculated option price and actual option price is > tolerance,
        // if so then iterate until the difference is less than tolerance
        while diff.abs() > tolerance {
            inputs_copy.sigma = Some(sigma);
            diff = Inputs::calc_price(&inputs_copy)? - p;
            sigma -= diff / (Inputs::calc_vega(&inputs_copy)? * 100.0);

            if sigma.is_nan() || sigma.is_infinite() {
                return Err(ERR_FAILED_CONVERGE.clone());
            }
        }
        Ok(sigma)
    }

    /// Calculates the implied volatility of the option.
    /// # Requires
    /// s, k, r, t, p
    /// # Returns
    /// f64 of the implied volatility of the option.
    /// # Example:
    /// ```
    /// use blackscholes::{Inputs, OptionType, ImpliedVolatility};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, Some(0.2), 0.05, 0.05, 20.0/365.25, None);
    /// let iv = inputs.calc_rational_iv().unwrap();
    /// ```
    ///
    /// Uses the "Let's be rational" method from ["Let's be rational" (2016) by Peter Jackel](http://www.jaeckel.org/LetsBeRational.pdf)
    /// from Jackel's C++ implementation, imported through the C FFI.  The C++ implementation is available at [here](http://www.jaeckel.org/LetsBeRational.7z)
    /// Per Jackel's whitepaper, this method can solve for the implied volatility to f64 precision in 2 iterations.
    fn calc_rational_iv(&self) -> Result<f64, String> {
        // extract price, or return error
        let p = self.p.ok_or_else(|| ERR_OPTION_PRICE_REQUIRED.clone())?;

        // "let's be rational" works with the forward and undiscounted option price, so remove the discount
        let rate_inv_discount = (self.r * self.t).exp();
        let p = p * rate_inv_discount;

        // compute the forward price
        let f = self.s * rate_inv_discount;
        // The Black-Scholes-Merton formula takes into account dividend yield by setting S = S * e^{-qt}, do this here with the forward
        let f = f * (-self.q * self.t).exp();

        let sigma = implied_volatility_from_a_transformed_rational_guess(
            p as f64,
            f as f64,
            self.k as f64,
            self.t as f64,
            self.option_type,
        );

        if sigma.is_nan() || sigma.is_infinite() || sigma < 0.0 {
            return Err(ERR_IV_FAILED.clone());
        }
        Ok(sigma)
    }
}
