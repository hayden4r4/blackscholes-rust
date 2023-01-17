use crate::{constants::*, greeks::Greeks, pricing::Pricing, Inputs};
use num_traits::Float;
pub trait ImpliedVolatility<T>
where
    T: Float,
{
    fn calc_iv(&self, tolerance: T) -> Result<T, String>;
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
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, Some(0.2), 0.05, 0.2, 20.0/365.25, None);
    /// let iv = inputs.calc_iv(0.0001).unwrap();
    /// ```
    fn calc_iv(&self, tolerance: f32) -> Result<f32, String> {
        let mut inputs: Inputs = self.clone();

        let p = self
            .p
            .ok_or("inputs.p must contain Some(f32), found None".to_string())?;
        // Initialize estimation of sigma using Brenn and Subrahmanyam (1998) method of calculating initial iv estimation.
        let mut sigma: f32 = (2.0 * PI / inputs.t).sqrt() * (p / inputs.s);
        // Initialize diff to 100 for use in while loop
        let mut diff: f32 = 100.0;

        // Uses Newton Raphson algorithm to calculate implied volatility.
        // Test if the difference between calculated option price and actual option price is > tolerance,
        // if so then iterate until the difference is less than tolerance
        while diff.abs() > tolerance {
            inputs.sigma = Some(sigma);
            diff = Inputs::calc_price(&inputs)? - p;
            sigma -= diff / (Inputs::calc_vega(&inputs)? * 100.0);
        }
        Ok(sigma)
    }
}
