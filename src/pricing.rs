use std::f32::consts::E;

use num_traits::Float;

use crate::{lets_be_rational, Inputs, OptionType, *};

// Static error message to avoid allocation
static ERR_MISSING_SIGMA: &str = "Expected Some(f32) for self.sigma, received None";

pub trait Pricing<T>
where
    T: Float,
{
    fn calc_price(&self) -> Result<T, String>;
    fn calc_rational_price(&self) -> Result<f64, String>;
}

impl Pricing<f32> for Inputs {
    /// Calculates the price of the option.
    /// # Requires
    /// s, k, r, q, t, sigma.
    /// # Returns
    /// f32 of the price of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Pricing};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let price = inputs.calc_price().unwrap();
    /// ```
    fn calc_price(&self) -> Result<f32, String> {
        // Pre-calculate common exponential terms
        let discount_r = E.powf(-self.r * self.t);
        let discount_q = E.powf(-self.q * self.t);
        
        // Get normalized distribution values
        let (nd1, nd2): (f32, f32) = calc_nd1nd2(self)?;
        
        // Combine terms to minimize operations
        let term1 = nd1 * self.s * discount_q;
        let term2 = nd2 * self.k * discount_r;
        
        // Perform the option type-specific calculation
        let price: f32 = match self.option_type {
            OptionType::Call => f32::max(0.0, term1 - term2),
            OptionType::Put => f32::max(0.0, term2 - term1),
        };
        
        Ok(price)
    }

    /// Calculates the price of the option using the "Let's Be Rational" implementation.
    /// # Requires
    /// s, k, r, q, t, sigma.
    /// # Returns
    /// f64 of the price of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Pricing};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let price = inputs.calc_rational_price().unwrap();
    /// ```
    fn calc_rational_price(&self) -> Result<f64, String> {
        let sigma = self.sigma.ok_or(ERR_MISSING_SIGMA)?;

        // Pre-calculate the time component once - exp((r-q)*t)
        let time_factor = ((self.r - self.q) * self.t).exp();
        
        // Calculate forward price using the time factor
        let forward = self.s * time_factor;

        // Convert input parameters to f64 once, outside the function call
        let forward_f64 = forward as f64;
        let strike_f64 = self.k as f64;
        let sigma_f64 = sigma as f64;
        let time_f64 = self.t as f64;
        
        // Get undiscounted price
        let undiscounted_price = lets_be_rational::black(
            forward_f64,
            strike_f64,
            sigma_f64,
            time_f64,
            self.option_type,
        );

        // Pre-calculate discount factor as f64 to avoid type conversion in multiplication
        let discount_factor = (-self.r as f64 * time_f64).exp();
        
        // Apply discount 
        let price = undiscounted_price * discount_factor;
        Ok(price)
    }
}
