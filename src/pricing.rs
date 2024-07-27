use num_traits::Float;

use crate::{lets_be_rational, Inputs, OptionType, *};

pub trait Pricing<T>
where
    T: Float,
{
    fn calc_price(&self) -> Result<T, String>;
    fn calc_rational_price(&self) -> Result<f64, String>;
}

impl Pricing<f64> for Inputs {
    /// Calculates the price of the option.
    /// # Requires
    /// s, k, r, q, t, sigma.
    /// # Returns
    /// f64 of the price of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Pricing};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let price = inputs.calc_price().unwrap();
    /// ```
    fn calc_price(&self) -> Result<f64, String> {
        // Calculates the price of the option
        let (nd1, nd2): (f64, f64) = calc_nd1nd2(self)?;
        let price: f64 = match self.option_type {
            OptionType::Call => f64::max(
                0.0,
                nd1 * self.s * (-self.q * self.t).exp() - nd2 * self.k * (-self.r * self.t).exp(),
            ),
            OptionType::Put => f64::max(
                0.0,
                nd2 * self.k * (-self.r * self.t).exp() - nd1 * self.s * (-self.q * self.t).exp(),
            ),
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
        let sigma = self
            .sigma
            .ok_or("Expected Some(f64) for self.sigma, received None")?;

        // let's be rational wants the forward price, not the spot price.
        let forward = self.s * ((self.r - self.q) * self.t).exp();

        // price using `black`
        let undiscounted_price =
            lets_be_rational::black(forward, self.k, sigma, self.t, self.option_type);

        // discount the price
        let price = undiscounted_price * (-self.r * self.t).exp();
        Ok(price)
    }
}
