use num_traits::Float;

use crate::{error::BlackScholesError, lets_be_rational, Inputs, OptionType, *};

pub trait Pricing<T>
where
    T: Float,
{
    fn calc_price(&self) -> Result<T, BlackScholesError>;
    fn calc_rational_price(&self) -> Result<f64, BlackScholesError>;
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
    fn calc_price(&self) -> Result<f64, BlackScholesError> {
        // Calculates the price of the option
        let (nd1, nd2): (f64, f64) = calc_nd1nd2(self)?;

        Ok(self.calc_price_with_params(
            nd1,
            nd2,
            self.s,
            self.q,
            self.t,
            self.k,
            self.r,
            self.option_type,
        ))
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
    fn calc_rational_price(&self) -> Result<f64, BlackScholesError> {
        let sigma = self.sigma.ok_or(BlackScholesError::ExpectedSigma)?;

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
