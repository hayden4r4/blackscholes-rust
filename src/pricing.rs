use num_traits::Float;

use crate::{*, Inputs, lets_be_rational, OptionType};
use crate::error::BlackScholesError;

pub trait Pricing
{
    // fn calc_price(&self) -> Result<f64, String>;
    fn calc_rational_price(&self) -> Result<f64, BlackScholesError>;
}

impl Pricing for Inputs {
    // /// Calculates the price of the option.
    // /// # Requires
    // /// s, k, r, q, t, sigma.
    // /// # Returns
    // /// f64 of the price of the option.
    // /// # Example
    // /// ```
    // /// use blackscholes::{Inputs, OptionType, Pricing};
    // /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    // /// let price = inputs.calc_price()?;
    // /// ```
    // fn calc_price(&self) -> Result<f64, String> {
    //     // Calculates the price of the option
    //     let (nd1, nd2) = if let Ok((nd1, nd2)) = calc_nd1nd2(self) {
    //         (nd1, nd2)
    //     }
    //     else {
    //         println!("Error in calc_nd1nd2");
    //         return Err("Dupa".to_string());
    //         // return Err(BlackScholesError::ExpectedSigma);
    //     };
    //
    //     let price: f64 = match self.option_type {
    //         OptionType::Call => f64::max(
    //             0.0,
    //             nd1 * self.s * (-self.q * self.t).exp() - nd2 * self.k * (-self.r * self.t).exp(),
    //         ),
    //         OptionType::Put => f64::max(
    //             0.0,
    //             nd2 * self.k * (-self.r * self.t).exp() - nd1 * self.s * (-self.q * self.t).exp(),
    //         ),
    //     };
    //     Ok(price)
    // }

    /// Calculates the price of the option using the "Let's Be Rational" implementation.
    /// # Requires
    /// s, k, r, q, t, sigma.
    /// # Returns
    /// f64 of the price of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Pricing};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let price = inputs.calc_rational_price()?;
    /// ```
    fn calc_rational_price(&self) -> Result<f64, BlackScholesError> {
        let sigma = if let Some(sigma) = self.sigma {
            sigma
        } else { return Err(BlackScholesError::ExpectedSigma); };

        // let sigma = self.sigma.ok_or("0".to_string()).map_err(|_| BlackScholesError::ExpectedSigma)?;

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
