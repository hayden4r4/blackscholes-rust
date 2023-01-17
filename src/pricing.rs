use crate::{common::*, constants::*, Inputs, OptionType};
use num_traits::Float;
pub trait Pricing<T>
where
    T: Float,
{
    fn calc_price(&self) -> Result<T, String>;
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
        // Calculates the price of the option
        let (nd1, nd2): (f32, f32) = calc_nd1nd2(&self)?;
        let price: f32 = match self.option_type {
            OptionType::Call => f32::max(
                0.0,
                nd1 * self.s * E.powf(-self.q * self.t) - nd2 * self.k * E.powf(-self.r * self.t),
            ),
            OptionType::Put => f32::max(
                0.0,
                nd2 * self.k * E.powf(-self.r * self.t) - nd1 * self.s * E.powf(-self.q * self.t),
            ),
        };
        Ok(price)
    }
}
