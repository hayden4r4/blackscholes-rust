use crate::{
    OptionType,
    Pricing,
    Shift,
    valuators::black_76::{distributions::*, Inputs},
};

impl Pricing<f64> for Inputs {
    /// Calculates the price of the option.
    /// # Requires
    /// s, k, r, q, t, sigma.
    /// # Returns
    /// f64 of the price of the option.
    /// # Example
    /// ```
    /// use blackscholes::OptionType;
    /// use blackscholes::valuators::black_76::{Inputs, Pricing};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2), true);
    /// let price = inputs.calc_price().unwrap();
    /// ```
    fn calc_price(&self) -> Result<f64, String> {
        // Calculates the price of the option
        let (nd1, nd2): (f64, f64) = calc_nd1nd2(self)?;
        let (f, k) = if self.shifted {
            self.shift()
        } else {
            (self.f, self.k)
        };
        let price: f64 = match self.option_type {
            OptionType::Call => f64::max(
                0.0,
                nd1 * f * (-self.r * self.t).exp() - nd2 * k * (-self.r * self.t).exp(),
            ),
            OptionType::Put => f64::max(
                0.0,
                nd2 * k * (-self.r * self.t).exp() - nd1 * f * (-self.r * self.t).exp(),
            ),
        };
        Ok(price)
    }
}
