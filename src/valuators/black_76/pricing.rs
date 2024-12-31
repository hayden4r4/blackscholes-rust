use crate::{
    lets_be_rational::black,
    valuators::black_76::{distributions::*, Inputs},
    OptionType, Pricing, Shift,
};

impl Pricing<f64> for Inputs {
    /// Calculates the price of the option.
    /// # Requires
    /// s, k, r, q, t, sigma.
    /// # Returns
    /// f64 of the price of the option.
    /// # Example
    /// ```
    /// use blackscholes::{OptionType, Pricing};
    /// use blackscholes::valuators::black_76::Inputs;
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 20.0/365.25, Some(0.2), true);
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

impl Inputs {
    /// Calculates the price of the option using the "Let's Be Rational" implementation.
    /// # Requires
    /// s, k, r, q, t, sigma.
    /// # Returns
    /// f64 of the price of the option.
    /// # Example
    /// ```
    /// use blackscholes::{OptionType, Pricing};
    /// use blackscholes::valuators::black_76::Inputs;
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 20.0/365.25, Some(0.2), true);
    /// let price = inputs.calc_rational_price().unwrap();
    /// ```
    pub fn calc_rational_price(&self) -> Result<f64, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f64) for self.sigma, received None")?;

        // price using `black`
        let undiscounted_price = black(self.f, self.k, sigma, self.t, self.option_type);

        // discount the price
        let price = undiscounted_price * (-self.r * self.t).exp();
        Ok(price)
    }
}
