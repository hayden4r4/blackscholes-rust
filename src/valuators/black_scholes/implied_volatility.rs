use crate::{
    lets_be_rational::implied_volatility_from_a_transformed_rational_guess,
    valuators::black_scholes::Inputs,
    ImpliedVolatility,
};

impl ImpliedVolatility<f64> for Inputs {
    /// Calculates the implied volatility of the option.
    /// # Requires
    /// s, k, r, t, p
    /// # Returns
    /// f64 of the implied volatility of the option.
    /// # Example:
    /// ```
    /// use option_valuators::{valuators::black_scholes::Inputs, OptionType, ImpliedVolatility};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, Some(0.2), 0.05, 0.05, 20.0/365.25, None);
    /// let iv = inputs.calc_iv().unwrap();
    /// ```
    ///
    /// Uses the "Let's be rational" method from ["Let’s be rational" (2016) by Peter Jackel](http://www.jaeckel.org/LetsBeRational.pdf)
    /// Per Jackel's whitepaper, this method can solve for the implied volatility to f64 precision in 2 iterations.
    fn calc_iv(&self) -> Result<f64, String> {
        // extract price, or return error
        let p = self.p.ok_or("Option price is required".to_string())?;

        // "let's be rational" works with the forward and undiscounted option price, so remove the discount
        let rate_inv_discount = (self.r * self.t).exp();
        let p = p * rate_inv_discount;

        // compute the forward price
        let f = self.s * rate_inv_discount;
        // The Black-Scholes-Merton formula takes into account dividend yield by setting S = S * e^{-qt}, do this here with the forward
        let f = f * (-self.q * self.t).exp();

        let sigma = implied_volatility_from_a_transformed_rational_guess(
            p,
            f,
            self.k,
            self.t,
            self.option_type,
        );

        if sigma.is_nan() || sigma.is_infinite() || sigma < 0.0 {
            Err("Implied volatility failed to converge".to_string())?
        }
        Ok(sigma)
    }
}
