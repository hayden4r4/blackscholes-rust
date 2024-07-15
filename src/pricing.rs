use num_traits::{float::FloatConst, AsPrimitive, Float};

use crate::{lets_be_rational, Inputs, OptionType, *};

pub trait Pricing<T>
where
    T: Float,
{
    fn calc_price(&self) -> Result<T, String>;
    fn calc_rational_price(&self) -> Result<f64, String>;
}

impl<T> Pricing<T> for Inputs<T>
where
    T: Float + FromPrimitive + AsPrimitive<f64> + FloatConst,
{
    /// Calculates the price of the option.
    /// # Requires
    /// s, k, r, q, t, sigma.
    /// # Returns
    /// T of the price of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Pricing};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let price = inputs.calc_price().unwrap();
    /// ```
    fn calc_price(&self) -> Result<T, String> {
        // Calculates the price of the option
        let (nd1, nd2): (T, T) = calc_nd1nd2(self)?;
        let price: T = match self.option_type {
            OptionType::Call => T::max(
                T::zero(),
                nd1 * self.s * T::E().powf(-self.q * self.t)
                    - nd2 * self.k * T::E().powf(-self.r * self.t),
            ),
            OptionType::Put => T::max(
                T::zero(),
                nd2 * self.k * T::E().powf(-self.r * self.t)
                    - nd1 * self.s * T::E().powf(-self.q * self.t),
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
            .ok_or("Expected Some(T) for self.sigma, received None")?;

        // let's be rational wants the forward price, not the spot price.
        let forward = self.s * ((self.r - self.q) * self.t).exp();

        // price using `black`
        let undiscounted_price = lets_be_rational::black(
            forward.as_(),
            self.k.as_(),
            sigma.as_(),
            self.t.as_(),
            self.option_type,
        );

        // discount the price
        let price = undiscounted_price * (-self.r.as_() * self.t.as_()).exp();
        Ok(price)
    }
}
