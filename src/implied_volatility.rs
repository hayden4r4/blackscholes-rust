use num_traits::{float::FloatConst, AsPrimitive, Float, FromPrimitive};

use crate::{
    greeks::Greeks, lets_be_rational::implied_volatility_from_a_transformed_rational_guess,
    pricing::Pricing, Inputs, A, B, C, D, F, SQRT_2PI, _E,
};

pub trait ImpliedVolatility<T>: Pricing<T> + Greeks<T>
where
    T: Float,
{
    fn calc_iv(&self, tolerance: T) -> Result<T, String>;
    fn calc_rational_iv(&self) -> Result<T, String>;
}

impl<T> ImpliedVolatility<T> for Inputs<T>
where
    T: Float + FromPrimitive + AsPrimitive<f64> + FloatConst,
{
    /// Calculates the implied volatility of the option.
    /// Tolerance is the max error allowed for the implied volatility,
    /// the lower the tolerance the more iterations will be required.
    /// Recommended to be a value between 0.001 - 0.0001 for highest efficiency/accuracy.
    /// Initializes estimation of sigma using Brenn and Subrahmanyam (1998) method of calculating initial iv estimation.
    /// Uses Newton Raphson algorithm to calculate implied volatility.
    /// # Requires
    /// s, k, r, q, t, p
    /// # Returns
    /// T of the implied volatility of the option.
    /// # Example:
    /// ```
    /// use blackscholes::{Inputs, OptionType, ImpliedVolatility};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, Some(0.5), 0.05, 0.2, 20.0/365.25, None);
    /// let iv = inputs.calc_iv(0.0001).unwrap();
    /// ```
    /// Initial estimation of sigma using Modified Corrado-Miller from ["A MODIFIED CORRADO-MILLER IMPLIED VOLATILITY ESTIMATOR" (2007) by Piotr P√luciennik](https://sin.put.poznan.pl/files/download/37938) method of calculating initial iv estimation.
    /// A more accurate method is the "Let's be rational" method from ["Let’s be rational" (2016) by Peter Jackel](http://www.jaeckel.org/LetsBeRational.pdf)
    /// however this method is much more complicated, it is available as calc_rational_iv().
    #[allow(non_snake_case)]
    fn calc_iv(&self, tolerance: T) -> Result<T, String> {
        let mut inputs: Inputs<T> = self.clone();

        let p = self
            .p
            .ok_or("inputs.p must contain Some(T), found None".to_string())?;

        let half = T::from(0.5).unwrap();

        let X: T = inputs.k * T::E().powf(-inputs.r * inputs.t);
        let fminusX: T = inputs.s - X;
        let fplusX: T = inputs.s + X;
        let oneoversqrtT: T = T::one() / inputs.t.sqrt();

        let x: T = oneoversqrtT * (T::from(SQRT_2PI).unwrap() / (fplusX));
        let y: T = p - (inputs.s - inputs.k) * half
            + ((p - fminusX * half).powi(2) - fminusX.powi(2) / T::PI()).sqrt();

        let mut sigma: T = oneoversqrtT
            * (T::from(SQRT_2PI).unwrap() / fplusX)
            * (p - fminusX * half
                + ((p - fminusX * half).powi(2) - fminusX.powi(2) / T::PI()).sqrt())
            + T::from(A).unwrap()
            + T::from(B).unwrap() / x
            + T::from(C).unwrap() * y
            + T::from(D).unwrap() / x.powi(2)
            + T::from(_E).unwrap() * y.powi(2)
            + T::from(F).unwrap() * y / x;

        if sigma.is_nan() {
            Err("Failed to converge".to_string())?
        }

        // Initialize diff to 100 for use in while loop
        let mut diff: T = T::from(100.0).unwrap();

        // Uses Newton Raphson algorithm to calculate implied volatility.
        // Test if the difference between calculated option price and actual option price is > tolerance,
        // if so then iterate until the difference is less than tolerance
        while diff.abs() > tolerance {
            inputs.sigma = Some(sigma);
            diff = Inputs::calc_price(&inputs)? - p;
            sigma = sigma - diff / (Inputs::calc_vega(&inputs)? * T::from(100.0).unwrap());

            if sigma.is_nan() || sigma.is_infinite() {
                Err("Failed to converge".to_string())?
            }
        }
        Ok(sigma)
    }

    /// Calculates the implied volatility of the option.
    /// # Requires
    /// s, k, r, t, p
    /// # Returns
    /// f64 of the implied volatility of the option.
    /// # Example:
    /// ```
    /// use blackscholes::{Inputs, OptionType, ImpliedVolatility};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, Some(0.2), 0.05, 0.05, 20.0/365.25, None);
    /// let iv = inputs.calc_rational_iv().unwrap();
    /// ```
    ///
    /// Uses the "Let's be rational" method from ["Let’s be rational" (2016) by Peter Jackel](http://www.jaeckel.org/LetsBeRational.pdf)
    /// from Jackel's C++ implementation, imported through the C FFI.  The C++ implementation is available at [here](http://www.jaeckel.org/LetsBeRational.7z)
    /// Per Jackel's whitepaper, this method can solve for the implied volatility to f64 precision in 2 iterations.
    fn calc_rational_iv(&self) -> Result<T, String> {
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

        if sigma.is_nan() || sigma.is_infinite() || sigma < T::zero() {
            Err("Implied volatility failed to converge".to_string())?
        }
        Ok(sigma)
    }
}
