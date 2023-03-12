use crate::{constants::*, greeks::Greeks, pricing::Pricing, Inputs};
use num_traits::Float;
pub trait ImpliedVolatility<T>: Pricing<T> + Greeks<T>
where
    T: Float,
{
    fn calc_iv(&self, tolerance: T) -> Result<T, String>;
}

impl ImpliedVolatility<f32> for Inputs {
    /// Calculates the implied volatility of the option.
    /// Tolerance is the max error allowed for the implied volatility,
    /// the lower the tolerance the more iterations will be required.
    /// Recommended to be a value between 0.001 - 0.0001 for highest efficiency/accuracy.
    /// Initializes estimation of sigma using Brenn and Subrahmanyam (1998) method of calculating initial iv estimation.
    /// Uses Newton Raphson algorithm to calculate implied volatility.
    /// # Requires
    /// s, k, r, q, t, p
    /// # Returns
    /// f32 of the implied volatility of the option.
    /// # Example:
    /// ```
    /// use blackscholes::{Inputs, OptionType, ImpliedVolatility};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, Some(0.2), 0.05, 0.2, 20.0/365.25, None);
    /// let iv = inputs.calc_iv(0.0001).unwrap();
    /// ```
    /// Initial estimation of sigma using Modified Corrado-Miller from ["A MODIFIED CORRADO-MILLER IMPLIED VOLATILITY ESTIMATOR" (2007) by Piotr P√luciennik](https://sin.put.poznan.pl/files/download/37938) method of calculating initial iv estimation.
    /// Note: While this method is more accurate than Brenn and Subrahmanyam (1998) it will still sometimes fail to converge.
    /// An example of failure to converge:
    /// ```
    /// use blackscholes::{Inputs, OptionType, ImpliedVolatility};
    /// let inputs = Inputs::new(OptionType::Call, 105.0, 100.0, Some(30.0), 0.05, 30.0 / 365.25, None).calc_iv(0.0001).unwrap();
    /// // This will fail to converge, the NaN sigma value is checked in the function and will return an error.
    /// assert(inputs.calc_iv(0.0001).is_err(), true);
    /// ```
    ///
    /// A more accurate method is the "Let's be rational" method from ["Let’s be rational" (2016) by Peter Jackel](http://www.jaeckel.org/LetsBeRational.pdf)
    /// however this method is much more complicated, it is available as calc_rational_iv() in the RationalImpliedVolatility trait.
    #[allow(non_snake_case)]
    fn calc_iv(&self, tolerance: f32) -> Result<f32, String> {
        let mut inputs: Inputs = self.clone();

        let p = self
            .p
            .ok_or("inputs.p must contain Some(f32), found None".to_string())?;
        // Initialize estimation of sigma using Brenn and Subrahmanyam (1998) method of calculating initial iv estimation.
        // Note: I have since found that this method is not accurate enough which will sometimes cause the algorithm to fail to converge, hence I have commented it out.
        // let mut sigma: f32 = (PI2 / inputs.t).sqrt() * (p / inputs.s);

        let X: f32 = inputs.k * E.powf(-inputs.r * inputs.t);
        let fminusX: f32 = inputs.s - X;
        let fplusX: f32 = inputs.s + X;
        let oneoversqrtT: f32 = 1.0 / inputs.t.sqrt();

        let x: f32 = oneoversqrtT * (SQRT_2PI / (fplusX));
        let y: f32 = p - (inputs.s - inputs.k) / 2.0
            + ((p - fminusX / 2.0).powf(2.0) - fminusX.powf(2.0) / PI).sqrt();

        let mut sigma: f32 = oneoversqrtT
            * (SQRT_2PI / fplusX)
            * (p - fminusX / 2.0 + ((p - fminusX / 2.0).powf(2.0) - fminusX.powf(2.0) / PI).sqrt())
            + A
            + B / x
            + C * y
            + D / x.powf(2.0)
            + _E * y.powf(2.0)
            + F * y / x;

        if sigma.is_nan() {
            Err("Failed to converge".to_string())?
        }
        // Initialize diff to 100 for use in while loop
        let mut diff: f32 = 100.0;

        // Uses Newton Raphson algorithm to calculate implied volatility.
        // Test if the difference between calculated option price and actual option price is > tolerance,
        // if so then iterate until the difference is less than tolerance
        while diff.abs() > tolerance {
            inputs.sigma = Some(sigma);
            diff = Inputs::calc_price(&inputs)? - p;
            sigma -= diff / (Inputs::calc_vega(&inputs)? * 100.0);
        }
        Ok(sigma)
    }
}
