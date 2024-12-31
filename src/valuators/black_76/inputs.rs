use crate::{OptionType, Shift};
use std::fmt::{Display, Formatter, Result as fmtResult};

/// The inputs to the Black-Scholes-Merton model.
#[derive(Debug, Clone, PartialEq)]
pub struct Inputs {
    /// The type of the option (call or put)
    pub option_type: OptionType,
    /// Forward price
    pub f: f64,
    /// Strike price
    pub k: f64,
    /// Option price
    pub p: Option<f64>,
    /// Risk-free rate
    pub r: f64,
    /// Time to maturity in years
    pub t: f64,
    /// Volatility
    pub sigma: Option<f64>,
    /// Whether to use the shifted Black76 model (in the case of negative futures/strike prices or rates)
    pub shifted: bool,
}

/// Methods for calculating the price, greeks, and implied volatility of an option.
impl Inputs {
    /// Creates instance ot the `Inputs` struct.
    /// # Arguments
    /// * `option_type` - The type of option to be priced.
    /// * `f` - The current price of the underlying asset.
    /// * `k` - The strike price of the option.
    /// * `p` - The dividend yield of the underlying asset.
    /// * `r` - The risk-free interest rate.
    /// * `t` - The time to maturity of the option in years.
    /// * `sigma` - The volatility of the underlying asset.
    /// * `shifted` - Whether to use the shifted Black76 model.
    /// # Example
    /// ```
    /// use blackscholes::OptionType;
    /// use blackscholes::valuators::black_76::Inputs;
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 20.0/365.25, Some(0.2), true);
    /// ```
    /// # Returns
    /// An instance of the `Inputs` struct.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        option_type: OptionType,
        f: f64,
        k: f64,
        p: Option<f64>,
        r: f64,
        t: f64,
        sigma: Option<f64>,
        shifted: bool,
    ) -> Self {
        Self {
            option_type,
            f,
            k,
            p,
            r,
            t,
            sigma,
            shifted
        }
    }
}

impl Shift<f64> for Inputs {
    // Shifts the forward and strike prices to avoid negative values.
    // Useful in the case of negative futures/strike prices or rates.
    // The shift amount is determined dynamically based on the minimum value of the forward price, strike price, and risk-free rate.
    fn shift(&self) -> (f64, f64) {
        let min = self.f.min(self.k).min(self.r);
        // Arbitrary shift factor
        let shift_factor = 0.01;
        let shift = shift_factor * min.max(0.0).abs();

        (self.f + shift, self.k + shift)
    }
}

impl Display for Inputs {
    fn fmt(&self, f: &mut Formatter) -> fmtResult {
        writeln!(f, "Option type: {}", self.option_type)?;
        writeln!(f, "Forward price: {:.2}", self.f)?;
        writeln!(f, "Strike price: {:.2}", self.k)?;
        match self.p {
            Some(p) => writeln!(f, "Option price: {:.2}", p)?,
            None => writeln!(f, "Option price: None")?,
        }
        writeln!(f, "Risk-free rate: {:.4}", self.r)?;
        writeln!(f, "Time to maturity: {:.4}", self.t)?;
        match self.sigma {
            Some(sigma) => writeln!(f, "Volatility: {:.4}", sigma)?,
            None => writeln!(f, "Volatility: None")?,
        }
        writeln!(f, "Shifted: {}", self.shifted)?;
        Ok(())
    }
}
