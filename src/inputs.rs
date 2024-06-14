use std::fmt::{Display, Formatter, Result as fmtResult};
use std::ops::Neg;

/// The type of option to be priced (call or put).
#[derive(Debug, Clone, Eq, PartialEq, Copy)]
#[repr(i32)]
pub enum OptionType {
    Call = 1,
    Put = -1,
}

impl Neg for OptionType {
    type Output = Self;

    fn neg(self) -> Self::Output {
        match self {
            OptionType::Call => OptionType::Put,
            OptionType::Put => OptionType::Call,
        }
    }
}

impl Display for OptionType {
    fn fmt(&self, f: &mut Formatter) -> fmtResult {
        match self {
            OptionType::Call => write!(f, "Call"),
            OptionType::Put => write!(f, "Put"),
        }
    }
}

impl From<OptionType> for f64 {
    fn from(val: OptionType) -> Self {
        val as i32 as f64
    }
}

/// The inputs to the Black-Scholes-Merton model.
#[derive(Debug, Clone, PartialEq)]
pub struct Inputs<T> {
    /// The type of the option (call or put)
    pub option_type: OptionType,
    /// Stock price
    pub s: T,
    /// Strike price
    pub k: T,
    /// Option price
    pub p: Option<T>,
    /// Risk-free rate
    pub r: T,
    /// Dividend yield
    pub q: T,
    /// Time to maturity in years
    pub t: T,
    /// Volatility
    pub sigma: Option<T>,
}

/// Methods for calculating the price, greeks, and implied volatility of an option.
impl<T> Inputs<T> {
    /// Creates instance ot the `Inputs` struct.
    /// # Arguments
    /// * `option_type` - The type of option to be priced.
    /// * `s` - The current price of the underlying asset.
    /// * `k` - The strike price of the option.
    /// * `p` - The dividend yield of the underlying asset.
    /// * `r` - The risk-free interest rate.
    /// * `q` - The dividend yield of the underlying asset.
    /// * `t` - The time to maturity of the option in years.
    /// * `sigma` - The volatility of the underlying asset.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// ```
    /// # Returns
    /// An instance of the `Inputs` struct.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        option_type: OptionType,
        s: T,
        k: T,
        p: Option<T>,
        r: T,
        q: T,
        t: T,
        sigma: Option<T>,
    ) -> Self {
        Self {
            option_type,
            s,
            k,
            p,
            r,
            q,
            t,
            sigma,
        }
    }
}

impl<T> Display for Inputs<T>
where
    T: std::fmt::Display + Copy,
{
    fn fmt(&self, f: &mut Formatter) -> fmtResult {
        writeln!(f, "Option type: {}", self.option_type)?;
        writeln!(f, "Stock price: {:.2}", self.s)?;
        writeln!(f, "Strike price: {:.2}", self.k)?;
        match self.p {
            Some(p) => writeln!(f, "Option price: {:.2}", p)?,
            None => writeln!(f, "Option price: None")?,
        }
        writeln!(f, "Risk-free rate: {:.4}", self.r)?;
        writeln!(f, "Dividend yield: {:.4}", self.q)?;
        writeln!(f, "Time to maturity: {:.4}", self.t)?;
        match self.sigma {
            Some(sigma) => writeln!(f, "Volatility: {:.4}", sigma)?,
            None => writeln!(f, "Volatility: None")?,
        }
        Ok(())
    }
}
