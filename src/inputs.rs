use std::{
    fmt::{Display, Formatter, Result as fmtResult},
    ops::Neg,
};

use num_traits::ConstZero;
use crate::calc_nd1nd2;
use crate::error::BlackScholesError;

/// The type of option to be priced (call or put).
#[derive(Debug, Clone, Eq, PartialEq, Copy)]
#[repr(i8)]
pub enum OptionType {
    Call = 1,
    Put = -1,
}

impl Neg for OptionType {
    type Output = Self;

    #[inline]
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

macro_rules! impl_option_type {
    ($type:ty) => {
        impl From<OptionType> for $type {
            #[inline]
            fn from(val: OptionType) -> Self {
                <$type>::from(val as i8)
            }
        }

        impl From<$type> for OptionType {
            #[inline]
            fn from(value: $type) -> Self {
                if value >= <$type>::ZERO {
                    OptionType::Call
                } else {
                    OptionType::Put
                }
            }
        }

        impl std::ops::Mul<OptionType> for $type {
            type Output = $type;

            #[inline]
            fn mul(self, rhs: OptionType) -> Self::Output {
                match rhs {
                    OptionType::Call => self,
                    OptionType::Put => -self,
                }
            }
        }

        impl std::ops::Mul<$type> for OptionType {
            type Output = $type;

            #[inline]
            fn mul(self, rhs: $type) -> Self::Output {
                match self {
                    OptionType::Call => rhs,
                    OptionType::Put => -rhs,
                }
            }
        }
    };
}

impl_option_type!(f32);
impl_option_type!(f64);
impl_option_type!(i8);
impl_option_type!(i16);
impl_option_type!(i32);
impl_option_type!(i64);
impl_option_type!(i128);
impl_option_type!(isize);

/// The inputs to the Black-Scholes-Merton model.
#[derive(Debug, Clone, PartialEq)]
pub struct Inputs {
    /// The type of the option (call or put)
    pub option_type: OptionType,
    /// Stock price
    pub s: f64,
    /// Strike price
    pub k: f64,
    /// Option price
    pub p: Option<f64>,
    /// Risk-free rate
    pub r: f64,
    /// Dividend yield
    pub q: f64,
    /// Time to maturity in years
    pub t: f64,
    /// Volatility
    pub sigma: Option<f64>,
}

/// Methods for calculating the price, greeks, and implied volatility of an option.
impl Inputs {
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
        s: f64,
        k: f64,
        p: Option<f64>,
        r: f64,
        q: f64,
        t: f64,
        sigma: Option<f64>,
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

    /// Calculates the price of the option.
    /// # Requires
    /// s, k, r, q, t, sigma.
    /// # Returns
    /// f64 of the price of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Pricing};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let price = inputs.calc_price()?;
    /// ```
    pub fn calc_price(&self) -> Result<f64, BlackScholesError> {
        // Calculates the price of the option
        let (nd1, nd2) = if let Ok((nd1, nd2)) = calc_nd1nd2(self) {
            (nd1, nd2)
        }
        else {
            // NOTE: Different perf behavior Enum vs String
            // Tt is not possible based on enum benches tests
            println!("Error in calc_nd1nd2");
            // return Err("ExpectedSigma".to_string());
            return Err(BlackScholesError::ExpectedSigma);
            // return Err(1);
        };

        let price: f64 = match self.option_type {
            OptionType::Call => f64::max(
                0.0,
                nd1 * self.s * (-self.q * self.t).exp() - nd2 * self.k * (-self.r * self.t).exp(),
            ),
            OptionType::Put => f64::max(
                0.0,
                nd2 * self.k * (-self.r * self.t).exp() - nd1 * self.s * (-self.q * self.t).exp(),
            ),
        };
        Ok(price)
    }
}

impl Display for Inputs {
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
