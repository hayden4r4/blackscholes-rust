use std::{
    fmt::{Display, Formatter, Result as fmtResult},
    ops::Neg,
};

use num_traits::ConstZero;

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