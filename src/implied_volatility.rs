use crate::{Greeks, Pricing};
use num_traits::Float;

pub trait ImpliedVolatility<T>: Pricing<T> + Greeks<T>
where
    T: Float,
{
    fn calc_iv(&self) -> Result<T, String>;
}
