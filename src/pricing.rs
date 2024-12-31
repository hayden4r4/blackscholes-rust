use num_traits::Float;

pub trait Pricing<T>
where
    T: Float,
{
    fn calc_price(&self) -> Result<T, String>;
}
