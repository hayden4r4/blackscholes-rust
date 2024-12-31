use num_traits::Float;

pub trait Shift<T>
where
    T: Float,
{
    fn shift(&self) -> (T, T);
}
