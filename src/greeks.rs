use crate::Pricing;
use num_traits::Float;
use std::collections::HashMap;

pub trait Greeks<T>: Pricing<T>
where
    T: Float,
{
    fn calc_delta(&self) -> Result<T, String>;
    fn calc_gamma(&self) -> Result<T, String>;
    fn calc_theta(&self) -> Result<T, String>;
    fn calc_vega(&self) -> Result<T, String>;
    fn calc_rho(&self) -> Result<T, String>;
    fn calc_epsilon(&self) -> Result<T, String>;
    fn calc_lambda(&self) -> Result<T, String>;
    fn calc_vanna(&self) -> Result<T, String>;
    fn calc_charm(&self) -> Result<T, String>;
    fn calc_veta(&self) -> Result<T, String>;
    fn calc_vomma(&self) -> Result<T, String>;
    fn calc_speed(&self) -> Result<T, String>;
    fn calc_zomma(&self) -> Result<T, String>;
    fn calc_color(&self) -> Result<T, String>;
    fn calc_ultima(&self) -> Result<T, String>;
    fn calc_dual_delta(&self) -> Result<T, String>;
    fn calc_dual_gamma(&self) -> Result<T, String>;
    fn calc_all_greeks(&self) -> Result<HashMap<String, T>, String>;
}
