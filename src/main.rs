use polars::prelude::*;
use std::fs::File;

const STOCK_PRICE: u32 = 50;
const EXERCISE_PRICE: u32 = 100;
// Risk-free interest rate
const R: f32 = 0.02;
// Time to expiration in days
const TTE: u32 = 180;
// Standard deviation of log returns (volatility)
const STD: f32 = 0.3;
// Option price
const OPTION_PRICE: f32 = 0.5;

fn load_returns() -> Result<DataFrame> {
    let file = File::open("SP500_returns.csv")?;

    CsvReader::new(file)
        .infer_schema(None)
        .has_header(true)
        .finish()
}

fn log_return(returns_df: DataFrame) -> f32 {
    let df: DataFrame = returns_df;
    let close_series:Series = df["Close"].clone();
    let log_return: f32 = close_series.str_value(0).to_string().parse::<f32>().unwrap()/close_series.str_value(close_series.len()-1).to_string().parse::<f32>().unwrap().ln() / close_series.len() as f32;
    log_return
}

fn main() {
    let returns_df: DataFrame = load_returns().unwrap();
    println!("{}", log_return(returns_df));
}
