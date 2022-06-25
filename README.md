# Black-Scholes Rust

A Black-Scholes pricing model built in Rust. Will calculate the price of calls or puts. The computation is parallelized using crossbeam. Simply pass an Inputs struct to the calc_price(), calc_\<greek>(), or calc_iv() function to return the calculated price, desired greek, or implied volatility of an option.
