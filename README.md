# Black-Scholes Rust

A Black-Scholes pricing model built in Rust. Will calculate the price of calls or puts. The computation is parallelized using crossbeam. Simply pass an Inputs struct to the calc_price() or calc_<greek>() function to return the calculated price or desired greek of an option.
