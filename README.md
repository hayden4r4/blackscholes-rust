# Black-Scholes Rust

A Black-Scholes pricing model built in Rust. Will calculate the price of calls or puts. The computation is parallelized using crossbeam. Simply pass an Inputs struct to the calc_price() function to return the calculated price of an option.
