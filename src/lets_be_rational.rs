/// FFI wrapper for the `lets_be_rational` library.
use libc::c_double;

#[link(name = "liblets_be_rational")]
extern "C" {
    // rename the ffi functions with `_ffi` so the wrapper functions with native types can have the name from the library

    #[link_name = "implied_volatility_from_a_transformed_rational_guess"]
    fn implied_volatility_from_a_transformed_rational_guess_ffi(
        price: c_double,
        F: c_double,
        K: c_double,
        T: c_double,
        q: c_double,
    ) -> c_double;

    #[link_name = "black"]
    //  double K, double sigma, double T, double q /* q=±1 */) -> c_double
    fn black_ffi(
        F: c_double,
        K: c_double,
        sigma: c_double,
        T: c_double,
        q: c_double,
    ) -> c_double;
}

/// This function returns the implied volatility of an option contract using a transformed rational approximation.
/// The function is a wrapper around the C function `implied_volatility_from_a_transformed_rational_guess`
/// from the library `lets_be_rational`.
/// # Requires
/// price, f, k, t, q.
/// # Returns
/// f64 of the implied volatility of the option.
pub fn implied_volatility_from_a_transformed_rational_guess(
    price: f64,
    f: f64,
    k: f64,
    t: f64,
    q: f64,
) -> f64 {
    let price: c_double = price.into();
    let f: c_double = f.into();
    let k: c_double = k.into();
    let t: c_double = t.into();
    let q: c_double = q.into();
    unsafe { implied_volatility_from_a_transformed_rational_guess_ffi(price, f, k, t, q) }
}

/// This function returns the Black price of an option contract.
/// The function is a wrapper around the C function `black` from the library `lets_be_rational`.
/// # Requires
/// f, k, sigma, t, q.
/// # Returns
/// f64 of the price of the option.
pub fn black(f: f64, k: f64, sigma: f64, t: f64, q: f64) -> f64 {
    let f: c_double = f.into();
    let k: c_double = k.into();
    let sigma: c_double = sigma.into();
    let t: c_double = t.into();
    let q: c_double = q.into();
    unsafe { black_ffi(f, k, sigma, t, q) }
}
