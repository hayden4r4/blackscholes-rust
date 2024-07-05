/// FFI wrapper for the `lets_be_rational` library.
use libc::c_double;
use log::error;

use crate::OptionType;

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
    fn black_ffi(F: c_double, K: c_double, sigma: c_double, T: c_double, q: c_double) -> c_double;
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
    q: OptionType,
) -> f64 {
    let price: c_double = price.into();
    let f: c_double = f.into();
    let k: c_double = k.into();
    let t: c_double = t.into();
    let q: c_double = match q {
        OptionType::Call => 1.0,
        OptionType::Put => -1.0,
    };
    unsafe {
        let result = implied_volatility_from_a_transformed_rational_guess_ffi(price, f, k, t, q);

        if result.is_nan() {
            error!("Implied volatility has failed in the calculation - NaN");
            0.0
        } else if result.is_sign_negative() {
            error!("Implied volatility has failed in the calculation - Negative value");
            result.into()
        } else {
            result.into()
        }
    }
}

/// This function returns the Black price of an option contract.
/// The function is a wrapper around the C function `black` from the library `lets_be_rational`.
/// # Requires
/// f, k, sigma, t, q.
/// # Returns
/// f64 of the price of the option.
pub fn black(f: f64, k: f64, sigma: f64, t: f64, q: OptionType) -> f64 {
    let f: c_double = f.into();
    let k: c_double = k.into();
    let sigma: c_double = sigma.into();
    let t: c_double = t.into();

    unsafe { black_ffi(f, k, sigma, t, q.into()) }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    #[test]
    fn test_black() {
        // arrange
        let f = 100.0;
        let k = 100.0;
        let sigma = 0.2;
        let t = 1.0;
        let q = OptionType::Call;

        // You will need to replace this with the expected result
        let expected_result = 7.9655674554057975;

        // act
        let result = black(f, k, sigma, t, q);

        // assert
        assert_eq!(result, expected_result);
    }

    #[test]
    fn test_black_zero_values() {
        let f = 0.0;
        let k = 0.0;
        let sigma = 0.0;
        let t = 0.0;
        let q = OptionType::Call;

        let result = black(f, k, sigma, t, q);

        // You will need to replace this with the expected result
        let expected_result = 0.0;

        assert_eq!(result, expected_result);
    }

    #[test]
    fn test_black_negative_values() {
        let f = -100.0;
        let k = -100.0;
        let sigma = -0.2;
        let t = -1.0;
        let q = OptionType::Put;

        let result = black(f, k, sigma, t, q);

        // You will need to replace this with the expected result
        let expected_result = 0.0;

        assert_eq!(result, expected_result);
    }

    #[test]
    fn test_black_high_sigma() {
        let f = 100.0;
        let k = 100.0;
        let sigma = 100.0;
        let t = 1.0;
        let q = OptionType::Call;

        let result = black(f, k, sigma, t, q);

        // You will need to replace this with the expected result
        let expected_result = 0.0;

        assert_eq!(result, expected_result);
    }

    #[test]
    fn test_black_high_t() {
        // arrange
        let f = 100.0;
        let k = 100.0;
        let sigma = 0.2;
        let t = 100.0;
        let q = OptionType::Call;

        // act
        let result = black(f, k, sigma, t, q);

        let expected_result = 68.26894921370861;

        // assert
        assert_eq!(result, expected_result);
    }

    fn prop_approx_eq(a: f64, b: f64, epsilon: f64) -> bool {
        (a - b).abs() < epsilon
    }

    #[test]
    fn test_black_function() {
        let f = 100.0;
        let k = 100.0;
        let sigma = 0.2;
        let t = 1.0;
        let q = OptionType::Call;
        let price = black(f, k, sigma, t, q);
        assert!(prop_approx_eq(price, 7.965567455405804, 1e-9));
    }

    #[test]
    fn test_implied_volatility() {
        let f = 100.0;
        let k = 100.0;
        let price = 7.965567455405804;
        let t = 1.0;
        let q = OptionType::Put;
        let sigma = implied_volatility_from_a_transformed_rational_guess(price, f, k, t, q);
        assert!(prop_approx_eq(sigma, 0.2, 1e-9));
    }

    #[test]
    fn test_black_function_typical_values() {
        let f = 100.0;
        let k = 100.0;
        let sigma = 0.2;
        let t = 1.0;
        let q = OptionType::Call;
        let price = black(f, k, sigma, t, q);
        assert!(prop_approx_eq(price, 7.965567455405804, 1e-9));
    }

    #[test]
    fn test_black_function_edge_case_high_volatility() {
        let f = 100.0;
        let k = 100.0;
        let sigma = 5.0;
        let t = 1.0;
        let q = OptionType::Put;
        let price = black(f, k, sigma, t, q);
        assert!(price > 0.0);
    }

    #[test]
    fn test_black_function_edge_case_zero_time() {
        let f = 100.0;
        let k = 100.0;
        let sigma = 0.2;
        let t = 0.0;
        let q = OptionType::Put;
        let price = black(f, k, sigma, t, q);
        assert_eq!(price, 0.0);
    }

    #[test]
    fn test_implied_volatility_typical_values() {
        let f = 100.0;
        let k = 100.0;
        let price = 7.965567455405804;
        let t = 1.0;
        let q = OptionType::Call;
        let result = implied_volatility_from_a_transformed_rational_guess(price, f, k, t, q);
        assert!(prop_approx_eq(result, 0.2, 1e-9));
    }

    #[test]
    fn test_implied_volatility_edge_case_high_price() {
        let f = 100.0;
        let k = 100.0;
        let price = 50.0;
        let t = 1.0;
        let q = OptionType::Put;
        let result = implied_volatility_from_a_transformed_rational_guess(price, f, k, t, q);
        assert!(result > 0.0);
    }

    #[test]
    fn test_implied_volatility_edge_case_zero_price() {
        let f = 100.0;
        let k = 100.0;
        let price = 0.01;
        let t = 1.0;
        let q = OptionType::Call;
        let result = implied_volatility_from_a_transformed_rational_guess(price, f, k, t, q);
        assert!(result > 0.0);
    }

    fn round_to(value: f64, places: u32) -> f64 {
        let factor = 10_f64.powi(places as i32);
        (value * factor).round() / factor
    }

    proptest! {
        #[test]
        fn test_black_function_random(f in 50.0_f64..150.0_f64, k in 50.0_f64..150.0_f64, sigma in 0.01_f64..2.0_f64, t in 0.1_f64..2.0_f64, q in any::<u8>().prop_map(|x| if x % 2 == 0 { OptionType::Call } else { OptionType::Put })) {
            let price = black(f, k, sigma, t, q);
            assert!(price >= 0.0);
        }

        // TODO: to discuss how to fix this test. Correspond to issue #8
        // #[test]
        fn test_implied_volatility_random(f in 50.0_f64..150.0_f64, k in 50.0_f64..150.0_f64, price in 0.01_f64..100.0_f64, t in 0.1_f64..2.0_f64, q in any::<u8>().prop_map(|x| if x % 2 == 0 { OptionType::Call } else { OptionType::Put })) {
            let iv = implied_volatility_from_a_transformed_rational_guess(price, f, k, t, q);

            let a = f64::total_cmp(&iv, &0.0);
            println!("Comparison result: {:?}", a);
            println!("Random Test - Implied Volatility: {:?}, f {:?}, k {:?} price {:?} t {:?} q {:?}", iv, f, k, price, t, q);
            assert!(!iv.is_sign_negative(), "Implied Volatility is negative");
            assert!(iv >= 0.0, "Implied Volatility is {}", iv);
        }
    }

    // TODO: to discuss how to fix this test. Correspond to issue #8
    // #[test]
    fn test_implied_volatility_specific() {
        let f = 137.9331040666909;
        let k = 50.0;
        let price = 0.01;
        let t = 0.1;
        let q = OptionType::Call;

        let result = implied_volatility_from_a_transformed_rational_guess(price, f, k, t, q);

        assert!(result.is_sign_negative() == false);
    }
}
