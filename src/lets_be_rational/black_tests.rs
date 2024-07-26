#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;

    use crate::lets_be_rational::black::{
        asymptotic_expansion_of_normalised_black_call,
        asymptotic_expansion_of_normalised_black_call_old,
        small_t_expansion_of_normalised_black_call, small_t_expansion_of_normalised_black_call_new,
    };

    #[test]
    fn compare_original_and_new_one() {
        let test_cases = [
            // correct data from the original implementation
            (-12.0, 0.1),
            (-20.0, 0.05),
            (-15.0, 0.2),
            (-30.0, 0.01),
            (-12.0, 0.1),
        ];

        for &(h, t) in &test_cases {
            let original = asymptotic_expansion_of_normalised_black_call_old(h, t).unwrap();
            let new_one = asymptotic_expansion_of_normalised_black_call(h, t)
                .expect("Wrong data or implementation - should happened here");

            println!(
                "OK: h: {}, t: {}, original: {}, optimized: {}, diff: {}",
                h,
                t,
                original,
                new_one,
                (original - new_one).abs()
            );
            assert_approx_eq!(original, new_one, 1e-10);
        }
    }

    #[test]
    fn test_small_t_expansion_comparison() {
        let test_cases = [
            (0.1, 0.1),
            (0.05, 0.05),
            (0.15, 0.15),
            (0.2, 0.2),
            (0.0, 0.0),
        ];

        for &(h, t) in test_cases.iter() {
            let result_new = Some(small_t_expansion_of_normalised_black_call(h, t).unwrap());
            let result_old = small_t_expansion_of_normalised_black_call_new(h, t).unwrap();

            match result_new {
                Some(value_new) => {
                    let value_old = result_old;
                    let tolerance = 1e-6;
                    println!(
                        "OK: h: {}, t: {}, original: {}, optimized: {}, diff: {}",
                        h,
                        t,
                        value_old,
                        value_new,
                        (value_old - value_new).abs()
                    );
                    assert_approx_eq!(value_new, value_old, tolerance);
                }
                None => {
                    assert!(t >= 0.21, "New implementation returned None for t < 0.21");
                }
            }
        }
    }

    #[test]
    fn test_t_above_threshold() {
        let h = 0.1;
        let t = 0.3;
        let result = small_t_expansion_of_normalised_black_call(h, t);
        assert!(result.is_ok());
    }
}
