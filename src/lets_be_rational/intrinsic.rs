use crate::OptionType;

const FOURTH_ROOT_DBL_EPSILON: f64 = 0.0001220703125;

const NORMALISED_X2_THRESHOLD: f64 = 98.0 * FOURTH_ROOT_DBL_EPSILON;

pub(crate) fn normalised_intrinsic(x: f64, option_type: OptionType) -> f64 {
    let q = option_type as i32 as f64;
    if q * x <= 0.0 {
        return 0.0;
    }
    let x2 = x * x;
    if x2 < NORMALISED_X2_THRESHOLD {
        return ((if q < 0.0 { -1.0 } else { 1.0 }) * x *
            (1.0 + x2 * ((1.0 / 24.0) + x2 * ((1.0 / 1920.0) + x2 * ((1.0 / 322560.0) + (1.0 / 92897280.0) * x2))))).max(0.0).abs();
    }
    let b_max = (0.5 * x).exp();
    let one_over_b_max = 1.0 / b_max;
    ((if q < 0.0 { -1.0 } else { 1.0 }) * (b_max - one_over_b_max)).max(0.0).abs()
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    #[test]
    fn test_normalised_intrinsic_near_zero() {
        let x = 0.0001;
        let option_type = OptionType::Call;
        let result = normalised_intrinsic(x, option_type);
        assert!(result > 0.0 && result < 0.0002, "Expected a very small positive intrinsic value for x near zero.");
    }

    #[test]
    fn test_normalised_intrinsic_with_large_x() {
        let x = f64::MAX;
        let option_type = OptionType::Call;
        let result = normalised_intrinsic(x, option_type);
        assert!(result > 0.0, "Expected a positive intrinsic value for a large x.");
    }

    #[test]
    fn test_normalised_intrinsic_at_zero() {
        let x = 0.0;
        let option_type = OptionType::Call;
        let result = normalised_intrinsic(x, option_type);
        assert_eq!(result, 0.0, "Expected intrinsic value to be 0 for x = 0.");
    }

    #[test]
    fn test_normalised_intrinsic_negative_x() {
        let x = -0.0001;
        let option_type = OptionType::Put;
        let result = normalised_intrinsic(x, option_type);
        assert!(result > 0.0, "Expected a positive intrinsic value for a small negative x.");
    }

    #[test]
    fn test_normalised_intrinsic_with_min_x() {
        let x = f64::MIN;
        let option_type = OptionType::Call;
        let result = normalised_intrinsic(x, option_type);
        assert_eq!(result, 0.0, "Expected intrinsic value to be 0 for f64::MIN.");
    }

    #[test]
    fn test_normalised_intrinsic_around_threshold() {
        // arrange
        let x = NORMALISED_X2_THRESHOLD.sqrt();
        const EXPECTED_RESULT: f64 = 0.1;
        let option_type = OptionType::Call;
        let result = normalised_intrinsic(x, option_type);

        assert!(prop_approx_eq(result, 0.2, 1e-9));
        assert!(result >= 0.0, "Expected non-negative intrinsic value around threshold.");
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(100000))]
        #[test]
        fn test_normalised_intrinsic_always_positive_for_valid_inputs(
            x in any::<f64>(),
            option_type in any::<u8>().prop_map(|x| if x % 2 == 0 { OptionType::Call } else { OptionType::Put })) {
            let result = normalised_intrinsic(x, option_type);
            prop_assert!(result >= 0.0, "Expected a non-negative intrinsic value. Got: {}, for x: {}, option_type: {:?}", result, x, option_type);
        }
    }
}