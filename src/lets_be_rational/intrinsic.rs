use crate::OptionType;

const FOURTH_ROOT_DBL_EPSILON: f64 = 0.0001220703125;

const NORMALISED_X2_THRESHOLD: f64 = 98.0 * FOURTH_ROOT_DBL_EPSILON;

pub(crate) fn normalised_intrinsic(x: f64, option_type: OptionType) -> f64 {
    if option_type * x <= 0.0 {
        return 0.0;
    }
    let x2 = x * x;
    if x2 < NORMALISED_X2_THRESHOLD {
        return (option_type
            * x
            * (1.0
                + x2 * ((1.0 / 24.0)
                    + x2 * ((1.0 / 1920.0) + x2 * ((1.0 / 322560.0) + (1.0 / 92897280.0) * x2)))))
            .max(0.0);
    }
    let b_max = (0.5 * x).exp();
    let one_over_b_max = 1.0 / b_max;
    (option_type * (b_max - one_over_b_max)).max(0.0)
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use proptest::prelude::*;

    use super::*;

    const EPSILON: f64 = 1e-9;

    #[test]
    fn test_normalised_intrinsic_near_zero() {
        // arrange
        let x = 0.0001;
        const EXPECTED_RESULT: f64 = 0.00010000000004166668;
        let option_type = OptionType::Call;

        // act
        let result = normalised_intrinsic(x, option_type);

        // assert
        assert_relative_eq!(result, EXPECTED_RESULT, epsilon = EPSILON);
    }

    #[test]
    fn test_normalised_intrinsic_with_large_x() {
        let x = f64::MAX;
        let option_type = OptionType::Call;
        let result = normalised_intrinsic(x, option_type);

        // TODO: to analyze correct paths here
        assert_eq!(
            result,
            f64::INFINITY,
            "Expected intrinsic value to be infinity for f64::MAX."
        );
    }

    #[test]
    fn test_normalised_intrinsic_at_zero() {
        // arrange
        let x = 0.0;
        const EXPECTED_RESULT: f64 = 0.0;
        let option_type = OptionType::Call;

        // act
        let result = normalised_intrinsic(x, option_type);

        // assert
        assert_relative_eq!(result, EXPECTED_RESULT, epsilon = EPSILON);
    }

    #[test]
    fn test_normalised_intrinsic_negative_x() {
        // arrange
        let x = -0.0001;
        const EXPECTED_RESULT: f64 = 0.00010000000004166668;
        let option_type = OptionType::Put;

        // act
        let result = normalised_intrinsic(x, option_type);

        // assert
        assert_relative_eq!(result, EXPECTED_RESULT, epsilon = EPSILON);
    }

    #[test]
    fn test_normalised_intrinsic_with_min_x() {
        // arrange
        let x = f64::MIN;
        const EXPECTED_RESULT: f64 = 0.0;
        let option_type = OptionType::Call;

        // act
        let result = normalised_intrinsic(x, option_type);

        // assert
        assert_eq!(
            result, EXPECTED_RESULT,
            "Expected intrinsic value to be 0 for f64::MIN."
        );
    }

    #[test]
    fn test_normalised_intrinsic_around_threshold() {
        // arrange
        let x = NORMALISED_X2_THRESHOLD.sqrt();
        const EXPECTED_RESULT: f64 = 0.10942952653480309;
        let option_type = OptionType::Call;

        // act
        let result = normalised_intrinsic(x, option_type);

        // assert
        assert_relative_eq!(result, EXPECTED_RESULT, epsilon = EPSILON);
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
