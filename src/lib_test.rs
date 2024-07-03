#[cfg(test)]
mod tests {
    use crate::{*};
    #[test]
    fn test_calc_npdf() {
        // arrange
        let test_cases = vec![
            (0.5, 0.352065325),
            (0.0, 0.3989422804),
            (-1.0, 0.241970733),
            (10.0, 1e-5),
            (-10.0, 1e-5),
        ];

        for (x, expected) in test_cases {
            // act
            let result = calc_npdf(x);

            // assert
            // floating point result for x larger than 10 is expected to be very small not comparable
            if x.abs() < 10.0 {
                assert_eq!(result, expected);
            } else {
                assert!(result < expected, "Expected a very small result for large x");
            }
        }
    }

    #[test]
    fn test_calc_nprimed1() {
        // TODO: after deciding of correct behavior
    }

    #[test]
    fn test_calc_nprimed2() {
        // TODO: after deciding of correct behavior
    }
}