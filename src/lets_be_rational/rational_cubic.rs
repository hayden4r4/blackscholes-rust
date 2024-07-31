use num_traits::Zero;

const MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE: f64 = -(1.0 - 1.4901161193847656e-8); // -(1.0 - f64::EPSILON.sqrt());
const MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE: f64 = 2.0 / (f64::EPSILON * f64::EPSILON);

pub fn rational_cubic_interpolation(
    x: f64,
    x_l: f64,
    x_r: f64,
    y_l: f64,
    y_r: f64,
    d_l: f64,
    d_r: f64,
    r: f64,
) -> Result<f64, String> {
    // TODO: verify is it necessary to check for infinite values
    // or ensure that they are not present in the input
    if x.is_infinite() || x_l.is_infinite() || x_r.is_infinite() {
        return Err("Infinite value found in input".to_string());
    }

    let h = x_r - x_l;
    // NOTE: This is an optimization as reuse calculated value of `h`.
    // Based on `additional_benches_to_verify::abs_vs_equality_benchmarks` as an evidence
    if h.abs() <= f64::EPSILON {
        return Ok(0.5 * (y_l + y_r));
    }

    let t = (x - x_l) / h;
    let omt = 1.0 - t;

    // r should be greater than -1.
    // We do not use assert(r > -1) here in order to allow values such as NaN to be propagated as they should.
    if r < MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE {
        let t_power = t * t;
        let omt2 = omt * omt;

        let numerator = y_r * t_power * t
            + (r * y_r - h * d_r) * t_power * omt
            + (r * y_l + h * d_l) * t * omt2
            + y_l * omt2 * omt;

        let denominator = 1.0 + (r - 3.0) * t * omt;

        // Formula (2.4) divided by formula (2.5)
        Ok(numerator / denominator)
    } else {
        // Linear interpolation without over-or underflow.
        // NOTE: `mul_add` is used to avoid the overhead of a separate multiplication and addition.
        Ok(y_r.mul_add(t, y_l * omt))
    }
}

fn minimum_rational_cubic_control_parameter(
    d_l: f64,
    d_r: f64,
    s: f64,
    prefer_shape_preservation_over_smoothness: bool,
) -> f64 {
    let monotonic = d_l * s >= 0.0 && d_r * s >= 0.0;
    let convex = d_l <= s && s <= d_r;
    let concave = d_l >= s && s >= d_r;
    if !monotonic && !convex && !concave {
        return MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE;
    }
    let d_r_m_d_l = d_r - d_l;
    let d_r_m_s = d_r - s;
    let s_m_d_l = s - d_l;

    // If monotonicity on this interval is possible,
    // set r1 to satisfy the monotonicity condition (3.8).
    let r1 = match (monotonic, s.is_zero()) {
        (true, false) => (d_r + d_l) / s,
        (true, true) if prefer_shape_preservation_over_smoothness => {
            MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
        }
        _ => f64::MIN,
    };

    let r2 = match (
        convex || concave,
        monotonic,
        prefer_shape_preservation_over_smoothness,
    ) {
        (true, _, _) if !s_m_d_l.is_zero() && !d_r_m_s.is_zero() => {
            f64::max((d_r_m_d_l / d_r_m_s).abs(), (d_r_m_d_l / s_m_d_l).abs())
        }
        (_, true, true) => MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE,
        _ => f64::MIN,
    };

    f64::max(
        MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE,
        f64::max(r1, r2),
    )
}

pub enum Side {
    Left,
    Right,
}

fn rational_cubic_control_parameter(
    x_l: f64,
    x_r: f64,
    y_l: f64,
    y_r: f64,
    d_l: f64,
    d_r: f64,
    second_derivative: f64,
    side: Side,
) -> f64 {
    let h = x_r - x_l;
    let numerator = 0.5 * second_derivative.mul_add(h, d_r - d_l);
    if numerator.is_zero() {
        return 0.0;
    }

    let denominator = match side {
        Side::Left => (y_r - y_l) / h - d_l,
        Side::Right => d_r - (y_r - y_l) / h,
    };

    match (
        denominator.is_zero() || denominator.is_infinite(),
        numerator > 0.0,
    ) {
        (true, true) => MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE,
        (true, false) => MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE,
        _ => numerator / denominator,
    }
}

pub fn convex_rational_cubic_control_parameter(
    x_l: f64,
    x_r: f64,
    y_l: f64,
    y_r: f64,
    d_l: f64,
    d_r: f64,
    second_derivative: f64,
    prefer_shape_preservation_over_smoothness: bool,
    side: Side,
) -> f64 {
    let r = rational_cubic_control_parameter(x_l, x_r, y_l, y_r, d_l, d_r, second_derivative, side);
    let r_min = minimum_rational_cubic_control_parameter(
        d_l,
        d_r,
        (y_r - y_l) / (x_r - x_l),
        prefer_shape_preservation_over_smoothness,
    );

    r.max(r_min)
}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;

    use super::*;

    const EPSILON: f64 = 1e-10;

    #[test]
    fn test_equal_x_values() {
        let result = rational_cubic_interpolation(1.0, 1.0, 1.0, 2.0, 3.0, 0.5, 0.5, 0.0)
            .expect("Test purpose");
        assert_approx_eq!(result, 2.5, 1e-6);
    }

    #[test]
    fn test_linear_interpolation() {
        let result = rational_cubic_interpolation(1.5, 1.0, 2.0, 1.0, 2.0, 1.0, 1.0, f64::MAX)
            .expect("Test purpose");
        assert_approx_eq!(result, 1.5, 1e-6);
    }

    #[test]
    fn test_cubic_interpolation() {
        let result = rational_cubic_interpolation(1.5, 1.0, 2.0, 1.0, 2.0, 0.0, 0.0, 0.0)
            .expect("Test purpose");
        assert_approx_eq!(result, 1.5, 1e-6);
    }

    #[test]
    fn test_extreme_r_value() {
        let result = rational_cubic_interpolation(
            1.5,
            1.0,
            2.0,
            1.0,
            2.0,
            1.0,
            1.0,
            MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE - 1e-10,
        )
        .expect("Test purpose");
        assert!(result.is_finite());
    }

    #[test]
    fn test_boundary_conditions() {
        let result_left = rational_cubic_interpolation(1.0, 1.0, 2.0, 1.0, 2.0, 0.5, 0.5, 0.0)
            .expect("Test purpose");
        assert_approx_eq!(result_left, 1.0, 1e-6);

        let result_right = rational_cubic_interpolation(2.0, 1.0, 2.0, 1.0, 2.0, 0.5, 0.5, 0.0)
            .expect("Test purpose");
        assert_approx_eq!(result_right, 2.0, 1e-6);
    }

    #[test]
    fn test_nan_propagation() {
        let result = rational_cubic_interpolation(f64::NAN, 1.0, 2.0, 1.0, 2.0, 0.5, 0.5, 0.0)
            .expect("Test purpose");
        assert!(result.is_nan());
    }

    #[test]
    fn test_infinity_handling() {
        let result = rational_cubic_interpolation(f64::INFINITY, 1.0, 2.0, 1.0, 2.0, 0.5, 0.5, 0.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_monotonic_increasing() {
        let result = minimum_rational_cubic_control_parameter(1.0, 2.0, 1.5, false);
        assert_approx_eq!(result, 2.0, EPSILON);
    }

    #[test]
    fn test_monotonic_decreasing() {
        let result = minimum_rational_cubic_control_parameter(-2.0, -1.0, -1.5, false);
        assert_approx_eq!(result, 2.0, EPSILON);
    }

    #[test]
    fn test_convex() {
        let result = minimum_rational_cubic_control_parameter(1.0, 3.0, 2.0, false);
        assert_approx_eq!(result, 2.0, EPSILON);
    }

    #[test]
    fn test_concave() {
        let result = minimum_rational_cubic_control_parameter(3.0, 1.0, 2.0, false);
        assert_approx_eq!(result, 2.0, EPSILON);
    }

    #[test]
    fn test_non_monotonic_non_convex_non_concave() {
        let result = minimum_rational_cubic_control_parameter(1.0, -1.0, 2.0, false);
        assert_approx_eq!(
            result,
            MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE,
            EPSILON
        );
    }

    #[test]
    fn test_prefer_shape_preservation_monotonic() {
        let result = minimum_rational_cubic_control_parameter(1.0, 2.0, 0.0, true);
        assert_approx_eq!(
            result,
            MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE,
            EPSILON
        );
    }

    #[test]
    fn test_prefer_shape_preservation_convex() {
        let result = minimum_rational_cubic_control_parameter(1.0, 3.0, 2.0, true);
        assert_approx_eq!(result, 2.0, EPSILON);
    }

    #[test]
    fn test_zero_s_max() {
        let result = minimum_rational_cubic_control_parameter(1.0, 2.0, 0.0, true);
        assert_approx_eq!(
            result,
            MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE,
            EPSILON
        );
    }

    #[test]
    fn test_zero_s_min() {
        let result = minimum_rational_cubic_control_parameter(1.0, 2.0, 0.0, false);
        assert_approx_eq!(
            result,
            MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE,
            EPSILON
        );
    }

    #[test]
    fn test_equal_d_l_and_d_r() {
        let result = minimum_rational_cubic_control_parameter(2.0, 2.0, 2.0, false);
        assert_approx_eq!(result, 2.0, EPSILON);
    }
}
