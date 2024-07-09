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
) -> f64 {
    let h = x_r - x_l;
    if h.abs() <= 0.0 {
        return 0.5 * (y_l + y_r);
    }
    // r should be greater than -1.
    // We do not use assert(r > -1) here in order to allow values such as NaN to be propagated as they should.
    let t = (x - x_l) / h;
    if r < MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE {
        let omt = 1.0 - t;
        let t2 = t * t;
        let omt2 = omt * omt;
        // Formula (2.4) divided by formula (2.5)
        return (y_r * t2 * t
            + (r * y_r - h * d_r) * t2 * omt
            + (r * y_l + h * d_l) * t * omt2
            + y_l * omt2 * omt)
            / (1.0 + (r - 3.0) * t * omt);
    }
    // Linear interpolation without over-or underflow.
    y_r * t + y_l * (1.0 - t)
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
    let mut r1 = f64::NEG_INFINITY;
    let mut r2 = r1;

    // If monotonicity on this interval is possible, set r1 to satisfy the monotonicity condition (3.8).
    if monotonic {
        if !s.is_zero() {
            r1 = (d_r + d_l) / s;
        } else if prefer_shape_preservation_over_smoothness {
            r1 = MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE;
        }
    }
    if convex || concave {
        if !s_m_d_l.is_zero() && !d_r_m_s.is_zero() {
            r2 = f64::max((d_r_m_d_l / d_r_m_s).abs(), (d_r_m_d_l / s_m_d_l).abs());
        } else if prefer_shape_preservation_over_smoothness {
            r2 = MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE;
        }
    } else if monotonic && prefer_shape_preservation_over_smoothness {
        r2 = MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE;
    }
    f64::max(
        MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE,
        f64::max(r1, r2),
    )
}

fn rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(
    x_l: f64,
    x_r: f64,
    y_l: f64,
    y_r: f64,
    d_l: f64,
    d_r: f64,
    second_derivative_l: f64,
) -> f64 {
    let h = x_r - x_l;
    let numerator = 0.5 * h * second_derivative_l + (d_r - d_l);
    if numerator.is_zero() {
        return 0.0;
    }
    let denominator = (y_r - y_l) / h - d_l;
    if denominator.is_zero() {
        return if numerator > 0.0 {
            MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
        } else {
            MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
        };
    }
    numerator / denominator
}

fn rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(
    x_l: f64,
    x_r: f64,
    y_l: f64,
    y_r: f64,
    d_l: f64,
    d_r: f64,
    second_derivative_r: f64,
) -> f64 {
    let h = x_r - x_l;
    let numerator = 0.5 * h * second_derivative_r + (d_r - d_l);
    if numerator.is_zero() {
        return 0.0;
    }
    let denominator = d_r - (y_r - y_l) / h;
    if denominator.is_zero() {
        return if numerator > 0.0 {
            MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
        } else {
            MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE
        };
    }
    numerator / denominator
}

pub fn convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(
    x_l: f64,
    x_r: f64,
    y_l: f64,
    y_r: f64,
    d_l: f64,
    d_r: f64,
    second_derivative_l: f64,
    prefer_shape_preservation_over_smoothness: bool,
) -> f64 {
    let r = rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(
        x_l,
        x_r,
        y_l,
        y_r,
        d_l,
        d_r,
        second_derivative_l,
    );
    let r_min = minimum_rational_cubic_control_parameter(
        d_l,
        d_r,
        (y_r - y_l) / (x_r - x_l),
        prefer_shape_preservation_over_smoothness,
    );
    r.max(r_min)
}

pub fn convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(
    x_l: f64,
    x_r: f64,
    y_l: f64,
    y_r: f64,
    d_l: f64,
    d_r: f64,
    second_derivative_r: f64,
    prefer_shape_preservation_over_smoothness: bool,
) -> f64 {
    let r = rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(
        x_l,
        x_r,
        y_l,
        y_r,
        d_l,
        d_r,
        second_derivative_r,
    );
    let r_min = minimum_rational_cubic_control_parameter(
        d_l,
        d_r,
        (y_r - y_l) / (x_r - x_l),
        prefer_shape_preservation_over_smoothness,
    );
    r.max(r_min)
}
