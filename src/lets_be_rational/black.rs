use std::f64::consts::FRAC_1_SQRT_2;

use statrs::consts::SQRT_2PI;

use crate::{
    lets_be_rational::{
        cody::optimized::{erfc, erfcx},
        intrinsic::normalised_intrinsic,
        normal_distribution::standard_normal_cdf,
        DENORMALISATION_CUTOFF, ONE_OVER_SQRT_TWO_PI,
    },
    OptionType,
};

const H_LARGE: f64 = -10.0;

const ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD: f64 = -10.0;
const SIXTEENTH_ROOT_DBL_EPSILON: f64 = 0.10566243270259357;
const CODYS_THRESHOLD: f64 = 0.46875;

const SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD: f64 = 2.0 * SIXTEENTH_ROOT_DBL_EPSILON;

#[allow(dead_code)]
fn normalised_black_call_using_norm_cdf(x: f64, s: f64) -> f64 {
    let h = x / s;
    let t = 0.5 * s;
    let b_max = (0.5 * x).exp();
    let b = standard_normal_cdf(h + t) * b_max - standard_normal_cdf(h - t) / b_max;
    b.max(0.0)
}

fn normalised_black_call_with_optimal_use_of_codys_functions(x: f64, s: f64) -> f64 {
    let h = x / s;
    let t = 0.5 * s;
    let q1 = -FRAC_1_SQRT_2 * (h + t);
    let q2 = -FRAC_1_SQRT_2 * (h - t);
    let two_b;

    if q1 < CODYS_THRESHOLD {
        if q2 < CODYS_THRESHOLD {
            two_b = (0.5 * x).exp() * erfc(q1) - (-0.5 * x).exp() * erfc(q2);
        } else {
            two_b = (0.5 * x).exp() * erfc(q1) - (-0.5 * (h * h + t * t)).exp() * erfcx(q2);
        }
    } else if q2 < CODYS_THRESHOLD {
        two_b = (-0.5 * (h * h + t * t)).exp() * erfcx(q1) - (-0.5 * x).exp() * erfc(q2);
    } else {
        two_b = (-0.5 * (h * h + t * t)).exp() * (erfcx(q1) - erfcx(q2));
    }

    (0.5 * two_b).max(0.0)
}

/// Computes the normalized Black call option value using a small-t expansion.
///
/// # Arguments
///
/// * `h` - The normalized log-moneyness (x / σ)
/// * `t` - Half the total volatility (σ / 2)
///
/// # Returns
///
/// * `Some(f64)` containing the computed value if t < 0.21132486540518713
/// * `None` otherwise, indicating the approximation is not valid
pub fn small_t_expansion_of_normalised_black_call(h: f64, t: f64) -> Option<f64> {
    if t >= SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD {
        return None;
    }

    let t2 = t * t;
    let h2 = h * h;

    let a = 1.0 + h * (0.5 * SQRT_2PI) * erfcx(-FRAC_1_SQRT_2 * h);

    // 12th order Taylor expansion with 6 terms
    let y_diff = {
        let term1 = a.mul_add(h2 + 3.0, -1.0) / 6.0;
        let term2 = h2.mul_add(a.mul_add(h2 + 10.0, -1.0), a.mul_add(15.0, -7.0)) / 120.0;
        let term3 = h2.mul_add(
            h2.mul_add(a.mul_add(h2 + 21.0, -1.0), a.mul_add(105.0, -18.0)),
            a.mul_add(105.0, -57.0),
        ) / 5040.0;
        let term4 = h2.mul_add(
            h2.mul_add(
                h2.mul_add(a.mul_add(h2 + 36.0, -1.0), a.mul_add(378.0, -33.0)),
                a.mul_add(1260.0, -285.0),
            ),
            a.mul_add(945.0, -561.0),
        ) / 362880.0;

        let term5 = h2.mul_add(
            h2.mul_add(
                h2.mul_add(
                    h2.mul_add(a.mul_add(h2 + 55.0, -1.0), a.mul_add(990.0, -52.0)),
                    a.mul_add(6930.0, -840.0),
                ),
                a.mul_add(17325.0, -4680.0),
            ),
            a.mul_add(10395.0, -6555.0),
        ) / 39916800.0;
        let term6 = h2.mul_add(
            h2.mul_add(
                h2.mul_add(
                    h2.mul_add(
                        h2.mul_add(a.mul_add(h2 + 78.0, -1.0), a.mul_add(2145.0, -75.0)),
                        a.mul_add(25740.0, -1926.0),
                    ),
                    a.mul_add(135135.0, -20370.0),
                ),
                a.mul_add(270270.0, -82845.0),
            ),
            a.mul_add(135135.0, -89055.0),
        ) / 6227020800.0;

        let t2_squared = t2 * t2;
        let t2_cubed = t2_squared * t2;
        let t2_quartic = t2_squared * t2_squared;
        let t2_quintic = t2_cubed * t2_squared;
        let t2_sextic = t2_cubed * t2_cubed;

        2.0 * t
            * (a.mul_add(
                1.0,
                term1.mul_add(
                    t2,
                    term2.mul_add(
                        t2_squared,
                        term3.mul_add(
                            t2_cubed,
                            term4.mul_add(t2_quartic, term5.mul_add(t2_quintic, term6 * t2_sextic)),
                        ),
                    ),
                ),
            ))
    };

    let exp_term = (-0.5 * (h2 + t2)).exp();
    let black_value = ONE_OVER_SQRT_TWO_PI * exp_term * y_diff;

    Some(black_value.max(0.0))
}

#[allow(dead_code)]
pub fn normalised_black_call_using_erfcx(h: f64, t: f64) -> f64 {
    let b = 0.5
        * (-0.5 * (h * h + t * t)).exp()
        * (erfcx(-FRAC_1_SQRT_2 * (h + t)) - erfcx(-FRAC_1_SQRT_2 * (h - t)));
    b.max(0.0)
}

pub(crate) fn normalised_black_call(x: f64, s: f64) -> f64 {
    if x > 0.0 {
        return normalised_intrinsic(x, OptionType::Call) + normalised_black_call(-x, s);
        // In the money.
    }
    if s <= x.abs() * DENORMALISATION_CUTOFF {
        return normalised_intrinsic(x, OptionType::Call); // sigma=0 -> intrinsic value.
    }

    // Denote h := x/s and t := s/2.
    // We evaluate the condition |h|>|η|, i.e., h<η  &&  t < τ+|h|-|η|  avoiding any divisions by s , where η = asymptotic_expansion_accuracy_threshold  and τ = small_t_expansion_of_normalised_black_threshold .
    if x < s * ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD
        && 0.5 * s * s + x
            < s * (SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD
                + ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD)
    {
        return asymptotic_expansion_of_normalised_black_call(x / s, 0.5 * s)
            .expect("Asymptotic expansion failed");
    }
    if 0.5 * s < SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD {
        return small_t_expansion_of_normalised_black_call(x / s, 0.5 * s)
            .expect("Small t expansion failed");
    }

    normalised_black_call_with_optimal_use_of_codys_functions(x, s)
}

pub(crate) fn normalised_black(x: f64, s: f64, option_type: OptionType) -> f64 {
    normalised_black_call(option_type * x, s) /* Reciprocal-strike call-put equivalence */
}
/// Computes the asymptotic expansion of the normalized Black call price.
///
/// # Arguments
///
/// * `h` - The normalized log-moneyness (x/σ)
/// * `t` - Half the normalized volatility (σ/2)
///
/// # Returns
///
/// * `Ok(f64)` - The computed price
/// * `Err(String)` - An error message if the input is out of the valid range
pub fn asymptotic_expansion_of_normalised_black_call(h: f64, t: f64) -> Result<f64, &'static str> {
    // Check if we're in the correct region
    // From section 6: "In the region of large negative h, we can realize these preferences by the aid of the formulation (6.10) for the normalized Black function"
    if h > H_LARGE
        || t >= (h.abs() - H_LARGE.abs() + SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD)
    {
        return Err("This asymptotic expansion is only valid for large negative h and small t");
    }

    let h_plus_t = h + t;
    let h_minus_t = h - t;
    let h_plus_t_sq = h_plus_t * h_plus_t;
    let h_minus_t_sq = h_minus_t * h_minus_t;

    // Initialize the series terms
    // From equation (6.13): "Y(z) ≈ 1/z − 1/z^3 + 1·3/z^5 + ..."
    let mut y_h_plus_t = 1.0 / h_plus_t - 1.0 / (h_plus_t_sq * h_plus_t)
        + 3.0 / (h_plus_t_sq * h_plus_t_sq * h_plus_t);
    let mut y_h_minus_t = 1.0 / h_minus_t - 1.0 / (h_minus_t_sq * h_minus_t)
        + 3.0 / (h_minus_t_sq * h_minus_t_sq * h_minus_t);

    // Initial factorial and sign values for the series expansion
    let mut factorial = 15.0; // 3! = 6, next factorial would be for 5!
    let mut sign = -1.0;

    let mut term_plus = 1.0 / (factorial * h_plus_t_sq * h_plus_t_sq * h_plus_t);
    let mut term_minus = 1.0 / (factorial * h_minus_t_sq * h_minus_t_sq * h_minus_t);

    // Compute the series expansion up to the 17th term
    // From section 6: "We found that for n = 17, the approximation series has a maximum relative error of 1.64 · 10^−16 for all z ≤ −10"
    for i in 3..=17 {
        y_h_plus_t += sign * term_plus;
        y_h_minus_t += sign * term_minus;

        // Update sign and factorial for the next term
        sign = -sign;
        let i_f64 = i as f64;
        let factor = (2.0 * i_f64 + 2.0) * (2.0 * i_f64 + 3.0);
        factorial *= factor;
        term_plus *= h_plus_t_sq / factor;
        term_minus *= h_minus_t_sq / factor;
    }

    // Calculate the difference of the series
    // From equation (6.12): "[Y(h + t) - Y(h - t)]"
    let diff_y = y_h_plus_t - y_h_minus_t;

    // Exponential term
    // From equation (6.10): "e^(−1/2(h^2+t^2))"
    let exp_term = f64::mul_add(-0.5, h * h, -0.5 * t * t).exp();

    // Return the final computed price
    // From equation (6.10): "b = 1/√(2π) · e^(−1/2(h^2+t^2)) · [Y(h + t) - Y(h - t)]"
    Ok((exp_term * diff_y / SQRT_2PI).abs())
}
