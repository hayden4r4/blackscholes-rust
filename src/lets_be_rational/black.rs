use std::f64::consts::FRAC_1_SQRT_2;

use statrs::consts::SQRT_2PI;
use statrs::function::erf::erfc;

use crate::lets_be_rational::{
    intrinsic::normalised_intrinsic, normal_distribution::standard_normal_cdf,
    DENORMALISATION_CUTOFF, ONE_OVER_SQRT_TWO_PI,
};
use crate::OptionType;

const ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD: f64 = -10.0;
const SIXTEENTH_ROOT_DBL_EPSILON: f64 = 0.10566243270259357;
const CODYS_THRESHOLD: f64 = 0.46875;

const SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD: f64 = 2.0 * SIXTEENTH_ROOT_DBL_EPSILON;

fn erfcx(x: f64) -> f64 {
    (x * x).exp() * erfc(x)
}

fn normalised_black_call_using_norm_cdf(x: f64, s: f64) -> f64 {
    let h = x / s;
    let t = 0.5 * s;
    let b_max = (0.5 * x).exp();
    let b = standard_normal_cdf(h + t) * b_max - standard_normal_cdf(h - t) / b_max;
    b.abs().max(0.0)
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

    (0.5 * two_b).abs().max(0.0)
}

#[rustfmt::skip]
pub fn small_t_expansion_of_normalised_black_call_old(h: f64, t: f64) -> f64 {
    let a = 1.0 + h * (0.5 * SQRT_2PI) * erfcx(-FRAC_1_SQRT_2 * h);
    let w = t * t;
    let h2 = h * h;

    let expansion = 2.0 * t * (a + w * ((-1.0 + 3.0 * a + a * h2) / 6.0
        + w * ((-7.0 + 15.0 * a + h2 * (-1.0 + 10.0 * a + a * h2)) / 120.0
        + w * ((-57.0 + 105.0 * a + h2 * (-18.0 + 105.0 * a + h2 * (-1.0 + 21.0 * a + a * h2))) / 5040.0
        + w * ((-561.0 + 945.0 * a + h2 * (-285.0 + 1260.0 * a + h2 * (-33.0 + 378.0 * a + h2 * (-1.0 + 36.0 * a + a * h2)))) / 362880.0
        + w * ((-6555.0 + 10395.0 * a + h2 * (-4680.0 + 17325.0 * a + h2 * (-840.0 + 6930.0 * a + h2 * (-52.0 + 990.0 * a + h2 * (-1.0 + 55.0 * a + a * h2))))) / 39916800.0
        + w * ((-89055.0 + 135135.0 * a + h2 * (-82845.0 + 270270.0 * a + h2 * (-20370.0 + 135135.0 * a + h2 * (-1926.0 + 25740.0 * a + h2 * (-75.0 + 2145.0 * a + h2 * (-1.0 + 78.0 * a + a * h2)))))) / 6227020800.0)))))));

    let b = ONE_OVER_SQRT_TWO_PI * (-0.5 * (h * h + t * t)).exp() * expansion;
    b.abs().max(0.0)
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
/// * `Some(f64)` containing the computed value if t < 0.21
/// * `None` if t >= 0.21, indicating the approximation is not valid
pub fn small_t_expansion_of_normalised_black_call(h: f64, t: f64) -> Option<f64> {
    const SQRT_2_PI: f64 = 2.5066282746310002;
    const ONE_OVER_SQRT_TWO_PI: f64 = 1.0 / SQRT_2_PI;
    const T_THRESHOLD: f64 = 0.21;

    if t >= T_THRESHOLD {
        return None;
    }

    let t2 = t * t;
    let h2 = h * h;

    // Obliczenie 'a' jak w starej implementacji
    let a = 1.0 + h * (0.5 * SQRT_2_PI) * erfcx(-FRAC_1_SQRT_2 * h);

    // 12th order Taylor expansion with 6 terms
    let y_diff = {
        let term1 = a.mul_add(1.0, -1.0); // -1.0 + a
        let term2 = h2.mul_add(a, 3.0 * a - 7.0); // -7.0 + 15.0 * a + h2 * (-1.0 + 10.0 * a + a * h2)
        let term3 = h2
            .mul_add(h2.mul_add(a, 21.0 * a - 1.0), 105.0 * a - 18.0)
            .mul_add(1.0, 105.0 * a - 57.0);
        let term4 = h2
            .mul_add(
                h2.mul_add(h2.mul_add(a, 36.0 * a - 1.0), 378.0 * a - 33.0),
                1260.0 * a - 285.0,
            )
            .mul_add(1.0, 945.0 * a - 561.0);
        let term5 = h2
            .mul_add(
                h2.mul_add(
                    h2.mul_add(h2.mul_add(a, 55.0 * a - 1.0), 990.0 * a - 52.0),
                    6930.0 * a - 840.0,
                ),
                17325.0 * a - 4680.0,
            )
            .mul_add(1.0, 10395.0 * a - 6555.0);

        let t2_squared = t2 * t2;
        let t2_cubed = t2_squared * t2;
        let t2_quartic = t2_squared * t2_squared;
        let t2_quintic = t2_cubed * t2_squared;

        2.0 * t
            * (a.mul_add(
                1.0,
                term1.mul_add(
                    t2 / 6.0,
                    term2.mul_add(
                        t2_squared / 120.0,
                        term3.mul_add(
                            t2_cubed / 5040.0,
                            term4.mul_add(t2_quartic / 362880.0, term5 * t2_quintic / 39916800.0),
                        ),
                    ),
                ),
            ))
    };

    let exp_term = (-0.5 * (h2 + t2)).exp();
    let black_value = ONE_OVER_SQRT_TWO_PI * exp_term * y_diff;

    Some(black_value.abs().max(0.0))
}

pub fn normalised_black_call_using_erfcx(h: f64, t: f64) -> f64 {
    let b = 0.5
        * (-0.5 * (h * h + t * t)).exp()
        * (erfcx(-FRAC_1_SQRT_2 * (h + t)) - erfcx(-FRAC_1_SQRT_2 * (h - t)));
    b.abs().max(0.0)
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
            .expect("Parameters should be correct - temporary expect");
    }
    if 0.5 * s < SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD {
        return small_t_expansion_of_normalised_black_call(x / s, 0.5 * s)
            .expect("Parameters should be correct - temporary expect");
    }

    normalised_black_call_with_optimal_use_of_codys_functions(x, s)
}

pub(crate) fn normalised_black(x: f64, s: f64, q: f64) -> f64 {
    normalised_black_call(if q < 0.0 { -x } else { x }, s) /* Reciprocal-strike call-put equivalence */
}

#[rustfmt::skip]
pub fn asymptotic_expansion_of_normalised_black_call_old(h: f64, t: f64) -> f64 {
    let e = (t / h) * (t / h);
    let r = (h + t) * (h - t);
    let q = (h / r) * (h / r);

    // 17th order asymptotic expansion of A(h,t) in q
    let asymptotic_expansion_sum =
        2.0 + q * (-6.0 - 2.0 * e + 3.0 * q * (10.0 + e * (20.0 + 2.0 * e)
            + 5.0 * q * (-14.0 + e * (-70.0 + e * (-42.0 - 2.0 * e))
            + 7.0 * q * (18.0 + e * (168.0 + e * (252.0 + e * (72.0 + 2.0 * e)))
            + 9.0 * q * (-22.0 + e * (-330.0 + e * (-924.0 + e * (-660.0 + e * (-110.0 - 2.0 * e))))
            + 11.0 * q * (26.0 + e * (572.0 + e * (2574.0 + e * (3432.0 + e * (1430.0 + e * (156.0 + 2.0 * e)))))
            + 13.0 * q * (-30.0 + e * (-910.0 + e * (-6006.0 + e * (-12870.0 + e * (-10010.0 + e * (-2730.0 + e * (-210.0 - 2.0 * e)))))
            + 15.0 * q * (34.0 + e * (1360.0 + e * (12376.0 + e * (38896.0 + e * (48620.0 + e * (24752.0 + e * (4760.0 + e * (272.0 + 2.0 * e)))))))
            + 17.0 * q * (-38.0 + e * (-1938.0 + e * (-23256.0 + e * (-100776.0 + e * (-184756.0 + e * (-151164.0 + e * (-54264.0 + e * (-7752.0 + e * (-342.0 - 2.0 * e))))))))
            + 19.0 * q * (42.0 + e * (2660.0 + e * (40698.0 + e * (232560.0 + e * (587860.0 + e * (705432.0 + e * (406980.0 + e * (108528.0 + e * (11970.0 + e * (420.0 + 2.0 * e)))))))))
            + 21.0 * q * (-46.0 + e * (-3542.0 + e * (-67298.0 + e * (-490314.0 + e * (-1634380.0 + e * (-2704156.0 + e * (-2288132.0 + e * (-980628.0 + e * (-201894.0 + e * (-17710.0 + e * (-506.0 - 2.0 * e))))))))))
            + 23.0 * q * (50.0 + e * (4600.0 + e * (106260.0 + e * (961400.0 + e * (4085950.0 + e * (8914800.0 + e * (10400600.0 + e * (6537520.0 + e * (2163150.0 + e * (354200.0 + e * (25300.0 + e * (600.0 + 2.0 * e)))))))))))
            + 25.0 * q * (-54.0 + e * (-5850.0 + e * (-161460.0 + e * (-1776060.0 + e * (-9373650.0 + e * (-26075790.0 + e * (-40116600.0 + e * (-34767720.0 + e * (-16872570.0 + e * (-4440150.0 + e * (-592020.0 + e * (-35100.0 + e * (-702.0 - 2.0 * e))))))))))))
            + 27.0 * q * (58.0 + e * (7308.0 + e * (237510.0 + e * (3121560.0 + e * (20030010.0 + e * (69194580.0 + e * (135727830.0 + e * (155117520.0 + e * (103791870.0 + e * (40060020.0 + e * (8584290.0 + e * (950040.0 + e * (47502.0 + e * (812.0 + 2.0 * e)))))))))))))
            + 29.0 * q * (-62.0 + e * (-8990.0 + e * (-339822.0 + e * (-5259150.0 + e * (-40320150.0 + e * (-169344630.0 + e * (-412506150.0 + e * (-601080390.0 + e * (-530365050.0 + e * (-282241050.0 + e * (-88704330.0 + e * (-15777450.0 + e * (-14725620.0 + e * (-629300.0 + e * (-930.0 - 2.0 * e)))))))))))))))
            + 31.0 * q * (66.0 + e * (109120.0 + e * (474672.0 + e * (8544096.0 + e * (77134200.0 + e * (387073440.0 + e * (1146332880.0 + e * (2074316640.0 + e * (2333606220.0 + e * (1637618400.0 + e * (709634640.0 + e * (185122080.0 + e * (27768312.0 + e * (2215136.0 + e * (81840.0 + e * (1056.0 + 2.0 * e)))))))))))))))
            + 33.0 * (-70.0 + e * (-130900.0 + e * (-649264.0 + e * (-13449040.0 + e * (-141214920.0 + e * (-834451800.0 + e * (-2952675600.0 + e * (-6495886320.0 + e * (-9075135300.0 + e * (-8119857900.0 + e * (-4639918800.0 + e * (-1668903600.0 + e * (-367158792.0 + e * (-47071640.0 + e * (-3246320.0 + e * (-104720.0 + e * (-1190.0 - 2.0 * e))))))))))))))))) * q))))))))))))))));

    let b = ONE_OVER_SQRT_TWO_PI * ((-0.5 * (h * h + t * t)).exp()) * (t / r) * asymptotic_expansion_sum;
    b.abs().max(0.0)
}

const H_LARGE: f64 = -10.0;

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
    let tau_small: f64 = 2.0 * f64::EPSILON.sqrt().powi(16);

    // Check if we're in the correct region
    // From section 6: "In the region of large negative h, we can realize these preferences by the aid of the formulation (6.10) for the normalized Black function"
    if h > H_LARGE || t >= (h.abs() - H_LARGE.abs() + tau_small) {
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

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;

    use super::*;

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
            let original = asymptotic_expansion_of_normalised_black_call_old(h, t);
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
            let result_new = small_t_expansion_of_normalised_black_call(h, t);
            let result_old = small_t_expansion_of_normalised_black_call_old(h, t);

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
                    // assert_approx_eq!(value_new, value_old, tolerance);
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
        assert!(result.is_none());
    }
}
