use std::f64::consts::FRAC_1_SQRT_2;

use statrs::function::erf::erfc;

use crate::lets_be_rational::{DENORMALISATION_CUTOFF, ONE_OVER_SQRT_TWO_PI, SQRT_TWO_PI};
use crate::lets_be_rational::intrinsic::normalised_intrinsic;
use crate::lets_be_rational::normal_distribution::standard_normal_cdf;
use crate::lets_be_rational::so_rational::asymptotic_expansion_of_normalised_black_call;
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

fn small_t_expansion_of_normalised_black_call(h: f64, t: f64) -> f64 {
    let a = 1.0 + h * (0.5 * SQRT_TWO_PI) * erfcx(-FRAC_1_SQRT_2 * h);
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

pub fn normalised_black_call_using_erfcx(h: f64, t: f64) -> f64 {
    let b = 0.5 * (-0.5 * (h * h + t * t)).exp() * (erfcx(-FRAC_1_SQRT_2 * (h + t)) - erfcx(-FRAC_1_SQRT_2 * (h - t)));
    b.abs().max(0.0)
}

pub(crate) fn normalised_black_call(x: f64, s: f64) -> f64 {
    if x > 0.0 {
        return normalised_intrinsic(x, OptionType::Call) + normalised_black_call(-x, s); // In the money.
    }
    if s <= x.abs() * DENORMALISATION_CUTOFF {
        return normalised_intrinsic(x, OptionType::Call); // sigma=0 -> intrinsic value.
    }

    // Denote h := x/s and t := s/2.
    // We evaluate the condition |h|>|η|, i.e., h<η  &&  t < τ+|h|-|η|  avoiding any divisions by s , where η = asymptotic_expansion_accuracy_threshold  and τ = small_t_expansion_of_normalised_black_threshold .
    if x < s * ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD
        && 0.5 * s * s + x < s * (SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD + ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD)
    {
        return asymptotic_expansion_of_normalised_black_call(x / s, 0.5 * s);
    }
    if 0.5 * s < SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD {
        return small_t_expansion_of_normalised_black_call(x / s, 0.5 * s);
    }

    normalised_black_call_with_optimal_use_of_codys_functions(x, s)
}

pub(crate) fn normalised_black(x: f64, s: f64, q: f64) -> f64 {
    normalised_black_call(if q < 0.0 { -x } else { x }, s) /* Reciprocal-strike call-put equivalence */
}