use std::f64::consts::FRAC_1_SQRT_2;

use statrs::function::erf::erfc;

use crate::lets_be_rational::{DENORMALISATION_CUTOFF, intrinsic::normalised_intrinsic, normal_distribution::standard_normal_cdf, ONE_OVER_SQRT_TWO_PI, SQRT_TWO_PI};
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


pub fn asymptotic_expansion_of_normalised_black_call(h: f64, t: f64) -> f64 {
    let e = (t / h) * (t / h);
    let r = (h + t) * (h - t);
    let q = (h / r) * (h / r);

    // 17th order asymptotic expansion of A(h,t) in q
    let asymptotic_expansion_sum = 2.0 + q * (-6.0 - 2.0 * e + 3.0 * q * (10.0 + e * (20.0 + 2.0 * e)
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