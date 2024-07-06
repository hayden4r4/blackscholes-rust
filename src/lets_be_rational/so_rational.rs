use crate::lets_be_rational::{DENORMALISATION_CUTOFF, ONE_OVER_SQRT_TWO_PI};
use crate::lets_be_rational::black::normalised_black_call;
use crate::lets_be_rational::intrinsic::normalised_intrinsic;
use crate::lets_be_rational::normal_distribution::{inverse_f_upper_map, inverse_normal_cdf, standard_normal_cdf};
use crate::lets_be_rational::rational_cubic::{convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side, convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side, rational_cubic_interpolation};
use crate::OptionType;

const VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM: f64 = f64::MAX;

const VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC: f64 = -f64::MAX;


const SQRT_DBL_MIN: f64 = 1.4916681462400413e-154;
const SQRT_DBL_MAX: f64 = 1.340_780_792_994_259_6e154;

const SQRT_THREE: f64 = 1.732_050_807_568_877_2; //3.0_f64.sqrt();

const TWO_PI_OVER_SQRT_TWENTY_SEVEN: f64 = 1.209_199_576_156_145_2; // f64::PI() * 2.0 / sqrt(27.0);

const SQRT_PI_OVER_TWO: f64 = 1.253_314_137_315_500_3; //sqrt(f64::PI() / 2.0_f64);

const SQRT_ONE_OVER_THREE: f64 = 0.577_350_269_189_625_7; // sqrt(1.0 / 3.0_f64);

const PI_OVER_SIX: f64 = std::f64::consts::PI / 6.0;


fn is_below_horizon(x: f64) -> bool {
    x.abs() < DENORMALISATION_CUTOFF
}


pub(crate) fn asymptotic_expansion_of_normalised_black_call(h: f64, t: f64) -> f64 {
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


pub(crate) fn normalised_vega(x: f64, s: f64) -> f64 {
    let ax = x.abs();
    if ax <= 0.0 {
        ONE_OVER_SQRT_TWO_PI * (-0.125 * s * s).exp()
    } else if s <= 0.0 || s <= ax * SQRT_DBL_MIN {
        0.0
    } else {
        ONE_OVER_SQRT_TWO_PI * (-0.5 * ((x / s).powi(2) + (0.5 * s).powi(2))).exp()
    }
}

fn householder_factor(newton: f64, halley: f64, hh3: f64) -> f64 {
    (1.0 + 0.5 * halley * newton) / (1.0 + newton * (halley + hh3 * newton / 6.0))
}

fn inverse_f_lower_map(x: f64, f: f64) -> f64 {
    if is_below_horizon(f) {
        0.0
    } else {
        (x / (SQRT_THREE * inverse_normal_cdf((f / (TWO_PI_OVER_SQRT_TWENTY_SEVEN * x.abs())).powf(1.0 / 3.0)))).abs()
    }
}

fn compute_f_upper_map_and_first_two_derivatives(x: f64, s: f64) -> (f64, f64, f64) {
    let f = standard_normal_cdf(-0.5 * s);
    let (fp, fpp);

    if is_below_horizon(x) {
        fp = -0.5;
        fpp = 0.0;
    } else {
        let w = (x / s).powi(2);
        fp = -0.5 * (0.5 * w).exp();
        fpp = SQRT_PI_OVER_TWO * (w + 0.125 * s * s).exp() * w / s;
    }

    (f, fp, fpp)
}

fn compute_f_lower_map_and_first_two_derivatives(x: f64, s: f64) -> (f64, f64, f64) {
    let ax = x.abs();
    let z = SQRT_ONE_OVER_THREE * ax / s;
    let y = z * z;
    let s2 = s * s;
    let phi = standard_normal_cdf(z);
    let phi_ = standard_normal_cdf(-z);
    let fpp = PI_OVER_SIX * y / (s2 * s) * phi_ * (8.0 * SQRT_THREE * s * ax + (3.0 * s2 * (s2 - 8.0) - 8.0 * x * x) * phi_ / phi) * (2.0 * y + 0.25 * s2).exp();

    let (fp, f);

    if is_below_horizon(s) {
        fp = 1.0;
        f = 0.0;
    } else {
        let phi2 = phi_ * phi_;
        fp = std::f64::consts::PI * 2.0 * y * phi2 * (y + 0.125 * s * s).exp();
        f = if is_below_horizon(x) {
            0.0
        } else {
            TWO_PI_OVER_SQRT_TWENTY_SEVEN * ax * (phi2 * phi_)
        };
    }

    (f, fp, fpp)
}

pub(crate) fn implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(
    market_price: f64, forward_price: f64, strike_price: f64, time_to_maturity: f64, mut option_type: OptionType, max_iteration: i32,
) -> f64 {
    let q = option_type as i32 as f64;
    let mut price = market_price;
    let intrinsic = (q as i32 as f64 * (strike_price - forward_price)).max(0.0).abs();
    if price < intrinsic {
        return VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC;
    }
    let max_price = if q < 0.0 { strike_price } else { forward_price };
    if price >= max_price {
        return VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM;
    }
    let x = (forward_price / strike_price).ln();
    if q as i32 as f64 * x > 0.0 {
        price = (price - intrinsic).max(0.0).abs();
        option_type = -option_type;
    }
    unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(
        price / (forward_price.sqrt() * strike_price.sqrt()), x, option_type, max_iteration,
    ) / time_to_maturity.sqrt()
}


fn normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(beta: f64, x: f64, option_type: OptionType, max_iteration: i32) -> f64 {
    // Map in-the-money to out-of-the-money
    let mut beta = beta;
    let mut q = option_type;
    if q as i32 as f64 * x > 0.0 {
        beta -= normalised_intrinsic(x, q);
        q = -q;
    }
    if beta < 0.0 {
        return VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC;
    }

    unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(beta, x, q, max_iteration)
}

pub(crate) fn unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(
    mut beta: f64, mut x: f64, option_type: OptionType, n: i32,
) -> f64 {
    if option_type as i32 as f64 * x > 0.0 {
        beta = (beta - normalised_intrinsic(x, option_type)).abs().max(0.0);
    }
    if option_type == OptionType::Put {
        x = -x;
    }
    if beta <= 0.0 {
        return 0.0;
    }
    if beta < DENORMALISATION_CUTOFF {
        return 0.0;
    }
    let b_max = (0.5 * x).exp();
    if beta >= b_max {
        return VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM;
    }

    let iterations = 0;
    let mut direction_reversal_count = 0;
    let mut f = -f64::MAX;
    let mut s = -f64::MAX;
    let mut ds = s;
    let mut ds_previous = 0.0;
    let mut s_left = f64::MIN;
    let mut s_right = f64::MAX;

    let s_c = (2.0 * x.abs()).sqrt();
    let b_c = normalised_black_call(x, s_c);
    let v_c = normalised_vega(x, s_c);

    if beta < b_c {
        let s_l = s_c - b_c / v_c;
        let b_l = normalised_black_call(x, s_l);
        if beta < b_l {
            let (f_lower_map_l, d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2) = compute_f_lower_map_and_first_two_derivatives(x, s_l);
            let r_ll = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(0.0, b_l, 0.0, f_lower_map_l, 1.0, d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2, true);
            f = rational_cubic_interpolation(beta, 0.0, b_l, 0.0, f_lower_map_l, 1.0, d_f_lower_map_l_d_beta, r_ll);
            if !(f > 0.0) {
                let t = beta / b_l;
                f = (f_lower_map_l * t + b_l * (1.0 - t)) * t;
            }
            s = inverse_f_lower_map(x, f);
            s_right = s_l;

            for _ in 0..n {
                if ds.abs() <= f64::EPSILON * s {
                    break;
                }
                if ds * ds_previous < 0.0 {
                    direction_reversal_count += 1;
                }
                if iterations > 0 && (direction_reversal_count == 3 || !(s > s_left && s < s_right)) {
                    s = 0.5 * (s_left + s_right);
                    if s_right - s_left <= f64::EPSILON * s {
                        break;
                    }
                    direction_reversal_count = 0;
                    ds = 0.0;
                }
                ds_previous = ds;
                let b = normalised_black_call(x, s);
                let bp = normalised_vega(x, s);
                if b > beta && s < s_right {
                    s_right = s;
                } else if b < beta && s > s_left {
                    s_left = s;
                }
                if b <= 0.0 || bp <= 0.0 {
                    ds = 0.5 * (s_left + s_right) - s;
                } else {
                    let ln_b = b.ln();
                    let ln_beta = beta.ln();
                    let bpob = bp / b;
                    let h = x / s;
                    let b_halley = h * h / s - s / 4.0;
                    let newton = (ln_beta - ln_b) * ln_b / ln_beta / bpob;
                    let halley = b_halley - bpob * (1.0 + 2.0 / ln_b);
                    let b_hh3 = b_halley * b_halley - 3.0 * (h / s).powi(2) - 0.25;
                    let hh3 = b_hh3 + 2.0 * bpob.powi(2) * (1.0 + 3.0 / ln_b * (1.0 + 1.0 / ln_b)) - 3.0 * b_halley * bpob * (1.0 + 2.0 / ln_b);
                    ds = newton * householder_factor(newton, halley, hh3);
                }
                s += ds.max(-0.5 * s);
            }
            return s;
        } else {
            let v_l = normalised_vega(x, s_l);
            let r_lm = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(b_l, b_c, s_l, s_c, 1.0 / v_l, 1.0 / v_c, 0.0, false);
            s = rational_cubic_interpolation(beta, b_l, b_c, s_l, s_c, 1.0 / v_l, 1.0 / v_c, r_lm);
            s_left = s_l;
            s_right = s_c;
        }
    } else {
        let s_h = if v_c > f64::MIN { s_c + (b_max - b_c) / v_c } else { s_c };
        let b_h = normalised_black_call(x, s_h);
        if beta <= b_h {
            let v_h = normalised_vega(x, s_h);
            let r_hm = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(b_c, b_h, s_c, s_h, 1.0 / v_c, 1.0 / v_h, 0.0, false);
            s = rational_cubic_interpolation(beta, b_c, b_h, s_c, s_h, 1.0 / v_c, 1.0 / v_h, r_hm);
            s_left = s_c;
            s_right = s_h;
        } else {
            let (f_upper_map_h, d_f_upper_map_h_d_beta, d2_f_upper_map_h_d_beta2) = compute_f_upper_map_and_first_two_derivatives(x, s_h);
            if d2_f_upper_map_h_d_beta2 > -SQRT_DBL_MAX && d2_f_upper_map_h_d_beta2 < SQRT_DBL_MAX {
                let r_hh = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(b_h, b_max, f_upper_map_h, 0.0, d_f_upper_map_h_d_beta, -0.5, d2_f_upper_map_h_d_beta2, true);
                f = rational_cubic_interpolation(beta, b_h, b_max, f_upper_map_h, 0.0, d_f_upper_map_h_d_beta, -0.5, r_hh);
            }
            if f <= 0.0 {
                let h = b_max - b_h;
                let t = (beta - b_h) / h;
                f = (f_upper_map_h * (1.0 - t) + 0.5 * h * t) * (1.0 - t);
            }
            s = inverse_f_upper_map(f);
            s_left = s_h;
            if beta > 0.5 * b_max {
                for _ in 0..n {
                    if ds.abs() <= f64::EPSILON * s {
                        break;
                    }
                    if ds * ds_previous < 0.0 {
                        direction_reversal_count += 1;
                    }
                    if iterations > 0 && (direction_reversal_count == 3 || !(s > s_left && s < s_right)) {
                        s = 0.5 * (s_left + s_right);
                        if s_right - s_left <= f64::EPSILON * s {
                            break;
                        }
                        direction_reversal_count = 0;
                        ds = 0.0;
                    }
                    ds_previous = ds;
                    let b = normalised_black_call(x, s);
                    let bp = normalised_vega(x, s);
                    if b > beta && s < s_right {
                        s_right = s;
                    } else if b < beta && s > s_left {
                        s_left = s;
                    }
                    if b >= b_max || bp <= f64::MIN {
                        ds = 0.5 * (s_left + s_right) - s;
                    } else {
                        let b_max_minus_b = b_max - b;
                        let g = ((b_max - beta) / b_max_minus_b).ln();
                        let gp = bp / b_max_minus_b;
                        let b_halley = (x / s).powi(2) / s - s / 4.0;
                        let b_hh3 = b_halley * b_halley - 3.0 * (x / (s * s)).powi(2) - 0.25;
                        let newton = -g / gp;
                        let halley = b_halley + gp;
                        let hh3 = b_hh3 + gp * (2.0 * gp + 3.0 * b_halley);
                        ds = newton * householder_factor(newton, halley, hh3);
                    }
                    s += ds.max(-0.5 * s);
                }
                return s;
            }
        }
    }

    for _ in 0..n {
        if ds.abs() <= f64::EPSILON * s {
            break;
        }
        if ds * ds_previous < 0.0 {
            direction_reversal_count += 1;
        }
        if iterations > 0 && (direction_reversal_count == 3 || !(s > s_left && s < s_right)) {
            s = 0.5 * (s_left + s_right);
            if s_right - s_left <= f64::EPSILON * s {
                break;
            }
            direction_reversal_count = 0;
            ds = 0.0;
        }
        ds_previous = ds;
        let b = normalised_black_call(x, s);
        let bp = normalised_vega(x, s);
        if b > beta && s < s_right {
            s_right = s;
        } else if b < beta && s > s_left {
            s_left = s;
        }
        let newton = (beta - b) / bp;
        let halley = (x / s).powi(2) / s - s / 4.0;
        let hh3 = halley * halley - 3.0 * (x / (s * s)).powi(2) - 0.25;
        ds = newton * householder_factor(newton, halley, hh3).max(-0.5 * s);
        s += ds;
    }

    s
}
