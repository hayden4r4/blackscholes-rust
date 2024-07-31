use std::f64::consts::PI;

use crate::{
    lets_be_rational::{
        black::normalised_black_call,
        intrinsic::normalised_intrinsic,
        normal_distribution::{inverse_f_upper_map, inverse_normal_cdf, standard_normal_cdf},
        rational_cubic::{
            convex_rational_cubic_control_parameter, rational_cubic_interpolation, Side,
        },
        DENORMALISATION_CUTOFF, ONE_OVER_SQRT_TWO_PI,
    },
    OptionType,
};

const VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM: f64 = f64::MAX;

const VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC: f64 = -f64::MAX;

const SQRT_DBL_MIN: f64 = 1.4916681462400413e-154;
const SQRT_DBL_MAX: f64 = 1.340_780_792_994_259_6e154;

const SQRT_THREE: f64 = 1.732_050_807_568_877_2; //3.0_f64.sqrt();

const TWO_PI_OVER_SQRT_TWENTY_SEVEN: f64 = 1.209_199_576_156_145_2; // f64::PI() * 2.0 / sqrt(27.0);

const SQRT_PI_OVER_TWO: f64 = 1.253_314_137_315_500_3; //sqrt(f64::PI() / 2.0_f64);

const SQRT_ONE_OVER_THREE: f64 = 0.577_350_269_189_625_7; // sqrt(1.0 / 3.0_f64);

const PI_OVER_SIX: f64 = PI / 6.0;

fn is_below_horizon(x: f64) -> bool {
    x.abs() < DENORMALISATION_CUTOFF
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
        (x / (SQRT_THREE
            * inverse_normal_cdf((f / (TWO_PI_OVER_SQRT_TWENTY_SEVEN * x.abs())).powf(1.0 / 3.0))))
        .abs()
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
    let fpp = PI_OVER_SIX * y / (s2 * s)
        * phi_
        * (8.0 * SQRT_THREE * s * ax + (3.0 * s2 * (s2 - 8.0) - 8.0 * x * x) * phi_ / phi)
        * (2.0 * y + 0.25 * s2).exp();

    let (fp, f);

    if is_below_horizon(s) {
        fp = 1.0;
        f = 0.0;
    } else {
        let phi2 = phi_ * phi_;
        fp = PI * 2.0 * y * phi2 * (y + 0.125 * s * s).exp();
        f = if is_below_horizon(x) {
            0.0
        } else {
            TWO_PI_OVER_SQRT_TWENTY_SEVEN * ax * (phi2 * phi_)
        };
    }

    (f, fp, fpp)
}

pub(crate) fn implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(
    market_price: f64,
    forward_price: f64,
    strike_price: f64,
    time_to_maturity: f64,
    option_type: OptionType,
    max_iteration: i32,
) -> f64 {
    let mut price = market_price;
    let intrinsic = (option_type * (forward_price - strike_price)).max(0.0);
    if price < intrinsic {
        return VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC;
    }
    let max_price = match option_type {
        OptionType::Call => strike_price,
        OptionType::Put => forward_price,
    };
    if price >= max_price {
        return VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM;
    }
    let x = (forward_price / strike_price).ln();
    let option_type = if option_type * x > 0.0 {
        price = (price - intrinsic).max(0.0);
        -option_type
    } else {
        option_type
    };
    unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(
        price / (forward_price.sqrt() * strike_price.sqrt()), x, option_type, max_iteration,
    ) / time_to_maturity.sqrt()
}

#[allow(dead_code)]
fn normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(
    beta: f64,
    x: f64,
    option_type: OptionType,
    max_iteration: i32,
) -> f64 {
    // Map in-the-money to out-of-the-money
    let mut beta = beta;
    let option_type = if option_type * x > 0.0 {
        beta -= normalised_intrinsic(x, option_type);
        -option_type
    } else {
        option_type
    };
    if beta < 0.0 {
        return VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC;
    }

    unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(beta, x, option_type, max_iteration)
}

pub(crate) fn unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(
    beta: f64,
    x: f64,
    option_type: OptionType,
    n: i32,
) -> f64 {
    let beta = if option_type * x > 0.0 {
        (beta - normalised_intrinsic(x, option_type)).max(0.0)
    } else {
        beta
    };
    let x = option_type * x;
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
            let (f_lower_map_l, d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2) =
                compute_f_lower_map_and_first_two_derivatives(x, s_l);
            let r_ll = convex_rational_cubic_control_parameter(
                0.0,
                b_l,
                0.0,
                f_lower_map_l,
                1.0,
                d_f_lower_map_l_d_beta,
                d2_f_lower_map_l_d_beta2,
                true,
                Side::Right,
            );
            // TODO: Expect terrible approach, handle it properly
            f = rational_cubic_interpolation(
                beta,
                0.0,
                b_l,
                0.0,
                f_lower_map_l,
                1.0,
                d_f_lower_map_l_d_beta,
                r_ll,
            )
            .expect("We should expect correct parameters");
            if f <= 0.0 {
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
                if iterations > 0 && (direction_reversal_count == 3 || !(s > s_left && s < s_right))
                {
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
                    let hh3 = b_hh3 + 2.0 * bpob.powi(2) * (1.0 + 3.0 / ln_b * (1.0 + 1.0 / ln_b))
                        - 3.0 * b_halley * bpob * (1.0 + 2.0 / ln_b);
                    ds = newton * householder_factor(newton, halley, hh3);
                }
                s += ds.max(-0.5 * s);
            }
            return s;
        } else {
            let v_l = normalised_vega(x, s_l);
            let r_lm = convex_rational_cubic_control_parameter(
                b_l,
                b_c,
                s_l,
                s_c,
                1.0 / v_l,
                1.0 / v_c,
                0.0,
                false,
                Side::Right,
            );
            // TODO: Expect terrible approach, handle it properly
            s = rational_cubic_interpolation(beta, b_l, b_c, s_l, s_c, 1.0 / v_l, 1.0 / v_c, r_lm)
                .expect("We should expect correct parameters");
            s_left = s_l;
            s_right = s_c;
        }
    } else {
        let s_h = if v_c > f64::MIN {
            s_c + (b_max - b_c) / v_c
        } else {
            s_c
        };
        let b_h = normalised_black_call(x, s_h);
        if beta <= b_h {
            let v_h = normalised_vega(x, s_h);
            let r_hm = convex_rational_cubic_control_parameter(
                b_c,
                b_h,
                s_c,
                s_h,
                1.0 / v_c,
                1.0 / v_h,
                0.0,
                false,
                Side::Left,
            );
            // TODO: Expect terrible approach, handle it properly
            s = rational_cubic_interpolation(beta, b_c, b_h, s_c, s_h, 1.0 / v_c, 1.0 / v_h, r_hm)
                .expect("We should expect correct parameters");
            s_left = s_c;
            s_right = s_h;
        } else {
            let (f_upper_map_h, d_f_upper_map_h_d_beta, d2_f_upper_map_h_d_beta2) =
                compute_f_upper_map_and_first_two_derivatives(x, s_h);
            if d2_f_upper_map_h_d_beta2 > -SQRT_DBL_MAX && d2_f_upper_map_h_d_beta2 < SQRT_DBL_MAX {
                let r_hh = convex_rational_cubic_control_parameter(
                    b_h,
                    b_max,
                    f_upper_map_h,
                    0.0,
                    d_f_upper_map_h_d_beta,
                    -0.5,
                    d2_f_upper_map_h_d_beta2,
                    true,
                    Side::Left,
                );
                // TODO: Expect terrible approach, handle it properly
                f = rational_cubic_interpolation(
                    beta,
                    b_h,
                    b_max,
                    f_upper_map_h,
                    0.0,
                    d_f_upper_map_h_d_beta,
                    -0.5,
                    r_hh,
                )
                .expect("We should expect correct parameters");
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
                    if iterations > 0
                        && (direction_reversal_count == 3 || !(s > s_left && s < s_right))
                    {
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
