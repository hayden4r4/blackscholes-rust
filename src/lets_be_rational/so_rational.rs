use crate::lets_be_rational::black::normalised_black_call;
use crate::lets_be_rational::intrinsic::normalised_intrinsic;
use crate::lets_be_rational::normal_distribution::{
    inverse_f_upper_map, inverse_normal_cdf, standard_normal_cdf,
};
use crate::lets_be_rational::rational_cubic::{
    convex_rational_cubic_control_parameter, rational_cubic_interpolation, Side,
};
use crate::lets_be_rational::{DENORMALISATION_CUTOFF, ONE_OVER_SQRT_TWO_PI};
use crate::OptionType;

use num_traits::{AsPrimitive, Float, FloatConst, FromPrimitive};

const VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM: f64 = f64::MAX;

const VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC: f64 = -f64::MAX;

const SQRT_DBL_MIN: f64 = 1.4916681462400413e-154;
const SQRT_DBL_MAX: f64 = 1.340_780_792_994_259_6e154;

const SQRT_THREE: f64 = 1.732_050_807_568_877_2; //3.0_f64.sqrt();

const TWO_PI_OVER_SQRT_TWENTY_SEVEN: f64 = 1.209_199_576_156_145_2; // f64::PI() * 2.0 / sqrt(27.0);

const SQRT_PI_OVER_TWO: f64 = 1.253_314_137_315_500_3; //sqrt(f64::PI() / 2.0_f64);

const SQRT_ONE_OVER_THREE: f64 = 0.577_350_269_189_625_7; // sqrt(1.0 / 3.0_f64);

const PI_OVER_SIX: f64 = std::f64::consts::PI / 6.0;

fn is_below_horizon<T>(x: T) -> bool
where
    T: Float + FromPrimitive,
{
    x.abs() < T::from(DENORMALISATION_CUTOFF).unwrap()
}

pub(crate) fn normalised_vega<T>(x: T, s: T) -> T
where
    T: Float + FromPrimitive,
{
    let ax = x.abs();
    if ax <= T::zero() {
        T::from(ONE_OVER_SQRT_TWO_PI).unwrap() * (T::from(-0.125).unwrap() * s * s).exp()
    } else if s <= T::zero() || s <= ax * T::from(SQRT_DBL_MIN).unwrap() {
        T::zero()
    } else {
        T::from(ONE_OVER_SQRT_TWO_PI).unwrap()
            * (T::from(-0.5).unwrap() * ((x / s).powi(2) + (T::from(0.5).unwrap() * s).powi(2)))
                .exp()
    }
}

fn householder_factor<T>(newton: T, halley: T, hh3: T) -> T
where
    T: Float + FromPrimitive,
{
    (T::one() + T::from(0.5).unwrap() * halley * newton)
        / (T::one() + newton * (halley + hh3 * newton / T::from(6.0).unwrap()))
}

fn inverse_f_lower_map<T>(x: T, f: T) -> T
where
    T: Float + FromPrimitive + AsPrimitive<f64>,
{
    if is_below_horizon(f) {
        T::zero()
    } else {
        (x / (T::from(SQRT_THREE).unwrap()
            * T::from(inverse_normal_cdf(
                (f.as_() / (TWO_PI_OVER_SQRT_TWENTY_SEVEN * x.as_().abs())).powf(1.0 / 3.0),
            ))
            .unwrap()))
        .abs()
    }
}

fn compute_f_upper_map_and_first_two_derivatives<T>(x: T, s: T) -> (T, T, T)
where
    T: Float + FromPrimitive + AsPrimitive<f64>,
{
    let f = T::from(standard_normal_cdf(-0.5 * s.as_())).unwrap();
    let (fp, fpp);

    if is_below_horizon(x) {
        fp = T::from(-0.5).unwrap();
        fpp = T::zero();
    } else {
        let w = (x / s).powi(2);
        fp = T::from(-0.5).unwrap() * (T::from(0.5).unwrap() * w).exp();
        fpp = T::from(SQRT_PI_OVER_TWO).unwrap() * (w + T::from(0.125).unwrap() * s * s).exp() * w
            / s;
    }

    (f, fp, fpp)
}

fn compute_f_lower_map_and_first_two_derivatives<T>(x: T, s: T) -> (T, T, T)
where
    T: Float + FromPrimitive + AsPrimitive<f64> + FloatConst,
{
    let ax = x.abs();
    let z = T::from(SQRT_ONE_OVER_THREE).unwrap() * ax / s;
    let y = z * z;
    let s2 = s * s;
    let phi = T::from(standard_normal_cdf(z.as_())).unwrap();
    let phi_ = T::from(standard_normal_cdf(-z.as_())).unwrap();
    let fpp = T::from(PI_OVER_SIX).unwrap() * y / (s2 * s)
        * phi_
        * (T::from(8.0).unwrap() * T::from(SQRT_THREE).unwrap() * s * ax
            + (T::from(3.0).unwrap() * s2 * (s2 - T::from(8.0).unwrap())
                - T::from(8.0).unwrap() * x * x)
                * phi_
                / phi)
        * (T::from(2.0).unwrap() * y + T::from(0.25).unwrap() * s2).exp();

    let (fp, f);

    if is_below_horizon(s) {
        fp = T::one();
        f = T::zero();
    } else {
        let phi2 = phi_ * phi_;
        fp = T::PI()
            * T::from(2.0).unwrap()
            * y
            * phi2
            * (y + T::from(0.125).unwrap() * s * s).exp();
        f = if is_below_horizon(x) {
            T::zero()
        } else {
            T::from(TWO_PI_OVER_SQRT_TWENTY_SEVEN).unwrap() * ax * (phi2 * phi_)
        };
    }

    (f, fp, fpp)
}

pub(crate) fn implied_volatility_from_a_transformed_rational_guess_with_limited_iterations<T>(
    market_price: T,
    forward_price: T,
    strike_price: T,
    time_to_maturity: T,
    mut option_type: OptionType,
    max_iteration: i32,
) -> T
where
    T: Float + FromPrimitive + AsPrimitive<f64> + FloatConst,
{
    let q = T::from(option_type as i32).unwrap();
    let mut price = market_price;
    let intrinsic = (q * (forward_price - strike_price)).max(T::zero()).abs();
    if price < intrinsic {
        return T::from(VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC).unwrap();
    }
    let max_price = if q < T::zero() {
        strike_price
    } else {
        forward_price
    };
    if price >= max_price {
        return T::from(VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM).unwrap();
    }
    let x = (forward_price / strike_price).ln();
    if q * x > T::zero() {
        price = (price - intrinsic).max(T::zero()).abs();
        option_type = -option_type;
    }
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

pub(crate) fn unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations<
    T,
>(
    mut beta: T,
    mut x: T,
    option_type: OptionType,
    n: i32,
) -> T
where
    T: Float + FromPrimitive + AsPrimitive<f64> + FloatConst,
{
    if T::from(option_type as i32).unwrap() * x > T::zero() {
        beta = (beta - normalised_intrinsic(x, option_type))
            .abs()
            .max(T::zero());
    }
    if option_type == OptionType::Put {
        x = -x;
    }
    if beta <= T::zero() {
        return T::zero();
    }
    if beta < T::from(DENORMALISATION_CUTOFF).unwrap() {
        return T::zero();
    }
    let b_max = (T::from(0.5).unwrap() * x).exp();
    if beta >= b_max {
        return T::from(VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM).unwrap();
    }

    let iterations = 0;
    let mut direction_reversal_count = 0;
    let mut f = -T::max_value();
    let mut s = -T::max_value();
    let mut ds = s;
    let mut ds_previous = T::zero();
    let mut s_left = T::min_value();
    let mut s_right = T::max_value();

    let s_c = (T::from(2.0).unwrap() * x.abs()).sqrt();
    let b_c = normalised_black_call(x, s_c);
    let v_c = normalised_vega(x, s_c);

    if beta < b_c {
        let s_l = s_c - b_c / v_c;
        let b_l = normalised_black_call(x, s_l);
        if beta < b_l {
            let (f_lower_map_l, d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2) =
                compute_f_lower_map_and_first_two_derivatives(x, s_l);
            let r_ll = convex_rational_cubic_control_parameter(
                T::zero(),
                b_l,
                T::zero(),
                f_lower_map_l,
                T::one(),
                d_f_lower_map_l_d_beta,
                d2_f_lower_map_l_d_beta2,
                true,
                Side::Right,
            );
            // TODO: Expect terrible approach, handle it properly
            f = rational_cubic_interpolation(
                beta,
                T::zero(),
                b_l,
                T::zero(),
                f_lower_map_l,
                T::one(),
                d_f_lower_map_l_d_beta,
                r_ll,
            )
            .expect("We should expect correct parameters");
            if f <= T::zero() {
                let t = beta / b_l;
                f = (f_lower_map_l * t + b_l * (T::one() - t)) * t;
            }
            s = inverse_f_lower_map(x, f);
            s_right = s_l;

            for _ in 0..n {
                if ds.abs() <= T::epsilon() * s {
                    break;
                }
                if ds * ds_previous < T::zero() {
                    direction_reversal_count += 1;
                }
                if iterations > 0 && (direction_reversal_count == 3 || !(s > s_left && s < s_right))
                {
                    s = T::from(0.5).unwrap() * (s_left + s_right);
                    if s_right - s_left <= T::epsilon() * s {
                        break;
                    }
                    direction_reversal_count = 0;
                    ds = T::zero();
                }
                ds_previous = ds;
                let b = normalised_black_call(x, s);
                let bp = normalised_vega(x, s);
                if b > beta && s < s_right {
                    s_right = s;
                } else if b < beta && s > s_left {
                    s_left = s;
                }
                if b <= T::zero() || bp <= T::zero() {
                    ds = T::from(0.5).unwrap() * (s_left + s_right) - s;
                } else {
                    let ln_b = b.ln();
                    let ln_beta = beta.ln();
                    let bpob = bp / b;
                    let h = x / s;
                    let b_halley = h * h / s - s / T::from(4.0).unwrap();
                    let newton = (ln_beta - ln_b) * ln_b / ln_beta / bpob;
                    let halley = b_halley - bpob * (T::one() + T::from(2.0).unwrap() / ln_b);
                    let b_hh3 = b_halley * b_halley
                        - T::from(3.0).unwrap() * (h / s).powi(2)
                        - T::from(0.25).unwrap();
                    let hh3 = b_hh3
                        + T::from(2.0).unwrap()
                            * bpob.powi(2)
                            * (T::one()
                                + T::from(3.0).unwrap() / ln_b * (T::one() + T::one() / ln_b))
                        - T::from(3.0).unwrap()
                            * b_halley
                            * bpob
                            * (T::one() + T::from(2.0).unwrap() / ln_b);
                    ds = newton * householder_factor(newton, halley, hh3);
                }
                s = s + ds.max(-T::from(0.5).unwrap() * s);
            }
            return s;
        } else {
            let v_l = normalised_vega(x, s_l);
            let r_lm = convex_rational_cubic_control_parameter(
                b_l,
                b_c,
                s_l,
                s_c,
                T::one() / v_l,
                T::one() / v_c,
                T::zero(),
                false,
                Side::Right,
            );
            // TODO: Expect terrible approach, handle it properly
            s = rational_cubic_interpolation(
                beta,
                b_l,
                b_c,
                s_l,
                s_c,
                T::one() / v_l,
                T::one() / v_c,
                r_lm,
            )
            .expect("We should expect correct parameters");
            s_left = s_l;
            s_right = s_c;
        }
    } else {
        let s_h = if v_c > T::min_value() {
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
                T::one() / v_c,
                T::one() / v_h,
                T::zero(),
                false,
                Side::Left,
            );
            // TODO: Expect terrible approach, handle it properly
            s = rational_cubic_interpolation(
                beta,
                b_c,
                b_h,
                s_c,
                s_h,
                T::one() / v_c,
                T::one() / v_h,
                r_hm,
            )
            .expect("We should expect correct parameters");
            s_left = s_c;
            s_right = s_h;
        } else {
            let (f_upper_map_h, d_f_upper_map_h_d_beta, d2_f_upper_map_h_d_beta2) =
                compute_f_upper_map_and_first_two_derivatives(x, s_h);
            if d2_f_upper_map_h_d_beta2 > T::from(-SQRT_DBL_MAX).unwrap()
                && d2_f_upper_map_h_d_beta2 < T::from(SQRT_DBL_MAX).unwrap()
            {
                let r_hh = convex_rational_cubic_control_parameter(
                    b_h,
                    b_max,
                    f_upper_map_h,
                    T::zero(),
                    d_f_upper_map_h_d_beta,
                    T::from(-0.5).unwrap(),
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
                    T::zero(),
                    d_f_upper_map_h_d_beta,
                    -T::from(0.5).unwrap(),
                    r_hh,
                )
                .expect("We should expect correct parameters");
            }
            if f <= T::zero() {
                let h = b_max - b_h;
                let t = (beta - b_h) / h;
                f = (f_upper_map_h * (T::one() - t) + T::from(0.5).unwrap() * h * t)
                    * (T::one() - t);
            }
            s = inverse_f_upper_map(f);
            s_left = s_h;
            if beta > T::from(0.5).unwrap() * b_max {
                for _ in 0..n {
                    if ds.abs() <= T::epsilon() * s {
                        break;
                    }
                    if ds * ds_previous < T::zero() {
                        direction_reversal_count += 1;
                    }
                    if iterations > 0
                        && (direction_reversal_count == 3 || !(s > s_left && s < s_right))
                    {
                        s = T::from(0.5).unwrap() * (s_left + s_right);
                        if s_right - s_left <= T::epsilon() * s {
                            break;
                        }
                        direction_reversal_count = 0;
                        ds = T::zero();
                    }
                    ds_previous = ds;
                    let b = normalised_black_call(x, s);
                    let bp = normalised_vega(x, s);
                    if b > beta && s < s_right {
                        s_right = s;
                    } else if b < beta && s > s_left {
                        s_left = s;
                    }
                    if b >= b_max || bp <= T::min_value() {
                        ds = T::from(0.5).unwrap() * (s_left + s_right) - s;
                    } else {
                        let b_max_minus_b = b_max - b;
                        let g = ((b_max - beta) / b_max_minus_b).ln();
                        let gp = bp / b_max_minus_b;
                        let b_halley = (x / s).powi(2) / s - s / T::from(4.0).unwrap();
                        let b_hh3 = b_halley * b_halley
                            - T::from(3.0).unwrap() * (x / (s * s)).powi(2)
                            - T::from(0.25).unwrap();
                        let newton = -g / gp;
                        let halley = b_halley + gp;
                        let hh3 = b_hh3
                            + gp * (T::from(2.0).unwrap() * gp + T::from(3.0).unwrap() * b_halley);
                        ds = newton * householder_factor(newton, halley, hh3);
                    }
                    s = s + ds.max(-T::from(0.5).unwrap() * s);
                }
                return s;
            }
        }
    }

    for _ in 0..n {
        if ds.abs() <= T::epsilon() * s {
            break;
        }
        if ds * ds_previous < T::zero() {
            direction_reversal_count += 1;
        }
        if iterations > 0 && (direction_reversal_count == 3 || !(s > s_left && s < s_right)) {
            s = T::from(0.5).unwrap() * (s_left + s_right);
            if s_right - s_left <= T::epsilon() * s {
                break;
            }
            direction_reversal_count = 0;
            ds = T::zero();
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
        let halley = (x / s).powi(2) / s - s / T::from(4.0).unwrap();
        let hh3 = halley * halley
            - T::from(3.0).unwrap() * (x / (s * s)).powi(2)
            - T::from(0.25).unwrap();
        ds = newton * householder_factor(newton, halley, hh3).max(-T::from(0.5).unwrap() * s);
        s = s + ds;
    }

    s
}
