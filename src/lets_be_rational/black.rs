use num_traits::{AsPrimitive, Float, FromPrimitive};
use std::f64::consts::FRAC_1_SQRT_2;

use statrs::function::erf::erfc;

use crate::lets_be_rational::{
    intrinsic::normalised_intrinsic, normal_distribution::standard_normal_cdf,
    DENORMALISATION_CUTOFF, ONE_OVER_SQRT_TWO_PI, SQRT_TWO_PI,
};
use crate::OptionType;

const ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD: f64 = -10.0;
const SIXTEENTH_ROOT_DBL_EPSILON: f64 = 0.10566243270259357;
const CODYS_THRESHOLD: f64 = 0.46875;

const SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD: f64 = 2.0 * SIXTEENTH_ROOT_DBL_EPSILON;

fn erfcx<T>(x: T) -> T
where
    T: Float + FromPrimitive + AsPrimitive<f64>,
{
    (x * x).exp() * T::from(erfc(x.as_())).unwrap()
}

#[allow(dead_code)]
fn normalised_black_call_using_norm_cdf(x: f64, s: f64) -> f64 {
    let h = x / s;
    let t = 0.5 * s;
    let b_max = (0.5 * x).exp();
    let b = standard_normal_cdf(h + t) * b_max - standard_normal_cdf(h - t) / b_max;
    b.abs().max(0.0)
}

fn normalised_black_call_with_optimal_use_of_codys_functions<T>(x: T, s: T) -> T
where
    T: Float + FromPrimitive + AsPrimitive<f64>,
{
    let h = x / s;
    let half = T::from(0.5).unwrap();
    let t = half * s;
    let q1 = -T::from(FRAC_1_SQRT_2).unwrap() * (h + t);
    let q2 = -T::from(FRAC_1_SQRT_2).unwrap() * (h - t);
    let two_b;

    if q1 < T::from(CODYS_THRESHOLD).unwrap() {
        if q2 < T::from(CODYS_THRESHOLD).unwrap() {
            two_b = (half * x).exp() * T::from(erfc(q1.as_())).unwrap()
                - (-half * x).exp() * T::from(erfc(q2.as_())).unwrap();
        } else {
            two_b = (half * x).exp() * T::from(erfc(q1.as_())).unwrap()
                - (-half * (h * h + t * t)).exp() * erfcx(q2);
        }
    } else if q2 < T::from(CODYS_THRESHOLD).unwrap() {
        two_b = (-half * (h * h + t * t)).exp() * erfcx(q1)
            - (-half * x).exp() * T::from(erfc(q2.as_())).unwrap();
    } else {
        two_b = (-half * (h * h + t * t)).exp() * (erfcx(q1) - erfcx(q2));
    }

    (half * two_b).abs().max(T::from(0.0).unwrap())
}

#[rustfmt::skip]
fn small_t_expansion_of_normalised_black_call<T>(h: T, t: T) -> T 
where
    T: Float + FromPrimitive + AsPrimitive<f64>,
{
    let a = T::one() + h * (T::from(0.5).unwrap() * T::from(SQRT_TWO_PI).unwrap()) * erfcx(-T::from(FRAC_1_SQRT_2).unwrap() * h);
    let w = t * t;
    let h2 = h * h;

    let expansion = T::from(2.0).unwrap() * t * (a + w * ((-T::from(1.0).unwrap() + T::from(3.0).unwrap() * a + a * h2) / T::from(6.0).unwrap()
        + w * ((-T::from(7.0).unwrap() + T::from(15.0).unwrap() * a + h2 * (-T::from(1.0).unwrap() + T::from(10.0).unwrap() * a + a * h2)) / T::from(120.0).unwrap()
        + w * ((-T::from(57.0).unwrap() + T::from(105.0).unwrap() * a + h2 * (-T::from(18.0).unwrap() + T::from(105.0).unwrap() * a + h2 * (-T::from(1.0).unwrap() + T::from(21.0).unwrap() * a + a * h2))) / T::from(5040.0).unwrap()
        + w * ((-T::from(561.0).unwrap() + T::from(945.0).unwrap() * a + h2 * (-T::from(285.0).unwrap() + T::from(1260.0).unwrap() * a + h2 * (-T::from(33.0).unwrap() + T::from(378.0).unwrap() * a + h2 * (-T::from(1.0).unwrap() + T::from(36.0).unwrap() * a + a * h2)))) / T::from(362880.0).unwrap()
        + w * ((-T::from(6555.0).unwrap() + T::from(10395.0).unwrap() * a + h2 * (-T::from(4680.0).unwrap() + T::from(17325.0).unwrap() * a + h2 * (-T::from(840.0).unwrap() + T::from(6930.0).unwrap() * a + h2 * (-T::from(52.0).unwrap() + T::from(990.0).unwrap() * a + h2 * (-T::from(1.0).unwrap() + T::from(55.0).unwrap() * a + a * h2))))) / T::from(39916800.0).unwrap()
        + w * ((-T::from(89055.0).unwrap() + T::from(135135.0).unwrap() * a + h2 * (-T::from(82845.0).unwrap() + T::from(270270.0).unwrap() * a + h2 * (-T::from(20370.0).unwrap() + T::from(135135.0).unwrap() * a + h2 * (-T::from(1926.0).unwrap() + T::from(25740.0).unwrap() * a + h2 * (-T::from(75.0).unwrap() + T::from(2145.0).unwrap() * a + h2 * (-T::from(1.0).unwrap() + T::from(78.0).unwrap() * a + a * h2)))))) / T::from(6227020800.0).unwrap())))))));

    let b = T::from(ONE_OVER_SQRT_TWO_PI).unwrap() * (-T::from(0.5).unwrap() * (h * h + t * t)).exp() * expansion;
    b.abs().max(T::zero())
}

#[allow(dead_code)]
pub fn normalised_black_call_using_erfcx(h: f64, t: f64) -> f64 {
    let b = 0.5
        * (-0.5 * (h * h + t * t)).exp()
        * (erfcx(-FRAC_1_SQRT_2 * (h + t)) - erfcx(-FRAC_1_SQRT_2 * (h - t)));
    b.abs().max(0.0)
}

pub(crate) fn normalised_black_call<T>(x: T, s: T) -> T
where
    T: Float + FromPrimitive + AsPrimitive<f64>,
{
    if x > T::zero() {
        return normalised_intrinsic(x, OptionType::Call) + normalised_black_call(-x, s);
        // In the money.
    }
    if s <= x.abs() * T::from(DENORMALISATION_CUTOFF).unwrap() {
        return normalised_intrinsic(x, OptionType::Call); // sigma=0 -> intrinsic value.
    }

    // Denote h := x/s and t := s/2.
    // We evaluate the condition |h|>|η|, i.e., h<η  &&  t < τ+|h|-|η|  avoiding any divisions by s , where η = asymptotic_expansion_accuracy_threshold  and τ = small_t_expansion_of_normalised_black_threshold .
    if x < s * T::from(ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD).unwrap()
        && T::from(0.5).unwrap() * s * s + x
            < s * (T::from(SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD).unwrap()
                + T::from(ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD).unwrap())
    {
        return asymptotic_expansion_of_normalised_black_call(x / s, T::from(0.5).unwrap() * s);
    }
    if T::from(0.5).unwrap() * s < T::from(SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD).unwrap()
    {
        return small_t_expansion_of_normalised_black_call(x / s, T::from(0.5).unwrap() * s);
    }

    normalised_black_call_with_optimal_use_of_codys_functions(x, s)
}

pub(crate) fn normalised_black<T>(x: T, s: T, q: T) -> T
where
    T: Float + FromPrimitive + AsPrimitive<f64>,
{
    normalised_black_call(if q < T::zero() { -x } else { x }, s) /* Reciprocal-strike call-put equivalence */
}

#[rustfmt::skip]
pub fn asymptotic_expansion_of_normalised_black_call<T>(h: T, t: T) -> T 
where
    T: Float + FromPrimitive,
{
    let e = (t / h) * (t / h);
    let r = (h + t) * (h - t);
    let q = (h / r) * (h / r);

    // 17th order asymptotic expansion of A(h,t) in q
    let asymptotic_expansion_sum = T::from(2.0).unwrap() + q * (-T::from(6.0).unwrap() - T::from(2.0).unwrap() * e + T::from(3.0).unwrap() * q * (T::from(10.0).unwrap() + e * (T::from(20.0).unwrap() + T::from(2.0).unwrap() * e)
        + T::from(5.0).unwrap() * q * (-T::from(14.0).unwrap() + e * (-T::from(70.0).unwrap() + e * (-T::from(42.0).unwrap() - T::from(2.0).unwrap() * e))
        + T::from(7.0).unwrap() * q * (T::from(18.0).unwrap() + e * (T::from(168.0).unwrap() + e * (T::from(252.0).unwrap() + e * (T::from(72.0).unwrap() + T::from(2.0).unwrap() * e)))
        + T::from(9.0).unwrap() * q * (-T::from(22.0).unwrap() + e * (-T::from(330.0).unwrap() + e * (-T::from(924.0).unwrap() + e * (-T::from(660.0).unwrap() + e * (-T::from(110.0).unwrap() - T::from(2.0).unwrap() * e))))
        + T::from(11.0).unwrap() * q * (T::from(26.0).unwrap() + e * (T::from(572.0).unwrap() + e * (T::from(2574.0).unwrap() + e * (T::from(3432.0).unwrap() + e * (T::from(1430.0).unwrap() + e * (T::from(156.0).unwrap() + T::from(2.0).unwrap() * e)))))
        + T::from(13.0).unwrap() * q * (-T::from(30.0).unwrap() + e * (-T::from(910.0).unwrap() + e * (-T::from(6006.0).unwrap() + e * (-T::from(12870.0).unwrap() + e * (-T::from(10010.0).unwrap() + e * (-T::from(2730.0).unwrap() + e * (-T::from(210.0).unwrap() - T::from(2.0).unwrap() * e)))))
        + T::from(15.0).unwrap() * q * (T::from(34.0).unwrap() + e * (T::from(1360.0).unwrap() + e * (T::from(12376.0).unwrap() + e * (T::from(38896.0).unwrap() + e * (T::from(48620.0).unwrap() + e * (T::from(24752.0).unwrap() + e * (T::from(4760.0).unwrap() + e * (T::from(272.0).unwrap() + T::from(2.0).unwrap() * e)))))))
        + T::from(17.0).unwrap() * q * (-T::from(38.0).unwrap() + e * (-T::from(1938.0).unwrap() + e * (-T::from(23256.0).unwrap() + e * (-T::from(100776.0).unwrap() + e * (-T::from(184756.0).unwrap() + e * (-T::from(151164.0).unwrap() + e * (-T::from(54264.0).unwrap() + e * (-T::from(7752.0).unwrap() + e * (-T::from(342.0).unwrap() - T::from(2.0).unwrap() * e))))))))
        + T::from(19.0).unwrap() * q * (T::from(42.0).unwrap() + e * (T::from(2660.0).unwrap() + e * (T::from(40698.0).unwrap() + e * (T::from(232560.0).unwrap() + e * (T::from(587860.0).unwrap() + e * (T::from(705432.0).unwrap() + e * (T::from(406980.0).unwrap() + e * (T::from(108528.0).unwrap() + e * (T::from(11970.0).unwrap() + e * (T::from(420.0).unwrap() + T::from(2.0).unwrap() * e)))))))))
        + T::from(21.0).unwrap() * q * (-T::from(46.0).unwrap() + e * (-T::from(3542.0).unwrap() + e * (-T::from(67298.0).unwrap() + e * (-T::from(490314.0).unwrap() + e * (-T::from(1634380.0).unwrap() + e * (-T::from(2704156.0).unwrap() + e * (-T::from(2288132.0).unwrap() + e * (-T::from(980628.0).unwrap() + e * (-T::from(201894.0).unwrap() + e * (-T::from(17710.0).unwrap() + e * (-T::from(506.0).unwrap() - T::from(2.0).unwrap() * e))))))))))
        + T::from(23.0).unwrap() * q * (T::from(50.0).unwrap() + e * (T::from(4600.0).unwrap() + e * (T::from(106260.0).unwrap() + e * (T::from(961400.0).unwrap() + e * (T::from(4085950.0).unwrap() + e * (T::from(8914800.0).unwrap() + e * (T::from(10400600.0).unwrap() + e * (T::from(6537520.0).unwrap() + e * (T::from(2163150.0).unwrap() + e * (T::from(354200.0).unwrap() + e * (T::from(25300.0).unwrap() + e * (T::from(600.0).unwrap() + T::from(2.0).unwrap() * e)))))))))))
        + T::from(25.0).unwrap() * q * (-T::from(54.0).unwrap() + e * (-T::from(5850.0).unwrap() + e * (-T::from(161460.0).unwrap() + e * (-T::from(1776060.0).unwrap() + e * (-T::from(9373650.0).unwrap() + e * (-T::from(26075790.0).unwrap() + e * (-T::from(40116600.0).unwrap() + e * (-T::from(34767720.0).unwrap() + e * (-T::from(16872570.0).unwrap() + e * (-T::from(4440150.0).unwrap() + e * (-T::from(592020.0).unwrap() + e * (-T::from(35100.0).unwrap() + e * (-T::from(702.0).unwrap() - T::from(2.0).unwrap() * e))))))))))))
        + T::from(27.0).unwrap() * q * (T::from(58.0).unwrap() + e * (T::from(7308.0).unwrap() + e * (T::from(237510.0).unwrap() + e * (T::from(3121560.0).unwrap() + e * (T::from(20030010.0).unwrap() + e * (T::from(69194580.0).unwrap() + e * (T::from(135727830.0).unwrap() + e * (T::from(155117520.0).unwrap() + e * (T::from(103791870.0).unwrap() + e * (T::from(40060020.0).unwrap() + e * (T::from(8584290.0).unwrap() + e * (T::from(950040.0).unwrap() + e * (T::from(47502.0).unwrap() + e * (T::from(812.0).unwrap() + T::from(2.0).unwrap() * e)))))))))))))
        + T::from(29.0).unwrap() * q * (-T::from(62.0).unwrap() + e * (-T::from(8990.0).unwrap() + e * (-T::from(339822.0).unwrap() + e * (-T::from(5259150.0).unwrap() + e * (-T::from(40320150.0).unwrap() + e * (-T::from(169344630.0).unwrap() + e * (-T::from(412506150.0).unwrap() + e * (-T::from(601080390.0).unwrap() + e * (-T::from(530365050.0).unwrap() + e * (-T::from(282241050.0).unwrap() + e * (-T::from(88704330.0).unwrap() + e * (-T::from(15777450.0).unwrap() + e * (-T::from(14725620.0).unwrap() + e * (-T::from(629300.0).unwrap() + e * (-T::from(930.0).unwrap() - T::from(2.0).unwrap() * e)))))))))))))))
        + T::from(31.0).unwrap() * q * (T::from(66.0).unwrap() + e * (T::from(109120.0).unwrap() + e * (T::from(474672.0).unwrap() + e * (T::from(8544096.0).unwrap() + e * (T::from(77134200.0).unwrap() + e * (T::from(387073440.0).unwrap() + e * (T::from(1146332880.0).unwrap() + e * (T::from(2074316640.0).unwrap() + e * (T::from(2333606220.0).unwrap() + e * (T::from(1637618400.0).unwrap() + e * (T::from(709634640.0).unwrap() + e * (T::from(185122080.0).unwrap() + e * (T::from(27768312.0).unwrap() + e * (T::from(2215136.0).unwrap() + e * (T::from(81840.0).unwrap() + e * (T::from(1056.0).unwrap() + T::from(2.0).unwrap() * e)))))))))))))))
        + T::from(33.0).unwrap() * (-T::from(70.0).unwrap() + e * (-T::from(130900.0).unwrap() + e * (-T::from(649264.0).unwrap() + e * (-T::from(13449040.0).unwrap() + e * (-T::from(141214920.0).unwrap() + e * (-T::from(834451800.0).unwrap() + e * (-T::from(2952675600.0).unwrap() + e * (-T::from(6495886320.0).unwrap() + e * (-T::from(9075135300.0).unwrap() + e * (-T::from(8119857900.0).unwrap() + e * (-T::from(4639918800.0).unwrap() + e * (-T::from(1668903600.0).unwrap() + e * (-T::from(367158792.0).unwrap() + e * (-T::from(47071640.0).unwrap() + e * (-T::from(3246320.0).unwrap() + e * (-T::from(104720.0).unwrap() + e * (-T::from(1190.0).unwrap() - T::from(2.0).unwrap() * e))))))))))))))))) * q))))))))))))))));

    let b = T::from(ONE_OVER_SQRT_TWO_PI).unwrap() * ((T::from(-0.5).unwrap() * (h * h + t * t)).exp()) * (t / r) * asymptotic_expansion_sum;
    b.abs().max(T::zero())
}
