//! Local constants
use core::f64;

pub(crate) use statrs::consts::SQRT_2PI;
pub(crate) const PI: f64 = std::f64::consts::PI;
pub(crate) const FRAC_1_SQRT_2: f64 = f64::consts::FRAC_1_SQRT_2;
pub(crate) const H_LARGE: f64 = -10.0;
pub(crate) const ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD: f64 = -10.0;
pub(crate) const SIXTEENTH_ROOT_DBL_EPSILON: f64 = 0.10566243270259357;
pub(crate) const CODYS_THRESHOLD: f64 = 0.46875;
pub(crate) const SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD: f64 =
    2.0 * SIXTEENTH_ROOT_DBL_EPSILON;
pub(crate) const IMPLIED_VOLATILITY_MAXIMUM_ITERATIONS: i32 = 2;
pub(crate) const DENORMALISATION_CUTOFF: f64 = 0.0;
pub(crate) const ONE_OVER_SQRT_TWO_PI: f64 = 1.0 / SQRT_2PI;
pub(crate) const VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM: f64 = f64::MAX;
pub(crate) const VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC: f64 = -f64::MAX;
pub(crate) const SQRT_DBL_MIN: f64 = 1.4916681462400413e-154;
pub(crate) const SQRT_DBL_MAX: f64 = 1.340_780_792_994_259_6e154;
pub(crate) const SQRT_THREE: f64 = 1.732_050_807_568_877_2; //3.0_f64.sqrt();
pub(crate) const TWO_PI_OVER_SQRT_TWENTY_SEVEN: f64 = 1.209_199_576_156_145_2; // f64::PI() * 2.0 / sqrt(27.0);
pub(crate) const SQRT_PI_OVER_TWO: f64 = 1.253_314_137_315_500_3; //sqrt(f64::PI() / 2.0_f64);
pub(crate) const SQRT_ONE_OVER_THREE: f64 = 0.577_350_269_189_625_7; // sqrt(1.0 / 3.0_f64);
pub(crate) const PI_OVER_SIX: f64 = PI / 6.0;
pub(crate) const FOURTH_ROOT_DBL_EPSILON: f64 = 0.0001220703125;
pub(crate) const NORMALISED_X2_THRESHOLD: f64 = 98.0 * FOURTH_ROOT_DBL_EPSILON;
pub(crate) const MINIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE: f64 =
    -(1.0 - 1.4901161193847656e-8); // -(1.0 - f64::EPSILON.sqrt());
pub(crate) const MAXIMUM_RATIONAL_CUBIC_CONTROL_PARAMETER_VALUE: f64 =
    2.0 / (f64::EPSILON * f64::EPSILON);
pub(crate) const FRAC_1_SQRT_PI: f64 = 0.564_189_583_547_756_3;
pub(crate) const THRESHOLD: f64 = 0.46875;
pub(crate) const XNEG: f64 = -26.6287357137514;
pub(crate) const XBIG: f64 = 26.543;
