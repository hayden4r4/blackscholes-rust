use assert_approx_eq::assert_approx_eq;
use blackscholes::{
    valuators::{black_76::Inputs as b76Inputs, black_scholes::Inputs as bsInputs},
    OptionType, Pricing,
};

const TOLERANCE_LARGE: f64 = 1e-3;
const TOLERANCE_SMALL: f64 = 1e-8;

const INPUTS_CALL_OTM: bsInputs = bsInputs {
    option_type: OptionType::Call,
    s: 100.0,
    k: 110.0,
    p: None,
    r: 0.05,
    q: 0.05,
    t: 20.0 / 365.25,
    sigma: Some(0.2),
};
const INPUTS_CALL_ITM: bsInputs = bsInputs {
    option_type: OptionType::Call,
    s: 100.0,
    k: 90.0,
    p: None,
    r: 0.05,
    q: 0.05,
    t: 20.0 / 365.25,
    sigma: Some(0.2),
};
const INPUTS_PUT_OTM: bsInputs = bsInputs {
    option_type: OptionType::Put,
    s: 100.0,
    k: 90.0,
    p: None,
    r: 0.05,
    q: 0.05,
    t: 20.0 / 365.25,
    sigma: Some(0.2),
};
const INPUTS_PUT_ITM: bsInputs = bsInputs {
    option_type: OptionType::Put,
    s: 100.0,
    k: 110.0,
    p: None,
    r: 0.05,
    q: 0.05,
    t: 20.0 / 365.25,
    sigma: Some(0.2),
};

// Black-Scholes
#[test]
fn price_call_otm_bs() {
    assert_approx_eq!(INPUTS_CALL_OTM.calc_price().unwrap(), 0.0376, TOLERANCE_LARGE);
}
#[test]
fn price_call_itm_bs() {
    assert_approx_eq!(INPUTS_CALL_ITM.calc_price().unwrap(), 9.9913, TOLERANCE_LARGE);
}

#[test]
fn price_put_otm_bs() {
    assert_approx_eq!(INPUTS_PUT_OTM.calc_price().unwrap(), 0.01867, TOLERANCE_LARGE);
}
#[test]
fn price_put_itm_bs() {
    assert_approx_eq!(INPUTS_PUT_ITM.calc_price().unwrap(), 10.0103, TOLERANCE_LARGE);
}