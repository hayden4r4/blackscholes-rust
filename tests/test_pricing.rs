use assert_approx_eq::assert_approx_eq;
use option_valuators::{
    valuators::{black_76::Inputs as b76Inputs, black_scholes::Inputs as bsInputs},
    OptionType, Pricing,
};

const TOLERANCE_LARGE: f64 = 1e-3;
const TOLERANCE_SMALL: f64 = 1e-5;


// B76 tests rely on r = q for these inputs
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

fn bs_inputs_to_b76_inputs(inputs: bsInputs) -> b76Inputs {
    b76Inputs {
        option_type: inputs.option_type,
        f: inputs.s,
        k: inputs.k,
        p: inputs.p,
        r: inputs.r,
        t: inputs.t,
        sigma: inputs.sigma,
        shifted: true,
    }
}

// Black-76
#[test]
fn price_call_otm_b76() {
    let b76_inputs: b76Inputs = bs_inputs_to_b76_inputs(INPUTS_CALL_OTM);
    assert_approx_eq!(b76_inputs.calc_price().unwrap(), INPUTS_CALL_OTM.calc_price().unwrap(), TOLERANCE_SMALL);
}

#[test]
fn price_call_itm_b76() {
    let b76_inputs: b76Inputs = bs_inputs_to_b76_inputs(INPUTS_CALL_ITM);
    assert_approx_eq!(b76_inputs.calc_price().unwrap(), INPUTS_CALL_ITM.calc_price().unwrap(), TOLERANCE_SMALL);
}

#[test]
fn price_put_otm_b76() {
    let b76_inputs: b76Inputs = bs_inputs_to_b76_inputs(INPUTS_PUT_OTM);
    assert_approx_eq!(b76_inputs.calc_price().unwrap(), INPUTS_PUT_OTM.calc_price().unwrap(), TOLERANCE_SMALL);
}

#[test]
fn price_put_itm_b76() {
    let b76_inputs: b76Inputs = bs_inputs_to_b76_inputs(INPUTS_PUT_ITM);
    assert_approx_eq!(b76_inputs.calc_price().unwrap(), INPUTS_PUT_ITM.calc_price().unwrap(), TOLERANCE_SMALL);
}