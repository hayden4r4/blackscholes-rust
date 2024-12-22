use assert_approx_eq::assert_approx_eq;
use blackscholes::OptionType;
use blackscholes::valuators::black_scholes::{Inputs, Pricing};

const INPUTS_CALL_OTM: Inputs = Inputs {
    option_type: OptionType::Call,
    s: 100.0,
    k: 110.0,
    p: None,
    r: 0.05,
    q: 0.05,
    t: 20.0 / 365.25,
    sigma: Some(0.2),
};
const INPUTS_CALL_ITM: Inputs = Inputs {
    option_type: OptionType::Call,
    s: 100.0,
    k: 90.0,
    p: None,
    r: 0.05,
    q: 0.05,
    t: 20.0 / 365.25,
    sigma: Some(0.2),
};
const INPUTS_PUT_OTM: Inputs = Inputs {
    option_type: OptionType::Put,
    s: 100.0,
    k: 90.0,
    p: None,
    r: 0.05,
    q: 0.05,
    t: 20.0 / 365.25,
    sigma: Some(0.2),
};
const INPUTS_PUT_ITM: Inputs = Inputs {
    option_type: OptionType::Put,
    s: 100.0,
    k: 110.0,
    p: None,
    r: 0.05,
    q: 0.05,
    t: 20.0 / 365.25,
    sigma: Some(0.2),
};

const INPUTS_BRANCH_CUT: Inputs = Inputs {
    option_type: OptionType::Put,
    s: 100.0,
    k: 100.0,
    p: None,
    r: 0.0,
    q: 0.0,
    sigma: Some(0.421),
    t: 1.0,
};

#[test]
fn price_call_otm() {
    assert_approx_eq!(INPUTS_CALL_OTM.calc_price().unwrap(), 0.0376, 0.001);
}
#[test]
fn price_call_itm() {
    assert_approx_eq!(INPUTS_CALL_ITM.calc_price().unwrap(), 9.9913, 0.001);
}

#[test]
fn price_put_otm() {
    assert_approx_eq!(INPUTS_PUT_OTM.calc_price().unwrap(), 0.01867, 0.001);
}
#[test]
fn price_put_itm() {
    assert_approx_eq!(INPUTS_PUT_ITM.calc_price().unwrap(), 10.0103, 0.001);
}

#[test]
fn price_using_lets_be_rational() {
    // compare the results from calc_price() and calc_rational_price() for the options above
    assert_approx_eq!(
        INPUTS_CALL_OTM.calc_price().unwrap(),
        INPUTS_CALL_OTM.calc_rational_price().unwrap(),
        0.001
    );

    assert_approx_eq!(
        INPUTS_CALL_ITM.calc_price().unwrap(),
        INPUTS_CALL_ITM.calc_rational_price().unwrap(),
        0.001
    );

    assert_approx_eq!(
        INPUTS_PUT_OTM.calc_price().unwrap(),
        INPUTS_PUT_OTM.calc_rational_price().unwrap(),
        0.001
    );

    assert_approx_eq!(
        INPUTS_PUT_ITM.calc_price().unwrap(),
        INPUTS_PUT_ITM.calc_rational_price().unwrap(),
        0.001
    );
}

#[test]
fn test_rational_price_near_branch_cut() {
    assert_approx_eq!(
        INPUTS_BRANCH_CUT.calc_rational_price().unwrap(),
        16.67224,
        0.001
    );
}
