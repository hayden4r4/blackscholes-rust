use assert_approx_eq::assert_approx_eq;

use blackscholes::{Inputs, OptionType, Pricing};

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

#[test]
fn price_call_otm() {
    assert!((INPUTS_CALL_OTM.calc_price().unwrap() - 0.0376).abs() < 0.001);
}
#[test]
fn price_call_itm() {
    assert!((INPUTS_CALL_ITM.calc_price().unwrap() - 9.9913).abs() < 0.001);
}

#[test]
fn price_put_otm() {
    assert!((INPUTS_PUT_OTM.calc_price().unwrap() - 0.01867).abs() < 0.001);
}
#[test]
fn price_put_itm() {
    assert!((INPUTS_PUT_ITM.calc_price().unwrap() - 10.0103).abs() < 0.001);
}
#[test]
fn price_using_lets_be_rational() {
    // compare the results from calc_price() and calc_rational_price() for the options above
    assert_approx_eq!(
        INPUTS_CALL_OTM.calc_price().unwrap() as f64,
        INPUTS_CALL_OTM.calc_rational_price().unwrap(),
        0.001
    );

    assert_approx_eq!(
        INPUTS_CALL_ITM.calc_price().unwrap() as f64,
        INPUTS_CALL_ITM.calc_rational_price().unwrap(),
        0.001
    );

    assert_approx_eq!(
        INPUTS_PUT_OTM.calc_price().unwrap() as f64,
        INPUTS_PUT_OTM.calc_rational_price().unwrap(),
        0.001
    );

    assert_approx_eq!(
        INPUTS_PUT_ITM.calc_price().unwrap() as f64,
        INPUTS_PUT_ITM.calc_rational_price().unwrap(),
        0.001
    );
}
