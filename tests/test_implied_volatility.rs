mod tests {
    use assert_approx_eq::assert_approx_eq;
    use blackscholes::{ImpliedVolatility, Inputs, OptionType, Pricing};

    // Tolerance is a bit higher due to IV being an approximation
    const TOLERANCE: f64 = 1e-8;

    #[test]
    fn test_put_otm_rational_iv() {
        // arrange
        let sigma: Option<f64> = Some(0.25);
        let mut inputs_put_otm: Inputs = Inputs {
            option_type: OptionType::Put,
            s: 90.0,
            k: 100.0,
            p: None,
            r: 0.03,
            q: 0.02,
            t: 45.0 / 365.25,
            sigma,
        };

        let price = inputs_put_otm.calc_price().unwrap();

        inputs_put_otm.p = Some(price);
        inputs_put_otm.sigma = None;

        // act
        let iv = inputs_put_otm.calc_rational_iv().unwrap();

        // assert
        println!("Put OTM: {}", iv);
        assert_approx_eq!(iv, sigma.unwrap(), TOLERANCE);
    }

    #[test]
    fn test_call_itm_rational_iv() {
        // arrange
        let sigma: Option<f64> = Some(0.15);
        let mut inputs_call_itm: Inputs = Inputs {
            option_type: OptionType::Call,
            s: 120.0,
            k: 100.0,
            p: None,
            r: 0.01,
            q: 0.0,
            t: 60.0 / 365.25,
            sigma,
        };

        let price = inputs_call_itm.calc_price().unwrap();

        inputs_call_itm.p = Some(price);
        inputs_call_itm.sigma = None;

        // act
        let iv = inputs_call_itm.calc_rational_iv().unwrap();

        // assert
        println!("Call ITM: {}", iv);
        assert_approx_eq!(iv, sigma.unwrap(), TOLERANCE);
    }

    #[test]
    fn test_put_itm_rational_iv() {
        // arrange
        let sigma: Option<f64> = Some(0.18);
        let mut inputs_put_itm: Inputs = Inputs {
            option_type: OptionType::Put,
            s: 80.0,
            k: 100.0,
            p: None,
            r: 0.04,
            q: 0.03,
            t: 60.0 / 365.25,
            sigma,
        };

        let price = inputs_put_itm.calc_price().unwrap();

        inputs_put_itm.p = Some(price);
        inputs_put_itm.sigma = None;

        // act
        let iv = inputs_put_itm.calc_rational_iv().unwrap();

        // assert
        println!("Put ITM: {}", iv);
        assert_approx_eq!(iv, sigma.unwrap(), TOLERANCE);
    }

    #[test]
    fn test_call_atm_rational_iv() {
        // arrange
        let sigma: Option<f64> = Some(0.2);
        let mut inputs_call_atm: Inputs = Inputs {
            option_type: OptionType::Call,
            s: 100.0,
            k: 100.0,
            p: None,
            r: 0.05,
            q: 0.04,
            t: 90.0 / 365.25,
            sigma,
        };

        let price = inputs_call_atm.calc_price().unwrap();

        inputs_call_atm.p = Some(price);
        inputs_call_atm.sigma = None;

        // act
        let iv = inputs_call_atm.calc_rational_iv().unwrap();

        // assert
        println!("Call ATM: {}", iv);
        assert_approx_eq!(iv, sigma.unwrap(), TOLERANCE);
    }

    #[test]
    fn test_put_atm_rational_iv() {
        // arrange
        let sigma: Option<f64> = Some(0.22);
        let mut inputs_put_atm: Inputs = Inputs {
            option_type: OptionType::Put,
            s: 100.0,
            k: 100.0,
            p: None,
            r: 0.06,
            q: 0.01,
            t: 120.0 / 365.25,
            sigma,
        };

        let price = inputs_put_atm.calc_price().unwrap();

        inputs_put_atm.p = Some(price);
        inputs_put_atm.sigma = None;

        // act
        let iv = inputs_put_atm.calc_rational_iv().unwrap();

        // assert
        println!("Put ATM: {}", iv);
        assert_approx_eq!(iv, sigma.unwrap(), TOLERANCE);
    }
}
