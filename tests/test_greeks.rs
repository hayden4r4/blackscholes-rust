#[cfg(test)]
mod tests {
    use blackscholes::{Greeks, Inputs, OptionType};

    #[test]
    fn test_calc_delta_zero_stock_price() {
        // arrange
        let option_type = OptionType::Call;
        let s = 0.0; // extreme value
        let k = 100.0;
        let p = None;
        let r = 0.05;
        let q = 0.0;
        let sigma = Some(0.2);
        let t = 20.0 / 365.25;

        let inputs = Inputs::new(option_type, s, k, p, r, q, t, sigma);

        // act
        let result = inputs.calc_delta();

        // assert
        assert!(result.is_err());
    }

    #[test]
    fn test_calc_delta_zero_strike_price() {
        // arrange
        let option_type = OptionType::Call;
        let s = 100.0;
        let k = 0.0;
        let p = None;
        let r = 0.05;
        let q = 0.0;
        let sigma = Some(0.2);
        let t = 20.0 / 365.25;

        let inputs = Inputs::new(option_type, s, k, p, r, q, t, sigma);

        // act
        let result = inputs.calc_delta();

        // assert
        assert!(result.is_err());
    }

    #[test]
    fn test_calc_delta_zero_risk_free_rate() {
        // arrange
        let option_type = OptionType::Call;
        let s = 100.0;
        let k = 100.0;
        let p = None;
        let r = 0.0; // extreme value
        let q = 0.0;
        let sigma = Some(0.2);
        let t = 20.0 / 365.25;

        let inputs = Inputs::new(option_type, s, k, p, r, q, t, sigma);

        // act
        let result = inputs.calc_delta();

        // assert
        assert!(result.is_ok());
    }

    #[test]
    fn test_calc_delta_none_volatility() {
        // arrange
        let option_type = OptionType::Call;
        let s = 100.0;
        let k = 100.0;
        let p = None;
        let r = 0.05;
        let q = 0.0;
        let sigma = None; // extreme value
        let t = 20.0 / 365.25;

        let inputs = Inputs::new(option_type, s, k, p, r, q, t, sigma);

        // act
        let result = inputs.calc_delta();

        // assert
        assert!(result.is_err());
    }

    #[test]
    fn test_calc_delta_zero_time_to_maturity() {
        // arrange
        let option_type = OptionType::Call;
        let s = 100.0;
        let k = 100.0;
        let p = None;
        let r = 0.05;
        let q = 0.0;
        let sigma = Some(0.2);
        let t = 0.0; // extreme value

        let inputs = Inputs::new(option_type, s, k, p, r, q, t, sigma);

        // act
        let result = inputs.calc_delta();

        // assert
        assert!(result.is_err());
    }
}
