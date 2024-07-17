use num_traits::{AsPrimitive, Float, FloatConst, FromPrimitive};
use std::collections::HashMap;

use crate::{
    calc_d1d2, calc_nd1nd2, calc_nprimed1, calc_nprimed2, Inputs, OptionType, Pricing,
    DAYS_PER_YEAR,
};

pub trait Greeks<T>: Pricing<T>
where
    T: Float,
{
    fn calc_delta(&self) -> Result<T, String>;
    fn calc_gamma(&self) -> Result<T, String>;
    fn calc_theta(&self) -> Result<T, String>;
    fn calc_vega(&self) -> Result<T, String>;
    fn calc_rho(&self) -> Result<T, String>;
    fn calc_epsilon(&self) -> Result<T, String>;
    fn calc_lambda(&self) -> Result<T, String>;
    fn calc_vanna(&self) -> Result<T, String>;
    fn calc_charm(&self) -> Result<T, String>;
    fn calc_veta(&self) -> Result<T, String>;
    fn calc_vomma(&self) -> Result<T, String>;
    fn calc_speed(&self) -> Result<T, String>;
    fn calc_zomma(&self) -> Result<T, String>;
    fn calc_color(&self) -> Result<T, String>;
    fn calc_ultima(&self) -> Result<T, String>;
    fn calc_dual_delta(&self) -> Result<T, String>;
    fn calc_dual_gamma(&self) -> Result<T, String>;
    fn calc_all_greeks(&self) -> Result<HashMap<String, T>, String>;
}

impl<T> Greeks<T> for Inputs<T>
where
    T: Float + FromPrimitive + AsPrimitive<f64> + FloatConst,
{
    /// Calculates the delta of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the delta of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let delta = inputs.calc_delta().unwrap();
    /// ```
    fn calc_delta(&self) -> Result<T, String> {
        let (nd1, _): (T, T) = calc_nd1nd2(self)?;

        let option_type: T = T::from::<f64>(self.option_type.into()).unwrap();
        let delta = option_type * T::E().powf(-self.q * self.t) * nd1;

        Ok(delta)
    }

    /// Calculates the gamma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the gamma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let gamma = inputs.calc_gamma().unwrap();
    /// ```
    fn calc_gamma(&self) -> Result<T, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(T) for self.sigma, received None")?;

        let nprimed1: T = calc_nprimed1(self)?;
        let gamma: T = T::E().powf(-self.q * self.t) * nprimed1 / (self.s * sigma * self.t.sqrt());
        Ok(gamma)
    }

    /// Calculates the theta of the option.
    /// Uses 365.25 days in a year for calculations.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of theta per day (not per year).
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let theta = inputs.calc_theta().unwrap();
    /// ```
    fn calc_theta(&self) -> Result<T, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(T) for self.sigma, received None")?;

        let nprimed1: T = calc_nprimed1(self)?;
        let (nd1, nd2): (T, T) = calc_nd1nd2(self)?;

        // Calculation uses 365.25 for T: Time of days per year.
        let option_type: T = T::from::<f64>(self.option_type.into()).unwrap();
        let theta = (-(self.s * sigma * T::E().powf(-self.q * self.t) * nprimed1
            / (T::from(2.).unwrap() * self.t.sqrt()))
            - self.r * self.k * T::E().powf(-self.r * self.t) * nd2 * option_type
            + self.q * self.s * T::E().powf(-self.q * self.t) * nd1 * option_type)
            / T::from(DAYS_PER_YEAR).unwrap();

        Ok(theta)
    }

    /// Calculates the vega of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the vega of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let vega = inputs.calc_vega().unwrap();
    /// ```
    fn calc_vega(&self) -> Result<T, String> {
        let nprimed1: T = calc_nprimed1(self)?;
        let vega: T = T::from(0.01).unwrap()
            * self.s
            * T::E().powf(-self.q * self.t)
            * self.t.sqrt()
            * nprimed1;
        Ok(vega)
    }

    /// Calculates the rho of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the rho of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let rho = inputs.calc_rho().unwrap();
    /// ```
    fn calc_rho(&self) -> Result<T, String> {
        let (_, nd2): (T, T) = calc_nd1nd2(self)?;

        let option_type: T = T::from::<f64>(self.option_type.into()).unwrap();
        let rho = option_type
            * T::from(0.01).unwrap()
            * self.k
            * self.t
            * T::E().powf(-self.r * self.t)
            * nd2;

        Ok(rho)
    }

    // The formulas for the greeks below are from the wikipedia page for the Black-Scholes greeks
    // https://en.wikipedia.org/wiki/Greeks_(finance)#Black.E2.80.93Scholes_Greeks
    // Some sources I reviewed contain variations of these formulas and/or varying values, therefore the
    // values returned by this library may not match other libraries or sources.
    // These functions have not been throughouly tested.

    /// Calculates the epsilon of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the epsilon of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let epsilon = inputs.calc_epsilon().unwrap();
    /// ```
    fn calc_epsilon(&self) -> Result<T, String> {
        let (nd1, _) = calc_nd1nd2(self)?;
        let e_negqt = T::E().powf(-self.q * self.t);

        let option_type: T = T::from::<f64>(self.option_type.into()).unwrap();
        let epsilon: T = -self.s * self.t * e_negqt * nd1 * option_type;

        Ok(epsilon)
    }

    /// Calculates the lambda of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the lambda of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks, Pricing};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let lambda = inputs.calc_lambda().unwrap();
    /// ```
    fn calc_lambda(&self) -> Result<T, String> {
        let delta = self.calc_delta()?;
        Ok(delta * self.s / self.calc_price()?)
    }

    /// Calculates the vanna of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the vanna of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let vanna = inputs.calc_vanna().unwrap();
    /// ```
    fn calc_vanna(&self) -> Result<T, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(T) for self.sigma, received None")?;

        let nprimed1 = calc_nprimed1(self)?;
        let (_, d2) = calc_d1d2(self)?;
        let vanna: T =
            d2 * T::E().powf(-self.q * self.t) * nprimed1 * T::from(-0.01).unwrap() / sigma;
        Ok(vanna)
    }

    // /// Calculates the charm of the option.
    // /// # Requires
    // /// s, k, r, q, t, sigma
    // /// # Returns
    // /// T of the charm of the option.
    // /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let charm = inputs.calc_charm().unwrap();
    /// ```
    fn calc_charm(&self) -> Result<T, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(T) for self.sigma, received None")?;
        let nprimed1 = calc_nprimed1(self)?;
        let (nd1, _) = calc_nd1nd2(self)?;
        let (_, d2) = calc_d1d2(self)?;
        let e_negqt = T::E().powf(-self.q * self.t);

        let option_type: T = T::from::<f64>(self.option_type.into()).unwrap();
        let charm: T = option_type * self.q * e_negqt * nd1
            - e_negqt
                * nprimed1
                * (T::from(2.0).unwrap() * (self.r - self.q) * self.t - d2 * sigma * self.t.sqrt())
                / (T::from(2.0).unwrap() * self.t * sigma * self.t.sqrt());

        Ok(charm)
    }

    /// Calculates the veta of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the veta of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let veta = inputs.calc_veta().unwrap();
    /// ```
    fn calc_veta(&self) -> Result<T, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(T) for self.sigma, received None")?;
        let nprimed1 = calc_nprimed1(self)?;
        let (d1, d2) = calc_d1d2(self)?;
        let e_negqt = T::E().powf(-self.q * self.t);

        let veta = -self.s
            * e_negqt
            * nprimed1
            * self.t.sqrt()
            * (self.q + ((self.r - self.q) * d1) / (sigma * self.t.sqrt())
                - ((T::one() + d1 * d2) / (T::from(2.0).unwrap() * self.t)));
        Ok(veta)
    }

    /// Calculates the vomma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the vomma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let vomma = inputs.calc_vomma().unwrap();
    /// ```
    fn calc_vomma(&self) -> Result<T, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(T) for self.sigma, received None")?;
        let (d1, d2) = calc_d1d2(self)?;

        let vomma = Inputs::calc_vega(self)? * ((d1 * d2) / sigma);
        Ok(vomma)
    }

    /// Calculates the speed of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the speed of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let speed = inputs.calc_speed().unwrap();
    /// ```
    fn calc_speed(&self) -> Result<T, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(T) for self.sigma, received None")?;
        let (d1, _) = calc_d1d2(self)?;
        let gamma = Inputs::calc_gamma(self)?;

        let speed = -gamma / self.s * (d1 / (sigma * self.t.sqrt()) + T::one());
        Ok(speed)
    }

    /// Calculates the zomma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the zomma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let zomma = inputs.calc_zomma().unwrap();
    /// ```
    fn calc_zomma(&self) -> Result<T, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(T) for self.sigma, received None")?;
        let (d1, d2) = calc_d1d2(self)?;
        let gamma = Inputs::calc_gamma(self)?;

        let zomma = gamma * ((d1 * d2 - T::one()) / sigma);
        Ok(zomma)
    }

    /// Calculates the color of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the color of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let color = inputs.calc_color().unwrap();
    /// ```
    fn calc_color(&self) -> Result<T, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(T) for self.sigma, received None")?;
        let (d1, d2) = calc_d1d2(self)?;
        let nprimed1 = calc_nprimed1(self)?;
        let e_negqt = T::E().powf(-self.q * self.t);

        let color = -e_negqt
            * (nprimed1 / (T::from(2.0).unwrap() * self.s * self.t * sigma * self.t.sqrt()))
            * (T::from(2.0).unwrap() * self.q * self.t
                + T::one()
                + (T::from(2.0).unwrap() * (self.r - self.q) * self.t
                    - d2 * sigma * self.t.sqrt())
                    / (sigma * self.t.sqrt())
                    * d1);
        Ok(color)
    }

    /// Calculates the ultima of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the ultima of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let ultima = inputs.calc_ultima().unwrap();
    /// ```
    fn calc_ultima(&self) -> Result<T, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(T) for self.sigma, received None")?;
        let (d1, d2) = calc_d1d2(self)?;
        let vega = Inputs::calc_vega(self)?;

        let ultima = -vega / sigma.powf(T::from(2.0).unwrap())
            * (d1 * d2 * (T::one() - d1 * d2)
                + d1.powf(T::from(2.0).unwrap())
                + d2.powf(T::from(2.0).unwrap()));
        Ok(ultima)
    }

    /// Calculates the dual delta of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the dual delta of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let dual_delta = inputs.calc_dual_delta().unwrap();
    /// ```
    fn calc_dual_delta(&self) -> Result<T, String> {
        let (_, nd2) = calc_nd1nd2(self)?;
        let e_negqt = T::E().powf(-self.q * self.t);

        let dual_delta = match self.option_type {
            OptionType::Call => -e_negqt * nd2,
            OptionType::Put => e_negqt * nd2,
        };
        Ok(dual_delta)
    }

    /// Calculates the dual gamma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// T of the dual gamma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let dual_gamma = inputs.calc_dual_gamma().unwrap();
    /// ```
    fn calc_dual_gamma(&self) -> Result<T, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(T) for self.sigma, received None")?;
        let nprimed2 = calc_nprimed2(self)?;
        let e_negqt = T::E().powf(-self.q * self.t);

        let dual_gamma = e_negqt * (nprimed2 / (self.k * sigma * self.t.sqrt()));
        Ok(dual_gamma)
    }

    /// Calculates all Greeks of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// HashMap of type <String, T> of all Greeks of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let greeks = inputs.calc_all_greeks().unwrap();
    /// ```
    fn calc_all_greeks(&self) -> Result<HashMap<String, T>, String> {
        let mut greeks: HashMap<String, T> = HashMap::with_capacity(17);
        greeks.insert("delta".into(), self.calc_delta()?);
        greeks.insert("gamma".into(), self.calc_gamma()?);
        greeks.insert("theta".into(), self.calc_theta()?);
        greeks.insert("vega".into(), self.calc_vega()?);
        greeks.insert("rho".into(), self.calc_rho()?);
        greeks.insert("epsilon".into(), self.calc_epsilon()?);
        greeks.insert("lambda".into(), self.calc_lambda()?);
        greeks.insert("vanna".into(), self.calc_vanna()?);
        greeks.insert("charm".into(), self.calc_charm()?);
        greeks.insert("veta".into(), self.calc_veta()?);
        greeks.insert("vomma".into(), self.calc_vomma()?);
        greeks.insert("speed".into(), self.calc_speed()?);
        greeks.insert("zomma".into(), self.calc_zomma()?);
        greeks.insert("color".into(), self.calc_color()?);
        greeks.insert("ultima".into(), self.calc_ultima()?);
        greeks.insert("dual_delta".into(), self.calc_dual_delta()?);
        greeks.insert("dual_gamma".into(), self.calc_dual_gamma()?);
        Ok(greeks)
    }
}
