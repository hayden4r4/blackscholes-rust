use std::{
    collections::HashMap,
    fmt::{Display, Formatter},
};

use num_traits::Float;

use crate::{Inputs, OptionType, Pricing, *};

#[derive(Hash, Eq, PartialEq, Debug)]
pub enum Greek {
    Delta,
    Gamma,
    Theta,
    Vega,
    Rho,
    Epsilon,
    Lambda,
    Vanna,
    Charm,
    Veta,
    Vomma,
    Speed,
    Zomma,
    Color,
    Ultima,
    DualDelta,
    DualGamma,
}

impl Display for Greek {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

pub trait Greeks<T>: Pricing<T>
where
    T: Float,
{
    fn calc_delta(&self) -> Result<T, BlackScholesError>;
    fn calc_gamma(&self) -> Result<T, BlackScholesError>;
    fn calc_theta(&self) -> Result<T, BlackScholesError>;
    fn calc_vega(&self) -> Result<T, BlackScholesError>;
    fn calc_rho(&self) -> Result<T, BlackScholesError>;
    fn calc_epsilon(&self) -> Result<T, BlackScholesError>;
    fn calc_lambda(&self) -> Result<T, BlackScholesError>;
    fn calc_vanna(&self) -> Result<T, BlackScholesError>;
    fn calc_charm(&self) -> Result<T, BlackScholesError>;
    fn calc_veta(&self) -> Result<T, BlackScholesError>;
    fn calc_vomma(&self) -> Result<T, BlackScholesError>;
    fn calc_speed(&self) -> Result<T, BlackScholesError>;
    fn calc_zomma(&self) -> Result<T, BlackScholesError>;
    fn calc_color(&self) -> Result<T, BlackScholesError>;
    fn calc_ultima(&self) -> Result<T, BlackScholesError>;
    fn calc_dual_delta(&self) -> Result<T, BlackScholesError>;
    fn calc_dual_gamma(&self) -> Result<T, BlackScholesError>;
    fn calc_all_greeks(&self) -> Result<HashMap<Greek, T>, BlackScholesError>;
}

const PERCENTAGE_CONVERSION: f64 = 0.01;

impl Greeks<f64> for Inputs {
    /// Calculates the delta of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f64 of the delta of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let delta = inputs.calc_delta().unwrap();
    /// ```
    fn calc_delta(&self) -> Result<f64, BlackScholesError> {
        let (nd1, _): (f64, f64) = calc_nd1nd2(self)?;

        Ok(self.calc_delta_with_params(nd1, self.q, self.t, self.option_type))
    }

    /// Calculates the gamma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f64 of the gamma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let gamma = inputs.calc_gamma().unwrap();
    /// ```
    fn calc_gamma(&self) -> Result<f64, BlackScholesError> {
        let sigma = self.sigma.ok_or(BlackScholesError::ExpectedSigma)?;

        let nprimed1: f64 = calc_nprimed1(self)?;
        Ok(self.calc_gamma_with_params(nprimed1, self.s, sigma, self.q, self.t))
    }

    /// Calculates the theta of the option.
    /// Uses DAYS_PER_YEAR days in a year for calculations.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f64 of theta per day (not per year).
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let theta = inputs.calc_theta().unwrap();
    /// ```
    fn calc_theta(&self) -> Result<f64, BlackScholesError> {
        let sigma = self.sigma.ok_or(BlackScholesError::ExpectedSigma)?;

        let nprimed1: f64 = calc_nprimed1(self)?;
        let (nd1, nd2): (f64, f64) = calc_nd1nd2(self)?;

        Ok(self.calc_theta_with_params(
            sigma,
            nprimed1,
            nd1,
            nd2,
            self.s,
            self.q,
            self.t,
            self.r,
            self.k,
            self.option_type,
        ))
    }

    /// Calculates the vega of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f64 of the vega of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let vega = inputs.calc_vega().unwrap();
    /// ```
    fn calc_vega(&self) -> Result<f64, BlackScholesError> {
        let nprimed1: f64 = calc_nprimed1(self)?;

        Ok(self.calc_vega_with_params(nprimed1, self.s, self.q, self.t))
    }

    /// Calculates the rho of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f64 of the rho of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let rho = inputs.calc_rho().unwrap();
    /// ```
    fn calc_rho(&self) -> Result<f64, BlackScholesError> {
        let (_, nd2): (f64, f64) = calc_nd1nd2(self)?;

        Ok(self.calc_rho_with_params(self.k, self.t, self.r, nd2, self.option_type))
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
    /// f64 of the epsilon of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let epsilon = inputs.calc_epsilon().unwrap();
    /// ```
    fn calc_epsilon(&self) -> Result<f64, BlackScholesError> {
        let (nd1, _) = calc_nd1nd2(self)?;
        let e_negqt = self.calc_e_negqt();

        Ok(self.calc_epsilon_with_params(self.s, self.t, e_negqt, nd1, self.option_type))
    }

    /// Calculates the lambda of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f64 of the lambda of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks, Pricing};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let lambda = inputs.calc_lambda()?;
    /// ```
    fn calc_lambda(&self) -> Result<f64, BlackScholesError> {
        let (nd1, nd2) = calc_nd1nd2(self)?;
        let delta = self.calc_delta_with_params(nd1, self.q, self.t, self.option_type);
        let price = self.calc_price_with_params(
            nd1,
            nd2,
            self.s,
            self.q,
            self.t,
            self.k,
            self.r,
            self.option_type,
        );
        Ok(self.calc_lambda_with_params(delta, self.s, price))
    }

    /// Calculates the vanna of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f64 of the vanna of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let vanna = inputs.calc_vanna().unwrap();
    /// ```
    fn calc_vanna(&self) -> Result<f64, BlackScholesError> {
        let sigma = self.sigma.ok_or(BlackScholesError::ExpectedSigma)?;

        let nprimed1 = calc_nprimed1(self)?;
        let (_, d2) = calc_d1d2(self)?;
        let e_negqt = self.calc_e_negqt();

        Ok(self.calc_vanna_with_params(d2, e_negqt, nprimed1, sigma))
    }

    // /// Calculates the charm of the option.
    // /// # Requires
    // /// s, k, r, q, t, sigma
    // /// # Returns
    // /// f64 of the charm of the option.
    // /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let charm = inputs.calc_charm().unwrap();
    /// ```
    fn calc_charm(&self) -> Result<f64, BlackScholesError> {
        let sigma = self.sigma.ok_or(BlackScholesError::ExpectedSigma)?;
        let nprimed1 = calc_nprimed1(self)?;
        let (nd1, _) = calc_nd1nd2(self)?;
        let (_, d2) = calc_d1d2(self)?;
        let e_negqt = self.calc_e_negqt();

        Ok(self.calc_charm_with_params(
            sigma,
            nprimed1,
            nd1,
            d2,
            e_negqt,
            self.q,
            self.r,
            self.t,
            self.option_type,
        ))
    }

    /// Calculates the veta of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f64 of the veta of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let veta = inputs.calc_veta().unwrap();
    /// ```
    fn calc_veta(&self) -> Result<f64, BlackScholesError> {
        let sigma = self.sigma.ok_or(BlackScholesError::ExpectedSigma)?;
        let nprimed1 = calc_nprimed1(self)?;
        let (d1, d2) = calc_d1d2(self)?;
        let e_negqt = self.calc_e_negqt();

        Ok(self.calc_veta_with_params(
            self.s, e_negqt, nprimed1, self.t, self.q, self.r, d1, d2, sigma,
        ))
    }

    /// Calculates the vomma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f64 of the vomma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let vomma = inputs.calc_vomma().unwrap();
    /// ```
    fn calc_vomma(&self) -> Result<f64, BlackScholesError> {
        let sigma = self.sigma.ok_or(BlackScholesError::ExpectedSigma)?;
        let (d1, d2) = calc_d1d2(self)?;
        let nprimed1: f64 = calc_nprimed1(self)?;

        let vega = self.calc_vega_with_params(nprimed1, self.s, self.q, self.t);

        Ok(self.calc_vomma_with_params(vega, d1, d2, sigma))
    }

    /// Calculates the speed of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f64 of the speed of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let speed = inputs.calc_speed().unwrap();
    /// ```
    fn calc_speed(&self) -> Result<f64, BlackScholesError> {
        let sigma = self.sigma.ok_or(BlackScholesError::ExpectedSigma)?;
        let (d1, _) = calc_d1d2(self)?;
        let gamma = Inputs::calc_gamma(self)?;

        Ok(self.calc_speed_with_params(gamma, d1, sigma, self.t))
    }

    /// Calculates the zomma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f64 of the zomma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let zomma = inputs.calc_zomma().unwrap();
    /// ```
    fn calc_zomma(&self) -> Result<f64, BlackScholesError> {
        let sigma = self.sigma.ok_or(BlackScholesError::ExpectedSigma)?;
        let (d1, d2) = calc_d1d2(self)?;
        let gamma = Inputs::calc_gamma(self)?;

        Ok(self.calc_zomma_with_params(gamma, d1, d2, sigma))
    }

    /// Calculates the color of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f64 of the color of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let color = inputs.calc_color().unwrap();
    /// ```
    fn calc_color(&self) -> Result<f64, BlackScholesError> {
        let sigma = self.sigma.ok_or(BlackScholesError::ExpectedSigma)?;
        let (d1, d2) = calc_d1d2(self)?;
        let nprimed1 = calc_nprimed1(self)?;
        let e_negqt = self.calc_e_negqt();

        let color = self.calc_color_with_params(sigma, d1, d2, nprimed1, e_negqt);
        Ok(color)
    }

    /// Calculates the ultima of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f64 of the ultima of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let ultima = inputs.calc_ultima().unwrap();
    /// ```
    fn calc_ultima(&self) -> Result<f64, BlackScholesError> {
        let sigma = self.sigma.ok_or(BlackScholesError::ExpectedSigma)?;
        let (d1, d2) = calc_d1d2(self)?;
        let vega = Inputs::calc_vega(self)?;

        let ultima = self.calc_ultima_with_params(vega, sigma, d1, d2);

        Ok(ultima)
    }

    /// Calculates the dual delta of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f64 of the dual delta of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let dual_delta = inputs.calc_dual_delta().unwrap();
    /// ```
    fn calc_dual_delta(&self) -> Result<f64, BlackScholesError> {
        let (_, nd2) = calc_nd1nd2(self)?;
        let e_negqt = self.calc_e_negqt();

        let dual_delta = self.calc_dual_delta_with_params(nd2, e_negqt);

        Ok(dual_delta)
    }

    /// Calculates the dual gamma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f64 of the dual gamma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let dual_gamma = inputs.calc_dual_gamma().unwrap();
    /// ```
    fn calc_dual_gamma(&self) -> Result<f64, BlackScholesError> {
        let sigma = self.sigma.ok_or(BlackScholesError::ExpectedSigma)?;
        let nprimed2 = calc_nprimed2(self)?;
        let e_negqt = self.calc_e_negqt();

        let dual_gamma = self.calc_dual_gamma_with_params(sigma, nprimed2, e_negqt);

        Ok(dual_gamma)
    }

    /// Calculates all Greeks of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// HashMap of type <Greek, f64> of all Greeks of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let greeks = inputs.calc_all_greeks().unwrap();
    /// ```
    fn calc_all_greeks(&self) -> Result<HashMap<Greek, f64>, BlackScholesError> {
        let sigma = self.sigma.ok_or(BlackScholesError::ExpectedSigma)?;
        let (nprimed1, nprimed2) = calc_nprimed_1_2(self)?;
        let (nd1, nd2) = calc_nd1nd2(self)?;
        let e_negqt = self.calc_e_negqt();
        let (d1, d2) = calc_d1d2(self)?;

        let mut greeks: HashMap<Greek, f64> = HashMap::with_capacity(17);
        let delta = self.calc_delta_with_params(nd1, self.q, self.t, self.option_type);
        let gamma = self.calc_gamma_with_params(nprimed1, self.s, sigma, self.q, self.t);
        let vega = self.calc_vega_with_params(nprimed1, self.s, self.q, self.t);
        let price = self.calc_price_with_params(
            nd1,
            nd2,
            self.s,
            self.q,
            self.t,
            self.k,
            self.r,
            self.option_type,
        );

        greeks.insert(Greek::Delta, delta);
        greeks.insert(Greek::Gamma, gamma);
        greeks.insert(
            Greek::Theta,
            self.calc_theta_with_params(
                sigma,
                nprimed1,
                nd1,
                nd2,
                self.s,
                self.q,
                self.t,
                self.r,
                self.k,
                self.option_type,
            ),
        );
        greeks.insert(Greek::Vega, vega);
        greeks.insert(
            Greek::Rho,
            self.calc_rho_with_params(self.k, self.t, self.r, nd2, self.option_type),
        );
        greeks.insert(
            Greek::Epsilon,
            self.calc_epsilon_with_params(self.s, self.t, e_negqt, nd1, self.option_type),
        );
        greeks.insert(
            Greek::Lambda,
            self.calc_lambda_with_params(delta, self.s, price),
        );
        greeks.insert(
            Greek::Vanna,
            self.calc_vanna_with_params(d2, e_negqt, nprimed1, sigma),
        );
        greeks.insert(
            Greek::Charm,
            self.calc_charm_with_params(
                sigma,
                nprimed1,
                nd1,
                d2,
                e_negqt,
                self.q,
                self.r,
                self.t,
                self.option_type,
            ),
        );
        greeks.insert(
            Greek::Veta,
            self.calc_veta_with_params(
                self.s, e_negqt, nprimed1, self.t, self.q, self.r, d1, d2, sigma,
            ),
        );
        greeks.insert(
            Greek::Vomma,
            self.calc_vomma_with_params(vega, d1, d2, sigma),
        );
        greeks.insert(
            Greek::Speed,
            self.calc_speed_with_params(gamma, d1, sigma, self.t),
        );
        greeks.insert(
            Greek::Zomma,
            self.calc_zomma_with_params(gamma, d1, d2, sigma),
        );
        greeks.insert(
            Greek::Color,
            self.calc_color_with_params(sigma, d1, d2, nprimed1, e_negqt),
        );
        greeks.insert(
            Greek::Ultima,
            self.calc_ultima_with_params(vega, sigma, d1, d2),
        );
        greeks.insert(
            Greek::DualDelta,
            self.calc_dual_delta_with_params(nd2, e_negqt),
        );
        greeks.insert(
            Greek::DualGamma,
            self.calc_dual_gamma_with_params(sigma, nprimed2, e_negqt),
        );

        Ok(greeks)
    }
}

impl Inputs {
    fn calc_dual_gamma_with_params(&self, sigma: f64, nprimed2: f64, e_negqt: f64) -> f64 {
        e_negqt * (nprimed2 / (self.k * sigma * self.t.sqrt()))
    }

    fn calc_dual_delta_with_params(&self, nd2: f64, e_negqt: f64) -> f64 {
        match self.option_type {
            OptionType::Call => -e_negqt * nd2,
            OptionType::Put => e_negqt * nd2,
        }
    }

    fn calc_ultima_with_params(&self, vega: f64, sigma: f64, d1: f64, d2: f64) -> f64 {
        -vega / sigma.powi(2) * (d1 * d2 * (1.0 - d1 * d2) + d1.powi(2) + d2.powi(2))
    }

    fn calc_color_with_params(
        &self,
        sigma: f64,
        d1: f64,
        d2: f64,
        nprimed1: f64,
        e_negqt: f64,
    ) -> f64 {
        -e_negqt
            * (nprimed1 / (2.0 * self.s * self.t * sigma * self.t.sqrt()))
            * (2.0 * self.q * self.t
                + 1.0
                + (2.0 * (self.r - self.q) * self.t - d2 * sigma * self.t.sqrt())
                    / (sigma * self.t.sqrt())
                    * d1)
    }

    fn calc_zomma_with_params(&self, gamma: f64, d1: f64, d2: f64, sigma: f64) -> f64 {
        gamma * ((d1 * d2 - 1.0) / sigma)
    }
    fn calc_speed_with_params(&self, gamma: f64, d1: f64, sigma: f64, t: f64) -> f64 {
        -gamma / self.s * (d1 / (sigma * t.sqrt()) + 1.0)
    }

    fn calc_vega_with_params(&self, nprimed1: f64, s: f64, q: f64, t: f64) -> f64 {
        PERCENTAGE_CONVERSION * s * (-q * t).exp() * t.sqrt() * nprimed1
    }
    fn calc_delta_with_params(&self, nd1: f64, q: f64, t: f64, option_type: OptionType) -> f64 {
        option_type * (-q * t).exp() * nd1
    }

    fn calc_gamma_with_params(&self, nprimed1: f64, s: f64, sigma: f64, q: f64, t: f64) -> f64 {
        (-q * t).exp() * nprimed1 / (s * sigma * t.sqrt())
    }

    fn calc_theta_with_params(
        &self,
        sigma: f64,
        nprimed1: f64,
        nd1: f64,
        nd2: f64,
        s: f64,
        q: f64,
        t: f64,
        r: f64,
        k: f64,
        option_type: OptionType,
    ) -> f64 {
        let term1 = -(s * sigma * (-q * t).exp() * nprimed1 / (2.0 * t.sqrt()));
        let term2 = r * k * (-r * t).exp() * nd2 * option_type;
        let term3 = q * s * (-q * t).exp() * nd1 * option_type;

        (term1 - term2 + term3) / DAYS_PER_YEAR
    }
    fn calc_rho_with_params(
        &self,
        k: f64,
        t: f64,
        r: f64,
        nd2: f64,
        option_type: OptionType,
    ) -> f64 {
        option_type * k * t * (-r * t).exp() * nd2 / 100.0
    }

    fn calc_epsilon_with_params(
        &self,
        s: f64,
        t: f64,
        e_negqt: f64,
        nd1: f64,
        option_type: OptionType,
    ) -> f64 {
        -s * t * e_negqt * nd1 * option_type
    }
    fn calc_lambda_with_params(&self, delta: f64, s: f64, price: f64) -> f64 {
        delta * s / price
    }

    pub(crate) fn calc_price_with_params(
        &self,
        nd1: f64,
        nd2: f64,
        s: f64,
        q: f64,
        t: f64,
        k: f64,
        r: f64,
        option_type: OptionType,
    ) -> f64 {
        match option_type {
            OptionType::Call => f64::max(0.0, nd1 * s * (-q * t).exp() - nd2 * k * (-r * t).exp()),
            OptionType::Put => f64::max(0.0, nd2 * k * (-r * t).exp() - nd1 * s * (-q * t).exp()),
        }
    }
    fn calc_vanna_with_params(&self, d2: f64, e_negqt: f64, nprimed1: f64, sigma: f64) -> f64 {
        d2 * e_negqt * nprimed1 * -PERCENTAGE_CONVERSION / sigma
    }
    fn calc_charm_with_params(
        &self,
        sigma: f64,
        nprimed1: f64,
        nd1: f64,
        d2: f64,
        e_negqt: f64,
        q: f64,
        r: f64,
        t: f64,
        option_type: OptionType,
    ) -> f64 {
        option_type * q * e_negqt * nd1
            - e_negqt * nprimed1 * (2.0 * (r - q) * t - d2 * sigma * t.sqrt())
                / (2.0 * t * sigma * t.sqrt())
    }

    fn calc_veta_with_params(
        &self,
        s: f64,
        e_negqt: f64,
        nprimed1: f64,
        t: f64,
        q: f64,
        r: f64,
        d1: f64,
        d2: f64,
        sigma: f64,
    ) -> f64 {
        -s * e_negqt
            * nprimed1
            * t.sqrt()
            * (q + ((r - q) * d1) / (sigma * t.sqrt()) - ((1.0 + d1 * d2) / (2.0 * t)))
    }
    fn calc_vomma_with_params(&self, vega: f64, d1: f64, d2: f64, sigma: f64) -> f64 {
        vega * ((d1 * d2) / sigma)
    }

    fn calc_e_negqt(&self) -> f64 {
        (-self.q * self.t).exp()
    }
}
