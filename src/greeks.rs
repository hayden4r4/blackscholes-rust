use num_traits::Float;

use crate::{Inputs, OptionType, Pricing, *};

// Static error message to avoid allocation
static ERR_MISSING_SIGMA: &str = "Expected Some(f32) for self.sigma, received None";

/// Contains all calculated Greek values for an option
#[derive(Debug, Clone)]
pub struct GreekResults {
    pub delta: f32,
    pub gamma: f32,
    pub theta: f32,
    pub vega: f32,
    pub rho: f32,
    pub epsilon: f32,
    pub lambda: f32,
    pub vanna: f32,
    pub charm: f32,
    pub veta: f32,
    pub vomma: f32,
    pub speed: f32,
    pub zomma: f32,
    pub color: f32,
    pub ultima: f32,
    pub dual_delta: f32,
    pub dual_gamma: f32,
}

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
    fn calc_all_greeks(&self) -> Result<GreekResults, String>;
}

impl Greeks<f32> for Inputs {
    /// Calculates the delta of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the delta of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let delta = inputs.calc_delta().unwrap();
    /// ```
    #[inline]
    fn calc_delta(&self) -> Result<f32, String> {
        let (nd1, _): (f32, f32) = calc_nd1nd2(self)?;

        let option_type: f64 = self.option_type.into();
        let delta = option_type as f32 * E.powf(-self.q * self.t) * nd1;

        Ok(delta)
    }

    /// Calculates the gamma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the gamma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let gamma = inputs.calc_gamma().unwrap();
    /// ```
    #[inline]
    fn calc_gamma(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;

        let nprimed1: f32 = calc_nprimed1(self)?;
        let gamma: f32 = E.powf(-self.q * self.t) * nprimed1 / (self.s * sigma * self.t.sqrt());
        Ok(gamma)
    }

    /// Calculates the theta of the option.
    /// Uses 365.25 days in a year for calculations.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of theta per day (not per year).
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let theta = inputs.calc_theta().unwrap();
    /// ```
    #[inline]
    fn calc_theta(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;

        let nprimed1: f32 = calc_nprimed1(self)?;
        let (nd1, nd2): (f32, f32) = calc_nd1nd2(self)?;

        // Calculation uses 365.25 for f32: Time of days per year.
        let option_type: f64 = self.option_type.into();
        let theta = (-(self.s * sigma * E.powf(-self.q * self.t) * nprimed1
            / (2.0 * self.t.sqrt()))
            - self.r * self.k * E.powf(-self.r * self.t) * nd2 * option_type as f32
            + self.q * self.s * E.powf(-self.q * self.t) * nd1 * option_type as f32)
            / DAYS_PER_YEAR;

        Ok(theta)
    }

    /// Calculates the vega of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the vega of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let vega = inputs.calc_vega().unwrap();
    /// ```
    #[inline]
    fn calc_vega(&self) -> Result<f32, String> {
        let nprimed1: f32 = calc_nprimed1(self)?;
        let vega: f32 = 0.01 * self.s * E.powf(-self.q * self.t) * self.t.sqrt() * nprimed1;
        Ok(vega)
    }

    /// Calculates the rho of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the rho of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let rho = inputs.calc_rho().unwrap();
    /// ```
    #[inline]
    fn calc_rho(&self) -> Result<f32, String> {
        let (_, nd2): (f32, f32) = calc_nd1nd2(self)?;

        let option_type: f64 = self.option_type.into();
        let rho = option_type as f32 / 100.0 * self.k * self.t * E.powf(-self.r * self.t) * nd2;

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
    /// f32 of the epsilon of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let epsilon = inputs.calc_epsilon().unwrap();
    /// ```
    fn calc_epsilon(&self) -> Result<f32, String> {
        let (nd1, _) = calc_nd1nd2(self)?;
        let e_negqt = E.powf(-self.q * self.t);

        let option_type: f64 = self.option_type.into();
        let epsilon: f32 = -self.s * self.t * e_negqt * nd1 * option_type as f32;

        Ok(epsilon)
    }

    /// Calculates the lambda of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the lambda of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks, Pricing};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let lambda = inputs.calc_lambda().unwrap();
    /// ```
    fn calc_lambda(&self) -> Result<f32, String> {
        let delta = self.calc_delta()?;
        Ok(delta * self.s / self.calc_price()?)
    }

    /// Calculates the vanna of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the vanna of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let vanna = inputs.calc_vanna().unwrap();
    /// ```
    fn calc_vanna(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;

        let nprimed1 = calc_nprimed1(self)?;
        let (_, d2) = calc_d1d2(self)?;
        let vanna: f32 = d2 * E.powf(-self.q * self.t) * nprimed1 * -0.01 / sigma;
        Ok(vanna)
    }

    // /// Calculates the charm of the option.
    // /// # Requires
    // /// s, k, r, q, t, sigma
    // /// # Returns
    // /// f32 of the charm of the option.
    // /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let charm = inputs.calc_charm().unwrap();
    /// ```
    fn calc_charm(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;
        let nprimed1 = calc_nprimed1(self)?;
        let (nd1, _) = calc_nd1nd2(self)?;
        let (_, d2) = calc_d1d2(self)?;
        let e_negqt = E.powf(-self.q * self.t);

        let option_type: f64 = self.option_type.into();
        let charm: f32 = option_type as f32 * self.q * e_negqt * nd1
            - e_negqt * nprimed1 * (2.0 * (self.r - self.q) * self.t - d2 * sigma * self.t.sqrt())
                / (2.0 * self.t * sigma * self.t.sqrt());

        Ok(charm)
    }

    /// Calculates the veta of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the veta of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let veta = inputs.calc_veta().unwrap();
    /// ```
    fn calc_veta(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;
        let nprimed1 = calc_nprimed1(self)?;
        let (d1, d2) = calc_d1d2(self)?;
        let e_negqt = E.powf(-self.q * self.t);

        let veta = -self.s
            * e_negqt
            * nprimed1
            * self.t.sqrt()
            * (self.q + ((self.r - self.q) * d1) / (sigma * self.t.sqrt())
                - ((1.0 + d1 * d2) / (2.0 * self.t)));
        Ok(veta)
    }

    /// Calculates the vomma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the vomma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let vomma = inputs.calc_vomma().unwrap();
    /// ```
    fn calc_vomma(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;
        let (d1, d2) = calc_d1d2(self)?;

        let vomma = Inputs::calc_vega(self)? * ((d1 * d2) / sigma);
        Ok(vomma)
    }

    /// Calculates the speed of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the speed of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let speed = inputs.calc_speed().unwrap();
    /// ```
    fn calc_speed(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;
        let (d1, _) = calc_d1d2(self)?;
        let gamma = Inputs::calc_gamma(self)?;

        let speed = -gamma / self.s * (d1 / (sigma * self.t.sqrt()) + 1.0);
        Ok(speed)
    }

    /// Calculates the zomma of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the zomma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let zomma = inputs.calc_zomma().unwrap();
    /// ```
    fn calc_zomma(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;
        let (d1, d2) = calc_d1d2(self)?;
        let gamma = Inputs::calc_gamma(self)?;

        let zomma = gamma * ((d1 * d2 - 1.0) / sigma);
        Ok(zomma)
    }

    /// Calculates the color of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the color of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let color = inputs.calc_color().unwrap();
    /// ```
    fn calc_color(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;
        let (d1, d2) = calc_d1d2(self)?;
        let nprimed1 = calc_nprimed1(self)?;
        let e_negqt = E.powf(-self.q * self.t);

        let color = -e_negqt
            * (nprimed1 / (2.0 * self.s * self.t * sigma * self.t.sqrt()))
            * (2.0 * self.q * self.t
                + 1.0
                + (2.0 * (self.r - self.q) * self.t - d2 * sigma * self.t.sqrt())
                    / (sigma * self.t.sqrt())
                    * d1);
        Ok(color)
    }

    /// Calculates the ultima of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the ultima of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let ultima = inputs.calc_ultima().unwrap();
    /// ```
    fn calc_ultima(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;
        let (d1, d2) = calc_d1d2(self)?;
        let vega = Inputs::calc_vega(self)?;

        let ultima =
            -vega / sigma.powf(2.0) * (d1 * d2 * (1.0 - d1 * d2) + d1.powf(2.0) + d2.powf(2.0));
        Ok(ultima)
    }

    /// Calculates the dual delta of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// f32 of the dual delta of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let dual_delta = inputs.calc_dual_delta().unwrap();
    /// ```
    fn calc_dual_delta(&self) -> Result<f32, String> {
        let (_, nd2) = calc_nd1nd2(self)?;
        let e_negqt = E.powf(-self.q * self.t);

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
    /// f32 of the dual gamma of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let dual_gamma = inputs.calc_dual_gamma().unwrap();
    /// ```
    fn calc_dual_gamma(&self) -> Result<f32, String> {
        let sigma = self
            .sigma
            .ok_or("Expected Some(f32) for self.sigma, received None")?;
        let nprimed2 = calc_nprimed2(self)?;
        let e_negqt = E.powf(-self.q * self.t);

        let dual_gamma = e_negqt * (nprimed2 / (self.k * sigma * self.t.sqrt()));
        Ok(dual_gamma)
    }

    /// Calculates all Greeks of the option.
    /// # Requires
    /// s, k, r, q, t, sigma
    /// # Returns
    /// GreekResults struct containing all Greeks of the option.
    /// # Example
    /// ```
    /// use blackscholes::{Inputs, OptionType, Greeks};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 0.2, 20.0/365.25, Some(0.2));
    /// let greeks = inputs.calc_all_greeks().unwrap();
    /// ```
    fn calc_all_greeks(&self) -> Result<GreekResults, String> {
        // Validate and extract sigma once
        let sigma = self.sigma.ok_or(ERR_MISSING_SIGMA)?;
        
        // Pre-calculate common values used across multiple Greeks
        let sqrt_t = self.t.sqrt();
        let e_negrt = E.powf(-self.r * self.t);
        let e_negqt = E.powf(-self.q * self.t);
        
        // Calculate d1 and d2 once
        let (d1, d2) = calc_d1d2(self)?;
        
        // Calculate normal distribution values once
        let (nd1, nd2) = calc_nd1nd2(self)?;
        
        // Calculate derivatives once
        let nprimed1 = calc_npdf(d1);
        let nprimed2 = calc_npdf(d2);
        
        // Calculate common terms for multiple Greeks
        let option_type_sign = match self.option_type {
            OptionType::Call => 1.0f32,
            OptionType::Put => -1.0f32,
        };
        
        // Calculate Greeks using cached values
        
        // Delta
        let delta = option_type_sign * e_negqt * nd1;
        
        // Gamma
        let gamma = e_negqt * nprimed1 / (self.s * sigma * sqrt_t);
        
        // Vega (scaled by 0.01)
        let vega = 0.01 * self.s * e_negqt * sqrt_t * nprimed1;
        
        // Rho (scaled by 0.01)
        let rho = option_type_sign / 100.0 * self.k * self.t * e_negrt * nd2;
        
        // Theta (daily)
        let theta = (-(self.s * sigma * e_negqt * nprimed1 / (2.0 * sqrt_t))
            - self.r * self.k * e_negrt * nd2 * option_type_sign
            + self.q * self.s * e_negqt * nd1 * option_type_sign)
            / DAYS_PER_YEAR;
        
        // Epsilon
        let epsilon = -self.s * self.t * e_negqt * nd1 * option_type_sign;
        
        // Lambda
        let lambda = delta * self.s / Inputs::calc_price(self)?;
        
        // Vanna
        let vanna = -e_negqt * nprimed1 * d2 / sigma * 0.01;
        
        // Charm
        let charm = -(self.q * e_negqt * nd1 * option_type_sign
            - e_negqt * nprimed1 * (2.0 * (self.r - self.q) * sqrt_t - d2 * sigma * sqrt_t)
                / (2.0 * self.t * sigma * sqrt_t))
            / DAYS_PER_YEAR;
        
        // Veta
        let veta = -self.s * e_negqt * nprimed1 * sqrt_t
            * (self.q + (self.r - self.q) * d1 / (sigma * sqrt_t)
                - (1.0 + d1 * d2) / (2.0 * self.t))
            * 0.01 / DAYS_PER_YEAR;
        
        // Vomma
        let vomma = vega * (d1 * d2 / sigma);
        
        // Speed
        let speed = -gamma / self.s * (d1 / (sigma * sqrt_t) + 1.0);
        
        // Zomma
        let zomma = gamma * ((d1 * d2 - 1.0) / sigma);
        
        // Color
        let color = -e_negqt
            * (nprimed1 / (2.0 * self.s * self.t * sigma * sqrt_t))
            * (2.0 * self.q * self.t
                + 1.0
                + (2.0 * (self.r - self.q) * self.t - d2 * sigma * sqrt_t)
                    / (sigma * sqrt_t)
                    * d1);
        
        // Ultima
        let ultima = -vega / sigma.powf(2.0) * (d1 * d2 * (1.0 - d1 * d2) + d1.powf(2.0) + d2.powf(2.0));
        
        // Dual Delta
        let dual_delta = match self.option_type {
            OptionType::Call => -e_negqt * nd2,
            OptionType::Put => e_negqt * nd2,
        };
        
        // Dual Gamma
        let dual_gamma = e_negqt * (nprimed2 / (self.k * sigma * sqrt_t));

        Ok(GreekResults {
            delta,
            gamma,
            theta,
            vega,
            rho,
            epsilon,
            lambda,
            vanna,
            charm,
            veta,
            vomma,
            speed,
            zomma,
            color,
            ultima,
            dual_delta,
            dual_gamma,
        })
    }
}
