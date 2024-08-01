use thiserror::Error;

#[derive(Error, Debug)]
pub enum BlackScholesError {
    #[error("Expected Some(f64) for self.sigma, received None")]
    ExpectedSigma,

    #[error("Log from s/k is infinity")]
    LogSdivKInfinity,

    #[error("Expected Some(f64) for self.p, received None")]
    ExpectedPrice,

    #[error("Implied volatility failed to converge")]
    FailedToConverge,

    #[error("Time to maturity is 0")]
    ZeroTimeToMaturity,
}
