import unittest
import blackscholes_python as bs

class Tests(unittest.TestCase):
    def setUp(self):
        self.tolerance = 0.02
        self.inputs_put1 = bs.Inputs(option_type=bs.OptionType.Put, s=100, k=110, p=10, r=0.05, q=0.05, t=20/365.25, sigma=0.2)
        self.inputs_put2 = bs.Inputs(option_type=bs.OptionType.Put, s=110, k=100, p=10, r=0.05, q=0.05, t=20/365.25, sigma=0.2)
        self.inputs_call1 = bs.Inputs(option_type=bs.OptionType.Call, s=100, k=110, p=10, r=0.05, q=0.05, t=20/365.25, sigma=0.2)
        self.inputs_call2 = bs.Inputs(option_type=bs.OptionType.Call, s=110, k=100, p=10, r=0.05, q=0.05, t=20/365.25, sigma=0.2)

    def testCalcPrice(self):
        price_put1 = bs.Price.calc_price(self.inputs_put1)
        self.assertTrue(abs(10 - price_put1) <= self.tolerance)
        price_put2 = bs.Price.calc_price(self.inputs_put2)
        self.assertTrue(abs(0.04 - price_put2) <= self.tolerance)
        price_call1 = bs.Price.calc_price(self.inputs_call1)
        self.assertTrue(abs(0.04 - price_call1) <= self.tolerance)
        price_call2 = bs.Price.calc_price(self.inputs_call2)
        self.assertTrue(abs(10 - price_call2) <= self.tolerance)

    def testCalcGreeksDelta(self) -> None:
        delta_put1 = bs.Greeks.calc_delta(self.inputs_put1)
        self.assertTrue(abs(-0.975 - delta_put1) <= self.tolerance)
        delta_put2 = bs.Greeks.calc_delta(self.inputs_put2)
        self.assertTrue(abs(-0.02 - delta_put2) <= self.tolerance)
        delta_call1 = bs.Greeks.calc_delta(self.inputs_call1)
        self.assertTrue(abs(0.022 - delta_call1) <= self.tolerance)
        delta_call2 = bs.Greeks.calc_delta(self.inputs_call2)
        self.assertTrue(abs(0.978 - delta_call2) <= self.tolerance)

    def testCalcGreeksGamma(self) -> None:
        gamma_put1 = bs.Greeks.calc_gamma(self.inputs_put1)
        self.assertTrue(abs(0.011 - gamma_put1) <= self.tolerance)
        gamma_put2 = bs.Greeks.calc_gamma(self.inputs_put2)
        self.assertTrue(abs(0.009 - gamma_put2) <= self.tolerance)
        gamma_call1 = bs.Greeks.calc_gamma(self.inputs_call1)
        self.assertTrue(abs(0.011 - gamma_call1) <= self.tolerance)
        gamma_call2 = bs.Greeks.calc_gamma(self.inputs_call2)
        self.assertTrue(abs(0.009 - gamma_call2) <= self.tolerance)

    def testCalcGreeksVega(self) -> None:
        vega_put1 = bs.Greeks.calc_vega(self.inputs_put1)
        self.assertTrue(abs(0.01229 - vega_put1) <= self.tolerance)
        vega_put2 = bs.Greeks.calc_vega(self.inputs_put2)
        self.assertTrue(abs(0.01229 - vega_put2) <= self.tolerance)
        vega_call1 = bs.Greeks.calc_vega(self.inputs_call1)
        self.assertTrue(abs(0.01229 - vega_call1) <= self.tolerance)
        vega_call2 = bs.Greeks.calc_vega(self.inputs_call2)
        self.assertTrue(abs(0.01229 - vega_call2) <= self.tolerance)
    
    def testCalcGreeksTheta(self) -> None:
        theta_put1 = bs.Greeks.calc_theta(self.inputs_put1)
        self.assertTrue(abs(-0.0047 - theta_put1) <= self.tolerance)
        theta_put2 = bs.Greeks.calc_theta(self.inputs_put2)
        self.assertTrue(abs(-0.006 - theta_put2) <= self.tolerance)
        theta_call1 = bs.Greeks.calc_theta(self.inputs_call1)
        self.assertTrue(abs(-0.006 - theta_call1) <= self.tolerance)
        theta_call2 = bs.Greeks.calc_theta(self.inputs_call2)
        self.assertTrue(abs(-0.0047 - theta_call2) <= self.tolerance)
    
    def testCalcGreeksRho(self) -> None:
        rho_put1 = bs.Greeks.calc_rho(self.inputs_put1)
        self.assertTrue(abs(-0.0589 - rho_put1) <= self.tolerance)
        rho_put2 = bs.Greeks.calc_rho(self.inputs_put2)
        self.assertTrue(abs(-0.0012 - rho_put2) <= self.tolerance)
        rho_call1 = bs.Greeks.calc_rho(self.inputs_call1)
        self.assertTrue(abs(0.0012 - rho_call1) <= self.tolerance)
        rho_call2 = bs.Greeks.calc_rho(self.inputs_call2)
        self.assertTrue(abs(0.0534 - rho_call2) <= self.tolerance)

    def testCalcIV(self):
        iv_put1 = bs.Volatility.calc_iv(self.inputs_put1, tolerance=0.0001)
        self.assertTrue(abs(0.1907 - iv_put1) <= self.tolerance)
        iv_put2 = bs.Volatility.calc_iv(self.inputs_put2, tolerance=0.0001)
        self.assertTrue(abs(1.4854 - iv_put2) <= self.tolerance)
        iv_call1 = bs.Volatility.calc_iv(self.inputs_call1, tolerance=0.0001)
        self.assertTrue(abs(1.4854 - iv_call1) <= self.tolerance)
        iv_call2 = bs.Volatility.calc_iv(self.inputs_call2, tolerance=0.0001)
        self.assertTrue(abs(0.1907 - iv_call2) <= self.tolerance)


if __name__ == '__main__':
    unittest.main()