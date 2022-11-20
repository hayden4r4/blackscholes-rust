from enum import Enum, auto
from typing import Optional

class OptionType(Enum):
    Call = auto()
    Put = auto()
    
class Inputs:
    """
    A class containing the inputs to the Black-Scholes model
    """
    def __new__(cls, option_type: OptionType, s: float, k: float, p: Optional[float], r: float, q: float, t: float, sigma: Optional[float]) -> Inputs: ...
    def __str__(self) -> str: ...

    def calc_price(self) -> float: ...
    """
    Calculates the price of the option
    """

    def calc_delta(self) -> float: ...
    """
    Calculates the delta of the option
    """

    def calc_gamma(self) -> float: ...
    """
    Calculates the gamma of the option
    """
  
    def calc_vega(self) -> float: ...
    """
    Calculates the vega of the option
    """
  
    def calc_theta(self) -> float: ...
    """
    Calculates the theta of the option
    """
  
    def calc_rho(self) -> float: ...
    """
    Calculates the rho of the option
    """

    def calc_iv(self, tolerance: float) -> float: ...
    """
    Calculates the implied volatility of the option
    """