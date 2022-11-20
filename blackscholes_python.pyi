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
class Price:
    """
    A class containing methods to calculate price
    """
    @staticmethod
    def calc_price(inp: Inputs) -> float: ...
class Greeks:
    """
    A class containing methods to calculate greeks
    """
    @staticmethod
    def calc_delta(inp: Inputs) -> float: ...
    @staticmethod
    def calc_gamma(inp: Inputs) -> float: ...
    @staticmethod
    def calc_vega(inp: Inputs) -> float: ...
    @staticmethod
    def calc_theta(inp: Inputs) -> float: ...
    @staticmethod
    def calc_rho(inp: Inputs) -> float: ...
class Volatility:
    """
    A class containing methods to calculate implied volatility
    """
    @staticmethod
    def calc_iv(inp: Inputs, tolerance: float) -> float: ...