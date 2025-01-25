
from abc import ABC


class Variable(ABC):
    pass


class Coordinate(Variable):
    """Class to define a physical coordinate in the system.

    Parameters
    ----------
    symbol: ~sympy.core.symbol.Symbol
        Sympy symbol of the coordinate
    extent: tuple(float)
        The natural extent of the coordinate.

    Attributes
    ----------
    symbol: ~sympy.core.symbol.Symbol
        Sympy symbol of the coordinate
    extent: tuple(float)
        The natural extent of the coordinate.

    Warning
    -------
    Coordinates with infinite extent are not currently supported.

    """

    def __init__(self, symbol, extent):

        self.symbol = symbol
        self.extent = extent

