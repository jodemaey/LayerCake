
from abc import ABC


class Variable(ABC):

    def __init__(self, name, symbol, units=None):

        self.name = name
        self.symbol = symbol
        if units is None:
            self.units = ""
        else:
            self.units = units


