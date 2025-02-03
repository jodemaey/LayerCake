
from little.arithmetic.terms.base import ArithmeticTerm


class LinearTerm(ArithmeticTerm):

    def __init__(self):

        self.name = 'Linear term'
        self.inner_products = None
        self.variables = None