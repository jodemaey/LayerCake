
"""
    Inner products definition module
    ================================

    Module containing classes to define the `inner products`_ used by the model.

    .. _inner products: https://en.wikipedia.org/wiki/Inner_product_space
    
"""

from abc import ABC, abstractmethod


class InnerProductDefinition(ABC):
    """Base class to define the model's basis inner products.

    Parameters
    ----------
    variables: list(~sympy.core.symbol.Symbol)
        List of variables used as coordinates.
    optimizer: None or callable, optional
        A function to optimize the computation of the integrals or the integrand.
        If `None`, does not optimize.

    Attributes
    ----------
    variables: list(~sympy.core.symbol.Symbol)
        List of variables used as coordinates.
    optimizer: None or callable
        A function to optimize the computation of the integrals or the integrand.
        If `None`, does not optimize the computation.
    """

    def __init__(self, variables=None, optimizer=None):

        # must warn user if None
        self.variables = variables

        self.optimizer = None

        if optimizer is not None:
            self.set_optimizer(optimizer)
        else:
            self.set_optimizer(self._no_optimizer)

    def set_optimizer(self, optimizer):
        """Function to set the optimizer.

        Parameters
        ----------
        optimizer: callable
            A function to optimize the computation of the integrals or the integrand.
        """
        self.optimizer = optimizer

    @staticmethod
    def _no_optimizer(expr):
        return expr

    @abstractmethod
    def inner_product(self, S, G, symbolic_expr=False, integrand=False):
        """Symbolic definition of the inner product :math:`(S, G)`.

        Parameters
        ----------
        S: Sympy expression
            Left-hand side function of the product.
        G: Sympy expression
            Right-hand side function of the product.
        symbolic_expr: bool, optional
            If `True`, return the integral as a symbolic expression object. Else, return the integral performed symbolically.
        integrand: bool, optional
            If `True`, return the integrand of the integral and its integration limits as a list of symbolic expression object. Else, return the integral performed symbolically.

        Returns
        -------
        Sympy expression
            The symbolic result of the inner product.
        """
        pass
