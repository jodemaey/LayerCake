"""
    Fourier Basis definition module
    ===============================

    Classes and functions defining Fourier basis of functions (`Fourier modes`_) on a plane.

    Description of the classes
    --------------------------

    * :class:`PlanarChannelFourierBasis`: Fourier basis defined on a zonally perdiodic channel, with no-flux boundary conditions in the meridional direction :math:`y`.
    * :class:`PlanarBasinFourierBasis`: Fourier basis defined on a closed basin, with no-flux boundary conditions in both the zonal and meridional direction :math:`x` and :math:`y`.

    .. _Fourier modes: https://en.wikipedia.org/wiki/Fourier_series

"""
import numpy as np

from layercake.basis.base import SymbolicBasis
from layercake.variables.systems import PlanarCartesianCoordinateSystem
from sympy import symbols, sin, cos, sqrt


class PlanarChannelFourierBasis(SymbolicBasis):
    """Fourier basis defined on a zonally periodic channel on a 2D plane, with no-flux boundary conditions in the meridional
    direction :math:`y`.

    Parameters
    ----------
    spectral_blocks: ~numpy.ndarray(int)
        Spectral blocks detailing the modes :math:`x`- and :math:`y`-wavenumber.
        Array of shape (nblocks, 2), where `nblocks` is the number of spectral blocks.
    aspect_ratio: float
        Spatial domain aspect ratio, :math:`n = 2 L_y/L_x` .
    length: float, optional
        Length of the domain along the :math:`x` coordinate: :math:`L_x` .
    """

    _n = symbols('n', positive=True)

    def __init__(self, spectral_blocks, aspect_ratio, length=None):

        if length is None:
            length = 2 * np.pi / aspect_ratio

        coordinate_system = PlanarCartesianCoordinateSystem(extent=((0., aspect_ratio * length / 2), (0., length)))
        SymbolicBasis.__init__(self, coordinate_system)
        self.substitutions.append((self._n, aspect_ratio))

        awavenum = channel_wavenumbers(spectral_blocks)

        for wv in awavenum:

            mode_eq = fourier_functions(wv, self._n, self.coordinate_system)

            if mode_eq is not None:
                self.functions.append(mode_eq)


class PlanarBasinFourierBasis(SymbolicBasis):
    """Fourier basis defined on a closed basin, with no-flux boundary conditions in both the zonal and meridional
    direction :math:`x` and :math:`y`.

    Parameters
    ----------
    spectral_blocks: ~numpy.ndarray(int)
        Spectral blocks detailing the modes :math:`x`- and :math:`y`-wavenumber.
        Array of shape (nblocks, 2), where `nblocks` is the number of spectral blocks.
    aspect_ratio: float
        Spatial domain aspect ratio, :math:`n = 2 L_y/L_x` .
    length: float, optional
        Length of the domain along the :math:`x` coordinate: :math:`L_x` .
    """

    _n = symbols('n', positive=True)

    def __init__(self, spectral_blocks, aspect_ratio, length=None):

        if length is None:
            length = 2 * np.pi / aspect_ratio

        coordinate_system = PlanarCartesianCoordinateSystem(extent=((0., aspect_ratio * length / 2), (0., length)))
        SymbolicBasis.__init__(self, coordinate_system)
        self.substitutions.append((self._n, aspect_ratio))

        owavenum = basin_wavenumbers(spectral_blocks)

        for wv in owavenum:

            mode_eq = fourier_functions(wv, self._n, self.coordinate_system)

            if mode_eq is not None:
                self.functions.append(mode_eq)


def contiguous_basin_basis(nxmax, nymax, aspect_ratio, length=None):
    """Function that returns the basis for contiguous spectral blocks of modes on a closed basin.

    Parameters
    ----------
    nxmax: int
        Maximum x-wavenumber to fill the spectral block up to.
    nymax: int
        Maximum :math:`y`-wavenumber to fill the spectral block up to.
    aspect_ratio: float
        Spatial domain aspect ratio, :math:`n = 2 L_y/L_x` .
    length: float, optional
        Length of the domain along the :math:`x` coordinate: :math:`L_x` .

    Returns
    -------
    PlanarBasinFourierBasis
        The closed basin contiguous basis up to the specified spectral truncation.
    """

    spectral_blocks = np.zeros((nxmax * nymax, 2), dtype=int)
    i = 0
    for nx in range(1, nxmax + 1):
        for ny in range(1, nymax+1):
            spectral_blocks[i, 0] = nx
            spectral_blocks[i, 1] = ny
            i += 1

    return PlanarBasinFourierBasis(spectral_blocks, aspect_ratio, length)


def contiguous_channel_basis(nxmax, nymax, aspect_ratio, length=None):
    """Function that returns the basis for contiguous spectral blocks of modes on a channel.

    Parameters
    ----------
    nxmax: int
        Maximum x-wavenumber to fill the spectral block up to.
    nymax: int
        Maximum :math:`y`-wavenumber to fill the spectral block up to.
    aspect_ratio: float
        Spatial domain aspect ratio, :math:`n = 2 L_y/L_x` .
    length: float, optional
        Length of the domain along the :math:`x` coordinate: :math:`L_x` .

    Returns
    -------
    ChannelFourierBasis
        The channel contiguous basis up to the specified spectral truncation.
    """

    spectral_blocks = np.zeros((nxmax * nymax, 2), dtype=int)
    i = 0
    for nx in range(1, nxmax + 1):
        for ny in range(1, nymax+1):
            spectral_blocks[i, 0] = nx
            spectral_blocks[i, 1] = ny
            i += 1

    return PlanarChannelFourierBasis(spectral_blocks, aspect_ratio, length)


def fourier_functions(wave_number, n, coordinate_system):
    """Function that return Fourier modes expressions:

    * `'A'` for a function of the form :math:`F^A_{P} (x, y) =  \sqrt{2}\, \cos(P y) = \sqrt{2}\, \cos(n_y\, y)`
    * `'K'` for a function of the form :math:`F^K_{M,P} (x, y) =  2\cos(M nx)\, \sin(P y) = 2\cos(n_x\,  n\, x)\, \sin(n_y\, y)`
    * `'L'` for a function of the form :math:`F^L_{H,P} (x, y) = 2\sin(H nx)\, \sin(P y) = 2\sin(n_x\, n \,x)\, \sin(n_y\, y)`

    Parameters
    ----------
    wave_number: WaveNumber
        The wavenumber and type information of the mode to be returned.

    Returns
    -------
    `Sympy`_ expression
        Symbolic expression of the mode.

    .. _Sympy: https://www.sympy.org/

    """
    mode_eq = None
    x, y = coordinate_system.coordinates_symbol
    if wave_number.type == 'A':
        mode_eq = sqrt(2) * cos(wave_number.ny * y)
    elif wave_number.type == 'K':
        mode_eq = 2 * cos(wave_number.nx * n * x) * sin(wave_number.ny * y)
    elif wave_number.type == 'L':
        mode_eq = 2 * sin(wave_number.nx * n * x) * sin(wave_number.ny * y)

    return mode_eq


class WaveNumber(object):
    """Class to define model base functions wavenumber. The basis function available are:

    * `'A'` for a function of the form :math:`F^A_{P} (x, y) =  \sqrt{2}\, \cos(P y) = \sqrt{2}\, \cos(n_y\, y)`
    * `'K'` for a function of the form :math:`F^K_{M,P} (x, y) =  2\cos(M nx)\, \sin(P y) = 2\cos(n_x\,  n\, x)\, \sin(n_y\, y)`
    * `'L'` for a function of the form :math:`F^L_{H,P} (x, y) = 2\sin(H nx)\, \sin(P y) = 2\sin(n_x\, n \,x)\, \sin(n_y\, y)`

    where :math:`x` and :math:`y` are the nondimensional model's domain coordinates (see :ref:`files/model/oro_model:Projecting the equations on a set of basis functions`).

    Parameters
    ----------
    function_type: str
        One character string to define the type of basis function. It can be `'A'`, `'K'` or `'L'`.
    P: int
        The :math:`y` wavenumber integer.
    M: int
        The :math:`x` wavenumber integer.
    H: int
        The :math:`x` wavenumber integer.
    nx: float
        The :math:`x` wavenumber.
    ny: float
        The :math:`y` wavenumber.

    Attributes
    ----------
    type: str
        One character string to define the type of basis function. It can be `'A'`, `'K'` or `'L'`.
    P: int
        The :math:`y` wavenumber integer.
    M: int
        The :math:`x` wavenumber integer.
    H: int
        The :math:`x` wavenumber integer.
    nx: float
        The :math:`x` wavenumber.
    ny: float
        The :math:`y` wavenumber.

    """

    def __init__(self, function_type, P, M, H, nx, ny):
        self.type = function_type
        self.P = P
        self.M = M
        self.H = H
        self.nx = nx
        self.ny = ny

    def __repr__(self):
        return "type = {}, P = {}, M= {},H={}, nx= {}, ny={}".format(self.type, self.P, self.M, self.H, self.nx, self.ny)


def channel_wavenumbers(spectral_blocks):
    """Functions that returns the :class:`WaveNumber` objects corresponding to a given list of spectral blocks for a
    channel-like spatial domain.

    Parameters
    ----------
    spectral_blocks: ~numpy.ndarray(int)
        Spectral blocks detailing the modes :math:`x`- and :math:`y`-wavenumber.
        Array of shape (nblocks, 2), where `nblocks` is the number of spectral blocks.

    Returns
    -------
    ~numpy.ndarray(WaveNumber)
        The array of wavenumber objects corresponding to the given spectral blocks.
    """
    # initialization of the variables
    wave_numbers = list()

    # Atmospheric wavenumbers definition

    for i in range(spectral_blocks.shape[0]):  # function type is limited to AKL for the moment: atmosphere is a channel

        if spectral_blocks[i, 0] == 1:
            wave_numbers.append(WaveNumber('A', spectral_blocks[i, 1], 0, 0, 0, spectral_blocks[i, 1]))
            wave_numbers.append(WaveNumber('K', spectral_blocks[i, 1], spectral_blocks[i, 0], 0, spectral_blocks[i, 0], spectral_blocks[i, 1]))
            wave_numbers.append(WaveNumber('L', spectral_blocks[i, 1], 0, spectral_blocks[i, 0], spectral_blocks[i, 0], spectral_blocks[i, 1]))
        else:
            wave_numbers.append(WaveNumber('K', spectral_blocks[i, 1], spectral_blocks[i, 0], 0, spectral_blocks[i, 0], spectral_blocks[i, 1]))
            wave_numbers.append(WaveNumber('L', spectral_blocks[i, 1], 0, spectral_blocks[i, 0], spectral_blocks[i, 0], spectral_blocks[i, 1]))

    return np.array(wave_numbers)


def basin_wavenumbers(spectral_blocks):
    """Functions that returns the :class:`WaveNumber` objects corresponding to a given list of spectral blocks for a
    closed basin spatial domain.

    Parameters
    ----------
    spectral_blocks: ~numpy.ndarray(int)
        Spectral blocks detailing the modes :math:`x`- and :math:`y`-wavenumber.
        Array of shape (nblocks, 2), where `nblocks` is the number of spectral blocks.

    Returns
    -------
    ~numpy.ndarray(WaveNumber)
        The array of wavenumber objects corresponding to the given spectral blocks.
    """

    # initialization of the variables
    wave_numbers = list()

    # Oceanic wavenumbers definition

    for i in range(spectral_blocks.shape[0]):  # function type is limited to L for the moment: ocean is a closed basin
        wave_numbers.append(WaveNumber('L', spectral_blocks[i, 1], 0, spectral_blocks[i, 0], spectral_blocks[i, 0] / 2., spectral_blocks[i, 1]))

    return np.array(wave_numbers)
