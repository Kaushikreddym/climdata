"""Statistical downscaling and bias adjustment (SDBA).

Bias correction and statistical downscaling of climate model output, following
the ISIMIP3BASD method: :class:`BiasCorrection` adjusts a simulation towards an
observational reference on a coarse grid, :class:`StatisticalDownscaling` maps
the corrected coarse field onto a fine grid, and :class:`BCSD` chains the two.

Example:
    >>> from climdata.sdba import BCSD                     # doctest: +SKIP
    >>> out = BCSD(obs_fine, sim_hist, sim_fut).run()      # doctest: +SKIP
"""

from .bcsd import BCSD, BiasCorrection, StatisticalDownscaling, regrid_to_coarse
from .utils import compute_daily_climatology, smooth_fft

__all__ = [
    "BCSD",
    "BiasCorrection",
    "StatisticalDownscaling",
    "regrid_to_coarse",
    "compute_daily_climatology",
    "smooth_fft",
]
