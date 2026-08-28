"""Dataset providers.

One class per data source, each wrapping whatever that provider actually
offers — a REST API, a Google Drive folder, an open-data HTTP server, a cloud
Zarr catalogue — behind a common ``fetch`` → ``load`` → ``extract`` lifecycle.
:class:`~climdata.utils.wrapper_workflow.ClimateExtractor` drives them from a
Hydra configuration; they can also be used directly.

The lifecycle has per-provider deviations worth knowing, documented on each
class. Two recur:

* **Extract before load.** :class:`~climdata.datasets.MSWX.MSWXmirror`,
  :class:`~climdata.datasets.NEXGDDP.NEXGDDP` and
  :class:`~climdata.datasets.HOSTRADA.HOSTRADAmirror` apply the region of
  interest *while* opening each file, so calling ``extract()`` first is what
  keeps a multi-year request tractable.
* **Not a grid.** :class:`~climdata.datasets.DWD.DWDmirror` returns station
  measurements and :class:`~climdata.datasets.NASAPOWER.POWER` a single point,
  so neither has spatial axes to map or regrid.

Optional dependencies are imported defensively, so importing this package never
requires the extras of a provider you are not using; the error naming the
missing one is raised when the class is constructed.
"""

from .AGRI_ISIMIP import AGRI_ISIMIP
from .CMIPCloud import CMIPCloud
from .CMIP_W5E5 import CMIPW5E5
from .DWD import DWDmirror
from .ERA5 import ERA5Mirror
from .HOSTRADA import HOSTRADAmirror
from .HYRAS import HYRASmirror
from .MSWX import MSWXmirror
from .NASAPOWER import POWER
from .NEXGDDP import NEXGDDP
from .W5E5 import W5E5

__all__ = [
    "AGRI_ISIMIP",
    "CMIPCloud",
    "CMIPW5E5",
    "DWDmirror",
    "ERA5Mirror",
    "HOSTRADAmirror",
    "HYRASmirror",
    "MSWXmirror",
    "NEXGDDP",
    "POWER",
    "W5E5",
]
