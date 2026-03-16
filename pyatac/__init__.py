
#Define version based on setup script
from importlib.metadata import version, PackageNotFoundError
try:
    __version__ = version("NucleoATAC")
except PackageNotFoundError:
    __version__ = "0.0.0"
