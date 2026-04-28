#Define version based on setup script
from importlib.metadata import version, PackageNotFoundError
try:
    __version__ = version("NucleoATAC2")
except PackageNotFoundError:
    __version__ = "0.0.0"
