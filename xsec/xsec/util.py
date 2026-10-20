"""Small helpers."""
import importlib.util
import os
import sys

MEASUREMENTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "measurements")


def load_measurement_module(measurement, module="analysis"):
    """Import measurements/<measurement>/<module>.py under a unique module name."""
    name = f"_mva_{measurement}_{module}"
    if name in sys.modules:
        return sys.modules[name]
    path = os.path.join(MEASUREMENTS, measurement, module + ".py")
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod
