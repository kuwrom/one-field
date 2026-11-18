"""The instruments.

Each probe is a standalone experiment: run(report=...) -> dict,
runnable directly as python3 <probe>.py.  Its verdict is asserted
inside the probe or frozen by the pytest files in tests/.  The probes
cover the mechanisms, the edge cases and the excluded readings.
"""
from . import _bootstrap  # noqa: F401
