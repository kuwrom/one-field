"""The instruments.

Each probe is a standalone experiment: run(report=...) -> dict, with
its verdict asserted internally, runnable directly as
python3 <probe>.py.  Checking whether the theory is flawless is
testing: the probes cover the mechanisms, edge cases, and rabbit
holes, and the pytest freezes in tests/ pin their verdicts.
"""
from . import _bootstrap  # noqa: F401
