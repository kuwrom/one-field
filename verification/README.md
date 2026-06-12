# Verification reports

Independent replication is the most valuable contribution this project
can receive. To add a report:

1. Clone fresh, run:

   ```
   pip install -r requirements.txt
   pytest -q
   python interference/run.py
   python tools/scorecard.py --check
   ```

2. Add a file `reports/YYYY-MM-DD-<handle>.md` with: OS, Python and
   numpy versions, CPU architecture, the pytest summary line, and any
   numeric deviations at full precision (there should be none beyond
   float reproducibility).

Reports of failures or platform-dependent deviations are more valuable
than reports of success. File them as issues as well.
