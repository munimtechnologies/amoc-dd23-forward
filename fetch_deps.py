# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Fetch the v2 estimator modules, unchanged, from the repository commit that v2 ran with.

Pinned commit 74c83b30f81c4b32ee7ad03d1a8460ca05550aa9 of munimtechnologies/amoc-dd23-forward
(ddfit.py, ddrep.py, grade.py, hadisst.py).  The ERSSTv5 test uses exactly the estimator that
produced the HadISST1 v2 grades.
"""
import os
import sys
import time
import urllib.request

COMMIT = "74c83b30f81c4b32ee7ad03d1a8460ca05550aa9"
BASE = "https://raw.githubusercontent.com/munimtechnologies/amoc-dd23-forward/%s/" % COMMIT
FILES = ("ddfit.py", "ddrep.py", "grade.py", "hadisst.py")
MUST = {"ddfit.py": "def forecast_check", "ddrep.py": "def fit_dd23", "grade.py": "def grade_shift",
        "hadisst.py": "def monthly_anomalies"}


def main(dest=None):
    here = dest or os.path.dirname(os.path.abspath(__file__))
    os.makedirs(here, exist_ok=True)
    for name in FILES:
        err = None
        for k in range(5):
            try:
                req = urllib.request.Request(BASE + name, headers={"User-Agent": "mv2-research-job/1.0"})
                with urllib.request.urlopen(req, timeout=120) as r:
                    data = r.read()
                break
            except Exception as e:      # retry transient network errors
                err = e
                time.sleep(5 * (k + 1))
        else:
            raise RuntimeError("could not fetch %s: %r" % (name, err))
        text = data.decode("utf-8")
        if MUST[name] not in text:
            raise RuntimeError("%s does not look like the v2 module" % name)
        with open(os.path.join(here, name), "w", encoding="utf-8") as f:
            f.write(text)
        print("fetched %s (%d bytes) from commit %s" % (name, len(data), COMMIT[:10]), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
