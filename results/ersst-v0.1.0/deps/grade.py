# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Grading rules, frozen in PREREGISTRATION.md (2026-11-24) before the job runs."""

DD23_TC = 2065.0          # corrected (2025) MLE tipping year for SG - 2xGM, HadISST 1870-2020
DD23_CI95 = (2037.0, 2109.0)


def grade_recon(r, rms):
    """Our HadISST rebuild of DD23's SG-2xGM fingerprint vs their published file, 1870-2020."""
    if r is None or rms is None:
        return "RECON_FAIL"
    if r >= 0.98 and rms <= 0.05:
        return "RECON_OK"
    if r >= 0.95:
        return "RECON_WEAK"
    return "RECON_FAIL"


def grade_control(tc):
    """Our estimator on DD23's own published series must give their corrected 2065."""
    if tc is None:
        return "REPRO_FAIL"
    d = abs(tc - DD23_TC)
    if d <= 5.0:
        return "REPRO_PASS"
    if d <= 10.0:
        return "REPRO_NEAR"
    return "REPRO_FAIL"


def grade_forecast(q):
    """q = share of model paths whose 2021-latest mean lies below the observed mean."""
    if q is None:
        return "NOT_RUN"
    if q < 0.05:
        return "BELOW_FORECAST"
    if q > 0.95:
        return "ABOVE_FORECAST"
    return "WITHIN_FORECAST"


def grade_shift(tc_2020, tc_latest, unbounded=False):
    """Change of the point tipping year when 2021-latest data are added."""
    if unbounded:
        return "UNBOUNDED"
    if tc_2020 is None or tc_latest is None:
        return "NOT_RUN"
    d = tc_latest - tc_2020
    if d <= -5.0:
        return "EARLIER"
    if d >= 5.0:
        return "LATER"
    return "STABLE"


def grade_onestep(summary, crit=2.576):
    """Flags for the out-of-sample one-step residuals (expected mean 0, variance 1)."""
    flags = []
    if not summary or summary.get("n", 0) < 3:
        return ["NOT_RUN"]
    if summary["t_mean"] > crit:
        flags.append("BIAS_HIGH")
    elif summary["t_mean"] < -crit:
        flags.append("BIAS_LOW")
    if summary["z_var"] > crit:
        flags.append("VARIANCE_HIGH")
    elif summary["z_var"] < -crit:
        flags.append("VARIANCE_LOW")
    if summary.get("ks_p") is not None and summary["ks_p"] < 0.01:
        flags.append("KS_REJECT")
    return flags or ["CALIBRATED"]
