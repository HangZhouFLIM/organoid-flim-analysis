"""Reconstruct descriptive summaries and native-bin predictions from archived fits."""
import numpy as np
import pandas as pd
from scipy.stats import t

from workflows.reconvolution import response


CELL_METRICS = ["tau1_ns", "alpha1_percent", "tau2_ns", "alpha2_percent",
                "amplitude_weighted_lifetime_ns", "intensity_weighted_lifetime_ns"]
BINS = [f"bin_{i:03d}" for i in range(129)]


def summary(values, keys):
    stats = values.groupby(keys, observed=True)["value"].agg(["count", "mean", "std", "sem"])
    stats = stats.rename(columns={"std": "sd"})
    difference = t.ppf(0.975, stats["count"] - 1) * stats["sem"]
    stats["ci95_low"] = stats["mean"] - difference
    stats["ci95_high"] = stats["mean"] + difference
    return stats.drop(columns="count")


def compare_summary(calculated, folder, output, name, keys):
    saved = pd.read_csv(folder / name).set_index(keys).sort_index()
    calculated = calculated.sort_index()
    if not saved.index.is_unique or not saved.index.equals(calculated.index):
        raise ValueError(f"Group or metric mismatch in {name}.")
    for col in ["unit", "n_images", "n_cells", "mean", "sd", "sem", "ci95_low", "ci95_high"]:
        if col not in saved:
            continue
        if col == "unit":
            agrees = calculated[col].equals(saved[col])
        else:
            agrees = np.allclose(calculated[col], saved[col], rtol=1e-10, atol=1e-12, equal_nan=True)
        if not agrees:
            raise ValueError(f"Reconstructed {col} differs from {name}.")
    calculated.to_csv(output / name)


def cell_summaries(fits, folder, output, panel, secondary=False):
    keys = ["condition"]
    prefix = ""
    if secondary:
        splits = pd.read_csv(folder / "secondary_intensity_groups.csv")
        joined = splits.merge(fits, on=["cell_id", "image_id", "condition"], validate="one_to_one")
        if len(joined) != len(splits):
            raise ValueError("An intensity-group assignment has no corresponding fit.")
        fits = joined
        keys.append("intensity_group")
        prefix = "secondary_"
    rows = fits.melt(id_vars=keys + ["image_id"], value_vars=CELL_METRICS,
                     var_name="metric", value_name="value")
    keys = keys + ["metric"]
    counts = rows.groupby(keys).agg(n_cells=("value", "count"), n_images=("image_id", "nunique"))
    for level in (["cell"] if panel == "S3" else ["cell", "image"]):
        obs = rows
        if level == "image":
            obs = rows.groupby(keys + ["image_id"], as_index=False)["value"].mean()
        stats = summary(obs, keys).join(counts).reset_index()
        stats["unit"] = np.where(stats.metric.str.endswith("percent"), "%", "ns")
        compare_summary(stats.set_index(keys), folder, output, f"{prefix}{level}_summary.csv", keys)


def roi_summaries(fits, folder, output):
    keys = ["display_panel", "figure_panel", "model_components", "lifetime_component", "condition"]
    rows = fits.melt(id_vars=keys + ["sample_id"], value_vars=["lifetime_ns", "alpha_percent"],
                     var_name="metric", value_name="value")
    rows["metric"] = rows.metric.map({"lifetime_ns": "lifetime", "alpha_percent": "alpha"})
    keys = keys + ["metric"]
    calculated = summary(rows, keys)
    calculated["n_images"] = rows.groupby(keys).sample_id.nunique()
    calculated["unit"] = np.where(calculated.index.get_level_values("metric") == "lifetime", "ns", "%")
    compare_summary(calculated, folder, output, "image_summary.csv", keys)


def replay(fits, folder, output, panel, curve_id=None):
    """Evaluate the recorded coefficients and compare all 129 bins with the saved predictions."""
    roi = panel == "S14"
    field = "sample_id" if roi else "image_id"
    endpoint = "sample_id" if roi else "cell_id"
    index = [endpoint, "model_components"] if roi else [endpoint]
    settings = pd.read_csv(folder / "acquisition_settings.csv").set_index(field, verify_integrity=True)
    decay = pd.read_csv(folder / "native_decays.csv").set_index(endpoint, verify_integrity=True)
    saved = pd.read_csv(folder / "predicted_decays.csv.gz").set_index(index, verify_integrity=True)
    irfs = {}
    for key, group in pd.read_csv(folder / "irf_samples.csv").groupby(field):
        weights = group.normalized_weight.to_numpy()
        if not np.isfinite(weights).all() or (weights < 0).any() or not np.isclose(weights.sum(), 1):
            raise ValueError(f"Invalid normalized IRF for {key}.")
        positive = weights > 0
        irfs[key] = (group.time_ns.to_numpy()[positive], weights[positive])
    by_fit = fits.set_index(index).sort_index()
    if not roi and not by_fit.index.is_unique:
        raise ValueError(f"Duplicated fitted cell in {panel}.")
    if set(fits[field]) != set(irfs) or set(irfs) != set(settings.index):
        raise ValueError("Fitted acquisitions, settings and IRFs do not agree.")
    if set(decay.index) != set(fits[endpoint]) or set(saved.index) != set(by_fit.index):
        raise ValueError("The fitted, observed and predicted record identifiers do not agree.")
    measured = decay[BINS].to_numpy(float)
    if not np.isfinite(measured).all() or (measured < 0).any() or (measured != np.floor(measured)).any():
        raise ValueError("Native decay bins must contain nonnegative integer counts.")
    if not (settings.native_bins.eq(129) & settings.fit_start_bin_zero_based.eq(0) &
            settings.fit_end_bin_exclusive.eq(128) & settings.overflow_bin_zero_based.eq(128)).all():
        raise ValueError("Unexpected acquisition window in the released fits.")
    if curve_id is not None and curve_id not in decay.index:
        raise ValueError(f"Unknown {endpoint}: {curve_id}")
    reports = []
    for idx, expected in zip(saved.index, saved[BINS].to_numpy()):
        if roi:
            specimen, k = idx
            group = by_fit.loc[[idx]].sort_values("lifetime_component")
            if group.lifetime_component.tolist() != list(range(1, k + 1)):
                raise ValueError(f"Unexpected component numbers for {idx}.")
            row = group.iloc[0]
            amplitudes = group.component_amplitude.to_numpy()
            lifetimes = group.lifetime_ns.to_numpy()
            image = specimen
        else:
            specimen, k = idx, 2
            row = by_fit.loc[idx]
            image = row.image_id
            if decay.at[specimen, "image_id"] != image or decay.at[specimen, "source_label_id"] != row.source_label_id:
                raise ValueError(f"Native decay and fitted cell assignments differ for {specimen}.")
            amplitudes = [row.amplitude1_counts_per_ns, row.amplitude2_counts_per_ns]
            lifetimes = [row.tau1_ns, row.tau2_ns]
        timing = settings.loc[image]
        prediction = response(amplitudes, lifetimes, row.background_counts_per_bin, row.irf_shift_ns,
                              *irfs[image], timing.bin_width_ns, timing.pulse_period_ns, len(BINS))
        if not np.allclose(prediction, expected, rtol=1e-8, atol=1e-8):
            raise ValueError(f"Replayed counts differ from the saved {k}-component curve for {specimen}.")
        reports.append({endpoint: specimen, "model_components": k,
                        "max_absolute_count_difference": np.max(abs(prediction - expected))})
        if specimen == curve_id:
            pd.DataFrame({"bin": range(129), "included_in_fit": np.arange(129) < 128,
                          "observed_counts": decay.loc[specimen, BINS].to_numpy(),
                          "replayed_counts": prediction, "archived_counts": expected}).to_csv(
                              output / f"selected_curve_{k}_components.csv", index=False)
    pd.DataFrame(reports).to_csv(output / "curve_comparisons.csv", index=False)
    return len(reports)


def reconstruct(panel, root, output, curve_id=None):
    folder = root / "fitting" / panel
    fits = pd.read_csv(folder / "fit_parameters.csv")
    if panel == "S14":
        if len(fits) != 365 or fits.sample_id.nunique() != 73:
            raise ValueError("S14 requires 365 component rows from 73 image-specific ROIs.")
        roi_summaries(fits, folder, output)
    else:
        expected = {"S3": (1715, 12), "S7": (2298, 27)}[panel]
        if (len(fits), fits.image_id.nunique()) != expected:
            raise ValueError(f"Unexpected number of cells or images in {panel}.")
        cell_summaries(fits, folder, output, panel)
        if panel == "S7":
            cell_summaries(fits, folder, output, panel, secondary=True)
    return replay(fits, folder, output, panel, curve_id)
