"""Shortlist M15/M39 FITS frames using a known reference star.

Install: pip install astropy numpy scipy
First run with DRY_RUN = True, review CSV, then set False.
"""

from pathlib import Path
import csv
import shutil

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter
from scipy.optimize import least_squares

# ------------------------------------------------------
# SETTINGS
# ------------------------------------------------------
BASE = Path(
    r"D:/Documents/UNIVERSITY/Year 2 2026-2027/SXPS288_Remote_Experiments/03_raw_science_frames/03_raw_science_frames_ALL"
)
OUTPUT = BASE.parent / "SELECTED_REFERENCE_STAR_FRAMES"

TARGETS = {
    "M15": {
        "filename_contains": "PEGASUS_GLOB_00_M15",
        "xy": (870.63, 1943.97),
    },
    "M39": {
        "filename_contains": "M39_OPEN_00",
        "xy": (953.68, 1683.50),
    },
}

# "FITS": coordinates from AIJ's FITS X,Y display (1-based,
#         origin bottom-left).
# "IJ":   coordinates from AIJ's ImageJ X,Y display (0-based,
#         origin top-left). Use unflipped/unrotated display.
COORDINATE_SYSTEM = "IJ"

SEARCH_RADIUS = 12.0
MIN_FWHM = 2.5
MAX_FWHM = 4.5
MAX_PEAK = 60000.0

# Preview first: record measurements without copying FITS files.
DRY_RUN = True

# ------------------------------------------------------
# MEASUREMENT FUNCTIONS
# ------------------------------------------------------

def to_array_xy(x, y, image_height):
    """Convert AstroImageJ coordinates to zero-based FITS array indices."""
    if COORDINATE_SYSTEM.upper() == "FITS":
        return x - 1.0, y - 1.0
    if COORDINATE_SYSTEM.upper() == "IJ":
        return x, image_height - 1.0 - y
    raise ValueError("COORDINATE_SYSTEM must be 'FITS' or 'IJ'")


def from_array_xy(x, y, image_height):
    """Convert measured position back to the selected AIJ coordinate system."""
    if COORDINATE_SYSTEM.upper() == "FITS":
        return x + 1.0, y + 1.0
    return x, image_height - 1.0 - y


def measure_reference_star(image, expected_xy):
    """Locate the reference star and fit an elliptical 2D Gaussian.

    FWHM is the geometric mean of x and y Gaussian FWHMs.
    Peak is the highest *actual pixel value* within 3.5 pixels
    of the fitted star centre (not the fitted Gaussian peak).
    """
    height, width = image.shape
    x_expected, y_expected = to_array_xy(*expected_xy, height)

    if not (0 <= x_expected < width and 0 <= y_expected < height):
        raise ValueError("Expected star coordinates are outside the image")

    # Only look for a bright source within 10 pixels of the expected star.
    x1 = max(0, int(np.floor(x_expected - SEARCH_RADIUS - 1)))
    x2 = min(width, int(np.ceil(x_expected + SEARCH_RADIUS + 2)))
    y1 = max(0, int(np.floor(y_expected - SEARCH_RADIUS - 1)))
    y2 = min(height, int(np.ceil(y_expected + SEARCH_RADIUS + 2)))

    search = np.asarray(image[y1:y2, x1:x2], dtype=float)
    finite = np.isfinite(search)
    if not np.any(finite):
        raise ValueError("No valid pixels near expected reference star")

    filled = np.where(finite, search, np.median(search[finite]))
    smoothed = gaussian_filter(filled, sigma=1.0)
    sy_grid, sx_grid = np.mgrid[y1:y2, x1:x2]
    within = (sx_grid - x_expected)**2 + (sy_grid - y_expected)**2 <= SEARCH_RADIUS**2
    scores = np.where(within & finite, smoothed, -np.inf)
    if not np.isfinite(scores).any():
        raise ValueError("No valid pixels inside 10-pixel search radius")

    sy_local, sx_local = np.unravel_index(np.argmax(scores), scores.shape)
    cx, cy = x1 + int(sx_local), y1 + int(sy_local)

    # Fit a small cutout around that peak.
    fit_radius = 8
    xa, xb = max(0, cx - fit_radius), min(width, cx + fit_radius + 1)
    ya, yb = max(0, cy - fit_radius), min(height, cy + fit_radius + 1)
    cutout = np.asarray(image[ya:yb, xa:xb], dtype=float)
    valid = np.isfinite(cutout)
    if np.count_nonzero(valid) < 40:
        raise ValueError("Too few valid pixels for Gaussian fit")

    yy, xx = np.mgrid[ya:yb, xa:xb]
    edge = np.zeros(cutout.shape, dtype=bool)
    edge[0, :] = edge[-1, :] = True
    edge[:, 0] = edge[:, -1] = True
    background = float(np.median(cutout[edge & valid]))
    amplitude = float(np.max(cutout[valid]) - background)
    if amplitude <= 0:
        raise ValueError("No positive stellar signal detected")

    # Model: background + elliptical Gaussian.
    def model(p):
        a, x0, y0, sigma_x, sigma_y, bg = p
        return bg + a * np.exp(
            -0.5 * (((xx - x0) / sigma_x)**2 + ((yy - y0) / sigma_y)**2)
        )

    initial = [amplitude, cx, cy, 1.5, 1.5, background]
    lower = [0.0, max(xa, cx - 3), max(ya, cy - 3), 0.5, 0.5, -np.inf]
    upper = [np.inf, min(xb - 1, cx + 3), min(yb - 1, cy + 3), 6.0, 6.0, np.inf]

    fit = least_squares(
        lambda p: (model(p) - cutout)[valid],
        initial, bounds=(lower, upper), max_nfev=300,
        x_scale="jac",
    )
    if not fit.success:
        raise ValueError("Gaussian profile fit did not converge")

    _, x_fit, y_fit, sigma_x, sigma_y, _ = fit.x
    separation = float(np.hypot(x_fit - x_expected, y_fit - y_expected))

    # Calculate peak for the *measured star*, avoiding the wider field.
    star_core = ((xx - x_fit)**2 + (yy - y_fit)**2 <= 3.5**2) & valid
    if not star_core.any():
        raise ValueError("No valid pixels at fitted reference star")
    peak = float(np.max(cutout[star_core]))

    factor = 2.0 * np.sqrt(2.0 * np.log(2.0))
    fwhm = float(factor * np.sqrt(sigma_x * sigma_y))
    reported_x, reported_y = from_array_xy(x_fit, y_fit, height)

    return {
        "FWHM_pixels": round(fwhm, 3),
        "Peak_ADU": round(peak, 1),
        "Distance_pixels": round(separation, 3),
        "Star_X": round(float(reported_x), 3),
        "Star_Y": round(float(reported_y), 3),
    }


# ------------------------------------------------------
# SEARCH FILES, TEST LIMITS, AND COPY PASSES
# ------------------------------------------------------

def main():
    if not BASE.is_dir():
        raise FileNotFoundError(f"Science folder not found: {BASE}")

    OUTPUT.mkdir(parents=True, exist_ok=True)
    rows = []
    passed = 0
    copied = 0

    for path in BASE.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in {".fits", ".fit", ".fts"}:
            continue

        name_upper = path.name.upper()
        target = next(
            (name for name, config in TARGETS.items()
             if config["filename_contains"].upper() in name_upper),
            None,
        )
        if target is None:  # M52 and unrelated files are ignored.
            continue

        row = {"Cluster": target, "Filename": path.name, "Source": str(path)}
        try:
            with fits.open(path, memmap=False) as hdul:
                image_hdu = next(
                    (h for h in hdul if h.data is not None and h.data.ndim == 2),
                    None,
                )
                if image_hdu is None:
                    raise ValueError("No 2D FITS image found")
                measures = measure_reference_star(
                    image_hdu.data, TARGETS[target]["xy"]
                )
            row.update(measures)

            problems = []
            if not (MIN_FWHM <= measures["FWHM_pixels"] <= MAX_FWHM):
                problems.append("FWHM outside 2.5–4.5 px")
            if measures["Peak_ADU"] >= MAX_PEAK:
                problems.append("Peak >= 60000")
            if measures["Distance_pixels"] > SEARCH_RADIUS:
                problems.append("Star more than 10 px from expected position")

            if problems:
                row["Decision"] = "REVIEW"
                row["Reason"] = "; ".join(problems)
            else:
                row["Decision"] = "PASS"
                row["Reason"] = "All three limits passed"
                passed += 1
                if not DRY_RUN:
                    target_folder = OUTPUT / target
                    target_folder.mkdir(parents=True, exist_ok=True)
                    destination = target_folder / path.name
                    if destination.exists():
                        row["Reason"] += "; destination exists, not overwritten"
                    else:
                        shutil.copy2(path, destination)
                        copied += 1

        except Exception as exc:
            row["Decision"] = "REVIEW"
            row["Reason"] = f"Could not measure: {exc}"

        rows.append(row)
        print(f"{row['Decision']:6}  {target:3}  {path.name}  |  {row['Reason']}")

    csv_file = OUTPUT / "reference_star_results.csv"
    fieldnames = [
        "Cluster", "Filename", "Source", "FWHM_pixels", "Peak_ADU",
        "Distance_pixels", "Star_X", "Star_Y", "Decision", "Reason",
    ]
    with csv_file.open("w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nMatching FITS images examined: {len(rows)}")
    print(f"Images passing all criteria: {passed}")
    print(f"FITS files copied: {copied} {'(preview mode)' if DRY_RUN else ''}")
    print(f"Results CSV: {csv_file}")
    print("M52 was intentionally excluded.")


if __name__ == "__main__":
    main()
