
from pathlib import Path
from astropy.io import fits

# Main directory containing your FITS images
base_folder = Path(
    r"D:/Documents/UNIVERSITY/Year 2 2026-2027/SXPS288_Remote_Experiments/03_raw_science_frames/03_raw_science_frames_ALL - Copy/M15"
)

# True = preview changes without renaming
# False = actually rename the files
dry_run = False

renamed = 0
skipped = 0

# Find all FITS files, including those in subfolders
for file_path in base_folder.rglob("*"):

    if not file_path.is_file():
        continue

    if file_path.suffix.lower() not in [".fits", ".fit", ".fts"]:
        continue

    try:
        # Read the FITS header
        header = fits.getheader(file_path)

        # Extract exposure time
        exposure = header.get(
            "EXPTIME",
            header.get("EXPOSURE")
        )

        if exposure is None:
            print(f"No exposure time: {file_path.name}")
            skipped += 1
            continue

        # Convert exposure time to a number
        exposure = float(exposure)

        # Format exposure as 16, 32, 64, etc.
        # Decimals are retained when necessary
        exposure_text = f"{exposure:g}"

        prefix = f"{exposure_text}_"

        # Avoid adding the same prefix twice
        if file_path.name.startswith(prefix):
            print(f"Already renamed: {file_path.name}")
            skipped += 1
            continue

        # Construct new filename
        new_filename = prefix + file_path.name

        new_path = file_path.with_name(new_filename)

        # Prevent overwriting existing files
        if new_path.exists():
            print(f"File already exists: {new_filename}")
            skipped += 1
            continue

        # Preview or rename
        if dry_run:
            print(
                f"PREVIEW: {file_path.name} "
                f"-> {new_filename}"
            )
        else:
            file_path.rename(new_path)

            print(
                f"RENAMED: {file_path.name} "
                f"-> {new_filename}"
            )

        renamed += 1

    except Exception as e:
        print(f"ERROR: {file_path.name}: {e}")
        skipped += 1

print("\n-----------------------------")
print(f"Files {'identified' if dry_run else 'renamed'}: {renamed}")
print(f"Files skipped: {skipped}")
print("-----------------------------")
