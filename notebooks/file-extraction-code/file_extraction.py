'''READ ME 
- This code will auto extract your .fits files from their zip folder. Make sure all your .ZIP files are stored in a folder on their own.
- Ensure you change source_folder and destination_folder to the correct paths where your zip files are located and where you want the extracted files to be saved.
- If you're running this in the OCL you'll might need to pip install the below libraries. I wrote this for VS code on my local machine - it might need some tweaks to work in Jupyter Notebooks
'''
from pathlib import Path
import zipfile
import shutil

# Set your folder paths to where your zip files are located, make sure you use FORWARD slashes and quotation marks!
source_folder = Path(r"D:/Documents/UNIVERSITY/Year 2 2026-2027/SXPS288_Remote_Experiments/04_raw_calibration")

#this is where you want the unzipped files to be saved.
destination_folder = Path(r"D:/Documents/UNIVERSITY/Year 2 2026-2027/SXPS288_Remote_Experiments/04_raw_calibration")

destination_folder.mkdir(parents=True, exist_ok=True)

zip_files = list(source_folder.glob("*.zip"))
extracted = 0

for zip_path in zip_files:
    with zipfile.ZipFile(zip_path, "r") as archive:
        fits_files = [
            f for f in archive.namelist()
            if f.lower().endswith(".fits") and not f.endswith("/")
        ]

        for fits_file in fits_files:
            filename = Path(fits_file).name
            output_path = destination_folder / filename

            # Avoid overwriting files with identical names
            counter = 1
            while output_path.exists():
                output_path = destination_folder / (
                    f"{Path(filename).stem}_{counter}.fits"
                )
                counter += 1

            with archive.open(fits_file) as src:
                with open(output_path, "wb") as dst:
                    shutil.copyfileobj(src, dst)

            extracted += 1
            print(f"Extracted: {output_path.name}")

print(f"\nCompleted: {extracted} FITS files extracted! Now go to your folder they should be there.")

