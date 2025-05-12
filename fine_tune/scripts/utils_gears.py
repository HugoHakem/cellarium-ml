import os
import shutil
import sys
from pathlib import Path
from zipfile import ZipFile

import requests
from tqdm import tqdm

# Code Credits: https://github.com/snap-stanford/GEARS/blob/master/gears/utils.py

def print_sys(s: str):
    """
    system pring in error output
xR
    Args:
        s (str): the string to print
    """
    print(s, flush=True, file=sys.stderr)

def dataverse_download(url: str, save_path: Path):
    """
    Dataverse download helper with progress bar

    Args:
        url (str): the url of the dataset
        path (str): the path to save the dataset
    """

    if os.path.exists(save_path):
        print_sys('Found local copy...')
    else:
        print_sys("Downloading...")
        response = requests.get(url, stream=True)
        total_size_in_bytes= int(response.headers.get('content-length', 0))
        progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)
        with open(save_path, 'wb') as file:
            for data in response.iter_content(chunk_size=1024):
                progress_bar.update(len(data))
                file.write(data)
        progress_bar.close()

def zip_data_download_wrapper(url: str, save_path: Path, data_path: Path):
    """
    Wrapper for zip file download.
    Guarantee the zip extracted file is flatten and name after save_path.stem

    Args:
        url (str): the url of the dataset
        save_path (str): the path where the file is donwloaded
        data_path (str): the path to save the extracted dataset
    """

    if os.path.exists(data_path.joinpath(save_path.stem)):
        print_sys('Found local copy...')
    else:
        dataverse_download(url, save_path.with_suffix('.zip'))
        print_sys('Extracting zip file...')

        # Make tmp dir to extract zip
        tmp_dir = data_path / f"__temp_extract_{save_path.stem}"
        tmp_dir.mkdir(parents=True, exist_ok=True)

        with ZipFile(save_path.with_suffix('.zip'), 'r') as zip:
            zip.extractall(path = tmp_dir)

        # Flatten path
        flattened_path = flatten_single_child_dirs(tmp_dir)

        # Final output directory
        output_dir = data_path / save_path.stem
        output_dir.mkdir(parents=True, exist_ok=True)

        # Move contents of the final subfolder to output_dir
        for item in flattened_path.iterdir():
            shutil.move(item, output_dir)

        # Clean up temp_dir
        shutil.rmtree(tmp_dir)

        print_sys("Done!")

def flatten_single_child_dirs(start_path: Path) -> Path:
    """
    Given a directory, recursively descend while there's only one subdirectory and no files.
    Returns the deepest such directory.
    """
    current = start_path
    while True:
        children = list(current.iterdir())
        dirs = [d for d in children if d.is_dir()]
        files = [f for f in children if f.is_file()]
        if len(dirs) == 1 and not files:
            current = dirs[0]
        else:
            break
    return current
