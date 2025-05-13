import os
import shutil
import sys
import tarfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from zipfile import ZipFile

import requests
from tqdm import tqdm

# Adapted from: https://github.com/snap-stanford/GEARS/blob/master/gears/utils.py

def print_sys(s: str):
    """
    system pring in error output
xR
    Args:
        s (str): the string to print
    """
    print(s, flush=True, file=sys.stderr)

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

def downloader(url: str, save_path: Path):
    """
    download helper with progress bar

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

def data_download_wrapper(url: str, save_path: Path, data_path: Path, ext: str=".zip"):
    """
    Wrapper for zip file download.
    Guarantee the zip extracted file is flatten and name after save_path.stem

    Args:
        url (str): the url of the dataset
        save_path (Path): the path where the file is donwloaded
        data_path (Path): the path to save the extracted dataset
        ext (str): the extension file to expect. Should start with a '.'
    """

    if os.path.exists(data_path.joinpath(save_path.stem)):
        print_sys('Found local copy...')
    else:
        downloader(url, save_path.with_suffix(ext))
        print_sys(f'Extracting {ext} file...')

        # Make tmp dir to extract zip
        tmp_dir = data_path / f"__temp_extract_{save_path.stem}"
        tmp_dir.mkdir(parents=True, exist_ok=True)

        if ext == ".zip":
            with ZipFile(save_path.with_suffix('.zip'), 'r') as zip:
                zip.extractall(path = tmp_dir)
        if ext == ".tar":
            with tarfile.open(save_path.with_suffix('.tar'), 'r') as tar:
                tar.extractall(path=tmp_dir)
        else:
            raise ValueError(
                "Only .tar file and .zar are supported file extensions." \
                f"You provided ext={ext}" \
                "make sure it starts with a '.'"
            )

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

# === SECTION: multiprocessing download ===

def download_chunk(url, start, end, idx, tmp_dir):
    headers = {'Range': f'bytes={start}-{end}'}
    response = requests.get(url, headers=headers, stream=True)
    chunk_path = tmp_dir / f"part_{idx}.tmp"
    with open(chunk_path, 'wb') as f:
        for chunk in response.iter_content(1024):
            f.write(chunk)

def parallel_download(url: str, save_path: Path, num_threads: int = 4):
    response = requests.head(url)
    file_size = int(response.headers.get('content-length', 0))

    if 'accept-ranges' not in response.headers:
        raise RuntimeError("Server does not support partial content downloading (no 'accept-ranges' support)")

    tmp_dir = save_path.parent / f"__tmp_{save_path.stem}"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    chunk_size = file_size // num_threads
    futures = []

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        for i in range(num_threads):
            start = i * chunk_size
            end = (start + chunk_size - 1) if i != num_threads - 1 else file_size - 1
            futures.append(executor.submit(download_chunk, url, start, end, i, tmp_dir))

        # Wait for all chunks to complete
        for f in tqdm(futures, desc="Downloading chunks"):
            f.result()

    # Merge chunks
    with open(save_path, 'wb') as outfile:
        for i in range(num_threads):
            chunk_path = tmp_dir / f"part_{i}.tmp"
            with open(chunk_path, 'rb') as infile:
                outfile.write(infile.read())
            chunk_path.unlink()

    tmp_dir.rmdir()
    print(f"Saved to {save_path}")

