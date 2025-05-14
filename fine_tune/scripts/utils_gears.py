import shutil
import sys
import tarfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import unquote_to_bytes
from zipfile import ZipFile

import requests
from requests.structures import CaseInsensitiveDict
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

def get_file_name(headers: CaseInsensitiveDict):
    content_disposition = {
        key: val
        for field in headers["Content-Disposition"].replace(" ", "").split(';')
        if (
            (parts := field.split("=", 1)) and
            (key_star := parts[0]) and
            (key := key_star.replace("*", "")) and
            (val_encoded := parts[1].replace('"', '') if len(parts) > 1 else parts[0]) and
            (val_tuple := val_encoded.partition("''")) and
            (val := unquote_to_bytes(val_tuple[2]).decode(val_tuple[0]) if key_star.endswith("*") else val_tuple[0])
           )
    }
    return content_disposition["filename"]

def downloader(url: str, save_path: Path):
    """
    download helper with progress bar

    Args:
        url (str): the url of the dataset
        path (str): the path to save the dataset
    """

    print_sys(f"Downloading {save_path.name}...")
    response = requests.get(url, stream=True)
    total_size_in_bytes= int(response.headers.get('content-length', 0))
    progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)
    with open(save_path, 'wb') as file:
        for data in response.iter_content(chunk_size=1024):
            progress_bar.update(len(data))
            file.write(data)
    progress_bar.close()

def data_download_wrapper(
        url: str,
        out_path: Path,
        raw_path: Path|None=None,
        filename: str|None=None,
        max_workers: int|None=None
        ):
    """
    Wrapper for file download.
    Guarantee the zip/tar extracted file is flatten and name after raw_path.stem
    Supports single / multiple file download, but those must already be unziped.

    Args:
        url (str or list): the `url` of the dataset
        raw_path (Path): the path where the file is donwloaded
        data_path (Path): the path to save the extracted dataset
        ext (str or None): the extension file to expect. Should start with a '.'
            Must be in [".zip", ".tar", None]
        filename: (str or list or None): the desired filename.
            If `url` is a list and filename not `None`, filename must be a list with same length
            If `None`, will be inferred.
    """
    headers = requests.head(url).headers
    if not filename:
        try:
            filename = get_file_name(headers)
        except Exception:
            raise ValueError(
                "`filename` is not given and couldn't be retrieved from the `url`.\n" \
                "Please provide a `filename` with valid file extension at the end"
            )
    save_path = raw_path / filename if raw_path else out_path / filename
    # Check if already there else download
    if save_path.exists():
        print_sys(f"[{save_path.name}] Found local copy in {save_path.parent}")
    else:
        faster_download = 'Accept-Ranges' in headers
        if faster_download:
            parallel_download(
                url=url,
                save_path=save_path,
                max_workers=max_workers
            )
        else:
            downloader(
                url=url,
                save_path=save_path
                )

    if (ext := save_path.suffix) in ['.zip', '.tar']:
        print_sys(f'Extracting {save_path.name} file...')
        extr_path = out_path / save_path.stem

        # Check if already ther, else extract
        if extr_path.is_dir():
            if not any(extr_path.iterdir()):
                print_sys(f"[{extr_path.name}] Found extracted file in {extr_path.parent}")
        else:
            # Make tmp dir to extract zip
            tmp_dir = out_path / f"__temp_extract_{extr_path.stem}"
            tmp_dir.mkdir(parents=True, exist_ok=True)

            if ext == ".zip":
                with ZipFile(save_path.with_suffix('.zip'), 'r') as zip:
                    zip.extractall(path = tmp_dir)
            elif ext == ".tar":
                with tarfile.open(save_path.with_suffix('.tar'), 'r') as tar:
                    tar.extractall(path=tmp_dir)

            # Flatten extracted file and move content in extr_path
            flattened_path = flatten_single_child_dirs(tmp_dir)
            extr_path.mkdir(parents=True, exist_ok=True)
            for item in flattened_path.iterdir():
                shutil.move(item, extr_path)

            # Clean up tmp_dir
            shutil.rmtree(tmp_dir)
    elif len(ext := save_path.suffixes) > 1:
            raise Warning(
                "The downloaded file sounds like being compressed, but the decompression scheme" \
                "doesn't support this file extension.\n" \
                f"The file extension is: {ext}"
            )
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
    print_sys(f"Downloading in parallel {save_path.name}...")
    response = requests.head(url)
    file_size = int(response.headers.get('content-length', 0))

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

