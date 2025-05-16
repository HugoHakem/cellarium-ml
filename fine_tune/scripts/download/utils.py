import sys
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union
from urllib.parse import unquote_to_bytes

import requests
from requests.structures import CaseInsensitiveDict
from tqdm import tqdm
from tqdm.contrib.concurrent import thread_map


@dataclass
class FileEntry:
    url: str
    out_path: Path
    raw_path: Optional[Path] = None
    filename: Optional[str] = None

    def __post_init__(self):
        if not isinstance(self.out_path, Path):
            if isinstance(self.out_path, str):
                self.out_path = Path(self.out_path)
            else:
                raise ValueError(
                    f"`out_path` must be a Path or str, got {type(self.out_path)}"
                )

        if self.raw_path is not None and not isinstance(self.raw_path, Path):
            if isinstance(self.raw_path, str):
                self.raw_path = Path(self.raw_path)
            else:
                raise ValueError(
                    f"`raw_path` must be a Path or str if provided, got {type(self.raw_path)}"
                )

FileEntryGroup = dict[str, FileEntry]
FileMap = dict[str, Union[FileEntry, FileEntryGroup]]


# === SECTION: utilities ===

def print_sys(s: str):
    """
    Print a message to standard error with flush.

    Args:
        s (str): The string to print.
    """
    print(s, flush=True, file=sys.stderr)

def flatten_single_child_dirs(start_path: Path) -> Path:
    """
    Recursively descend into a directory as long as there is only one subdirectory and no files.
    Returns the deepest such directory.

    Args:
        start_path (Path): The starting directory path.

    Returns:
        Path: The deepest directory with only one subdirectory and no files.
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

def get_file_name(headers: CaseInsensitiveDict) -> str:
    """
    Extract the filename from HTTP response headers.

    Args:
        headers (CaseInsensitiveDict): The HTTP response headers.

    Returns:
        str: The filename extracted from the Content-Disposition header.
    """
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

# === SECTION: base and parallel downloader ===

def downloader(url: str, save_path: Path):
    """
    Download a file from a URL with a progress bar.

    Args:
        url (str): The URL of the dataset.
        save_path (Path): The path to save the downloaded file.
    """

    print_sys(f"Downloading {save_path.name}...")
    response = requests.get(url, stream=True)
    total_size_in_bytes= int(response.headers.get('content-length', 0))
    progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    with open(save_path, 'wb') as file:
        for data in response.iter_content(chunk_size=1024 * 1024):
            progress_bar.update(len(data))
            file.write(data)
    progress_bar.close()


def download_chunk(url, start, end, idx, tmp_dir):
    """
    Download a chunk of a file for parallel downloading.

    Args:
        url (str): The URL of the file.
        start (int): The start byte of the chunk.
        end (int): The end byte of the chunk.
        idx (int): The chunk index.
        tmp_dir (Path): The temporary directory to save the chunk.
    """
    headers = {'Range': f'bytes={start}-{end}'}
    response = requests.get(url, headers=headers, stream=True)
    chunk_path = tmp_dir / f"part_{idx}.tmp"
    with open(chunk_path, 'wb') as f:
        for chunk in response.iter_content(1024 * 1024):
            f.write(chunk)

def parallel_downloader(url: str, save_path: Path, max_workers: int|None=None):
    """
    Download a file in parallel using multiple threads.

    Args:
        url (str): The URL of the file.
        save_path (Path): The path to save the downloaded file.
        max_workers (int, optional): The number of parallel workers. If None, uses the default.
    """
    print_sys(f"Downloading in parallel {save_path.name}...")
    response = requests.head(url)
    file_size = int(response.headers.get('content-length', 0))

    tmp_dir = save_path.parent / f"__tmp_{save_path.stem}"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # If no max_workers provided, get ThreadPoolExecutor default
    if not max_workers:
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            max_workers = ex._max_workers

    chunk_size = file_size // max_workers
    ranges = [
        (url, i * chunk_size, (file_size - 1 if i == max_workers - 1 else (i + 1) * chunk_size - 1), i, tmp_dir)
        for i in range(max_workers)
    ]

    # Download in parallel using tqdm's thread_map
    thread_map(lambda args: download_chunk(*args), ranges, max_workers=max_workers, desc="Downloading chunks")

    # Merge chunks
    with open(save_path, 'wb') as outfile:
        for i in range(max_workers):
            chunk_path = tmp_dir / f"part_{i}.tmp"
            with open(chunk_path, 'rb') as infile:
                outfile.write(infile.read())
            chunk_path.unlink()

    tmp_dir.rmdir()
