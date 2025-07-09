import shutil
import tarfile
from pathlib import Path
from typing import Optional
from zipfile import ZipFile

import requests  # type: ignore

from .utils import (
    FileEntry,
    FileMap,
    downloader,
    flatten_single_child_dirs,
    get_file_name,
    parallel_downloader,
    print_sys,
)

# === SECTION: core downloader and iterator ===

def download_url(
        url: str,
        out_path: Path,
        raw_path: Optional[Path]=None,
        filename: Optional[str]=None,
        max_workers: Optional[int]=None,
    ):
    """
    Wrapper for file download and extraction.
    Ensures the zip/tar extracted file is flattened and named after raw_path.stem.
    Supports single or multiple file downloads, but those must already be unzipped.

    Args:
        url (str): The URL of the dataset.
        out_path (Path): The path to save the extracted dataset.
        raw_path (Path, optional): The path where the file is downloaded.
        filename (str, optional): The desired filename. If None, will be inferred.
        max_workers (int, optional): Number of parallel workers for download.
    """
    headers = requests.head(url).headers
    # Format filename
    inferred_filename = None
    try:
        inferred_filename = get_file_name(headers)
    except Exception:
        pass

    ## Use provided filename or fall back to inferred one
    filename = filename or inferred_filename
    if not filename:
        raise ValueError(
            "`filename` is not given and couldn't be retrieved from the `url`.\n"
            "Please provide a `filename` with a valid file extension."
        )

    ## Ensure file extension are presents
    filename_path = Path(filename)
    if not filename_path.suffixes:
        if not inferred_filename:
            raise ValueError(
                "`filename` has no file extension and it / they couldn't be infered from the `url`.\n" \
                "Please provide a `filename` with valid file extension at the end"
            )
        filename += "".join(Path(inferred_filename).suffixes)

    save_path = raw_path / filename if raw_path else out_path / filename

    # Check if already there else download
    if save_path.exists():
        print_sys(f"[{save_path.name}] Found local copy in {save_path.parent}")
    else:
        faster_download = False#'Accept-Ranges' in headers
        if faster_download:
            parallel_downloader(
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
            if any(extr_path.iterdir()):
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
    elif len(ext := save_path.suffixes) > 1: # type: ignore
        Warning(
            "The downloaded file sounds like being compressed, but the decompression scheme" \
            "doesn't support this file extension.\n" \
            f"The file extension is: {ext}"
        )
    print_sys("Done!")

def iter_download_url(file_map: FileMap, use_key_name: bool=False):
    """
    Iterate over a file map and download each file or group of files.

    This function processes a mapping of dataset names to FileEntry objects or
    dictionaries of FileEntry objects. For each entry, it calls `download_url`
    to download the file(s) to the specified locations. If `use_key_name` is True,
    the key name is used as the filename for the download. Note that if file.filename
    is specified, it will always prevail.

    Args:
        file_map (FileMap): A mapping of dataset names to FileEntry or dicts of FileEntry.
        use_key_name (bool, optional): If True, use the key name as the filename.
    """
    for _root, _value in file_map.items():
        if isinstance(file := _value, FileEntry):
            download_url(
                url=file.url,
                out_path=file.out_path,
                raw_path=file.raw_path,
                filename=file.filename if file.filename is not None or not use_key_name else _root,
                max_workers=None,
            )
        elif isinstance(folder := _value, dict):
            assert all(isinstance(v, FileEntry) for v in folder.values())

            for _name, file in folder.items():
                download_url(
                    url=file.url,
                    out_path=file.out_path / _root,
                    raw_path=file.raw_path / _root if file.raw_path else None,
                    filename = file.filename if file.filename is not None or not use_key_name else _name,
                    max_workers=None,
                )
