from dataclasses import dataclass, field
from pathlib import Path
from typing import Union

from jsonargparse import CLI
from jsonargparse.typing import register_type

from fine_tune.scripts.download.core import iter_download_url
from fine_tune.scripts.download.utils import FileEntry, FileMap


# === SECTION: Custom serializer and deserializer to register custom data types ===
def _serialize_path(val):
    return str(val) if isinstance(val, Path) else val

def file_entry_serializer(val: FileEntry) -> dict:
    if isinstance(val, FileEntry):
        return {key: _serialize_path(val.__getattribute__(key))
                 for key in ["url", "out_path", "raw_path", "filename"]}
    raise ValueError(f"Cannot parse FileEntry from: {val}")

def file_entry_deserializer(val: Union[FileEntry, dict]) -> FileEntry:
    if isinstance(val, FileEntry):
        return val
    if isinstance(val, dict) and all(k in val for k in ('url', 'out_path')):
        return FileEntry(**val)
    raise ValueError(f"Cannot parse FileEntry from: {val}")

def file_map_serializer(val: FileMap) -> dict:
    if not isinstance(val, dict):
        raise ValueError("Expected a dict for FileMap")
    parsed = {}
    for key, v in val.items():
        # Try single FileEntry
        if isinstance(v, FileEntry):
            parsed[key] = file_entry_serializer(v)
        # Try FileEntryGroup
        elif isinstance(v, dict):
            group = {}
            for subkey, subv in v.items():
                if isinstance(subv, FileEntry):
                    group[subkey] = file_entry_serializer(subv)
                else:
                    raise ValueError(f"Expected a FileEntry for '{key}': {v}")
            parsed[key] = group
        else:
            raise ValueError(f"Expected a FileEntry for '{key}': {v}")
    return parsed

def file_map_deserializer(val: dict) -> FileMap:
    if not isinstance(val, dict):
        raise ValueError("Expected a dict for FileMap")

    parsed = {}
    for key, v in val.items():
        # Accept a FileEntry directly
        if isinstance(v, FileEntry):
            parsed[key] = v
        # Accept a dict that can be turned into a FileEntry
        elif isinstance(v, dict) and all(k in v for k in ('url', 'out_path')):
            parsed[key] = FileEntry(**v)
        # Try FileEntryGroup
        elif isinstance(v, dict):
            group = {}
            for subkey, subv in v.items():
                if isinstance(subv, FileEntry):
                    group[subkey] = subv
                elif isinstance(subv, dict) and all(k in subv for k in ('url', 'out_path')):
                    group[subkey] = FileEntry(**subv)
                else:
                    raise ValueError(f"Invalid value for subkey '{subkey}': {subv}")
            parsed[key] = group
        else:
            raise ValueError(f"Invalid value for key '{key}': {v}")
    return parsed


register_type(FileEntry, serializer=file_entry_serializer, deserializer=file_entry_deserializer)
register_type(FileMap, serializer=file_map_serializer, deserializer=file_map_deserializer)


# === SECTION: Core of the CLI ===

@dataclass
class DownloadConfig:
    file_map: FileMap = field( # Ensure that each instance of a dataclass has its own default value
        default_factory=lambda: {
            "example_file": FileEntry(
                url="https://example.com/file.csv",
                out_path=Path("./data/file.csv"),
                raw_path=Path("./data/.raw/file.csv"),
                filename="file.csv"
            ),
            "example_group": {
                "subfile1": FileEntry(
                    url="https://example.com/file1.csv",
                    out_path=Path("./data/group/file1.csv")
                ),
                "subfile2": FileEntry(
                    url="https://example.com/file2.csv",
                    out_path=Path("./data/group/file2.csv")
                )
            }
        }
    )
    use_key_name: bool = False

def cli_iter_download_url(cfg: DownloadConfig):
    iter_download_url(file_map=cfg.file_map, use_key_name=cfg.use_key_name)


if __name__ == "__main__":
    CLI(cli_iter_download_url)
