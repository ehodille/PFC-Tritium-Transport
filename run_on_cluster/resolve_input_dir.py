"""Resolve an input folder by searching inside PFC-TT and one level above it.

When given a bare folder name (no path separators), this module checks both
the repository root directory and its parent directory.  If the folder exists
in both locations (and they are distinct), execution is aborted with a clear
error message.  Absolute and relative paths with separators are used directly.
"""

import os
import sys


def resolve_input_dir(name: str, repo_root: str = None) -> str:
    """Return the absolute path to the input folder *name*.

    Parameters
    ----------
    name : str
        Folder argument provided by the user.  May be:
        * An absolute path (``/home/…/my_folder``)
        * A relative path containing separators (``../my_folder``)
        * A bare folder name (``my_folder``) – triggers the dual-location
          search described below.
    repo_root : str, optional
        Root of the PFC-Tritium-Transport repository.  Defaults to ``CWD``.

    Dual-location search (bare folder name)
    ----------------------------------------
    1. ``<repo_root>/<name>``  – "inside PFC-TT"
    2. ``<repo_root>/../<name>`` – "outside PFC-TT" (one level up)

    If found in **both** locations the function exits with an error.
    """
    if repo_root is None:
        repo_root = os.getcwd()

    # ---- absolute path → use directly ----
    if os.path.isabs(name):
        if os.path.isdir(name):
            return name
        print(f"Error: Input directory '{name}' not found!")
        sys.exit(1)

    # ---- relative path with separators → resolve from CWD ----
    if os.sep in name or "/" in name:
        if os.path.isdir(name):
            return os.path.abspath(name)
        print(f"Error: Input directory '{name}' not found!")
        sys.exit(1)

    # ---- bare folder name → search repo_root AND parent ----
    inside = os.path.join(repo_root, name)
    outside = os.path.join(repo_root, "..", name)

    found_inside = os.path.isdir(inside)
    found_outside = os.path.isdir(outside)

    inside_abs = os.path.abspath(inside) if found_inside else None
    outside_abs = os.path.abspath(outside) if found_outside else None

    # Guard against both resolving to the same physical directory
    if found_inside and found_outside and inside_abs != outside_abs:
        print(f"Error: Input directory '{name}' found in BOTH locations:")
        print(f"  Inside PFC-TT:  {inside_abs}")
        print(f"  Outside PFC-TT: {outside_abs}")
        print(f"  Provide a full or relative path to disambiguate.")
        sys.exit(1)

    if found_inside:
        print(f"  Resolved input directory: {inside_abs}")
        return inside_abs
    if found_outside:
        print(f"  Resolved input directory (outside PFC-TT): {outside_abs}")
        return outside_abs

    print(f"Error: Input directory '{name}' not found!")
    print(f"  Searched: {os.path.abspath(inside)}")
    print(f"  Searched: {os.path.abspath(outside)}")
    sys.exit(1)
