"""Resolve an input folder by searching simulations/, inside PFC-TT, and one level above it.

When given a bare folder name (no path separators), this module checks three
locations in priority order:
  1. ``<repo_root>/simulations/<name>``  – recommended location for working inputs
  2. ``<repo_root>/<name>``             – directly inside PFC-TT
  3. ``<repo_root>/../<name>``          – parent directory of PFC-TT

If the bare name resolves to more than one distinct location, execution is
aborted with a clear error message.  Absolute and relative paths with
separators are used directly without any search.
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
        * A bare folder name (``my_folder``) – triggers the three-location
          search described below.
    repo_root : str, optional
        Root of the PFC-Tritium-Transport repository.  Defaults to ``CWD``.

    Three-location search (bare folder name, in priority order)
    -----------------------------------------------------------
    1. ``<repo_root>/simulations/<name>``  – recommended location for working inputs
    2. ``<repo_root>/<name>``             – directly inside PFC-TT
    3. ``<repo_root>/../<name>``          – parent directory of PFC-TT

    If found in more than one distinct location the function exits with an error.
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

    # ---- bare folder name → search repo_root/simulations/, repo_root/, and parent ----
    in_simulations = os.path.join(repo_root, "simulations", name)
    inside = os.path.join(repo_root, name)
    outside = os.path.join(repo_root, "..", name)

    found_simulations = os.path.isdir(in_simulations)
    found_inside = os.path.isdir(inside)
    found_outside = os.path.isdir(outside)

    simulations_abs = os.path.abspath(in_simulations) if found_simulations else None
    inside_abs = os.path.abspath(inside) if found_inside else None
    outside_abs = os.path.abspath(outside) if found_outside else None

    # Collect distinct locations where the folder was found
    found_locs = {loc for loc in [simulations_abs, inside_abs, outside_abs] if loc and os.path.isdir(loc)}
    # Deduplicate by resolved path
    unique_locs = list(dict.fromkeys(found_locs))
    if len(unique_locs) > 1:
        print(f"Error: Input directory '{name}' found in multiple locations:")
        for loc in unique_locs:
            print(f"  {loc}")
        print(f"  Provide a full or relative path to disambiguate.")
        sys.exit(1)

    if found_simulations:
        print(f"  Resolved input directory (simulations/): {simulations_abs}")
        return simulations_abs
    if found_inside:
        print(f"  Resolved input directory: {inside_abs}")
        return inside_abs
    if found_outside:
        print(f"  Resolved input directory (outside PFC-TT): {outside_abs}")
        return outside_abs

    print(f"Error: Input directory '{name}' not found!")
    print(f"  Searched: {os.path.abspath(in_simulations)}")
    print(f"  Searched: {os.path.abspath(inside)}")
    print(f"  Searched: {os.path.abspath(outside)}")
    sys.exit(1)
