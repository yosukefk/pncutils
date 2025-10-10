#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unfortunately i am not sure if this is working: more debugging needed

ioapi_sanitize.py — validate & (optionally) fix IOAPI-like metadata.

Features
- Works on any "netCDF-like" object (PseudoNetCDF Dataset or netCDF4 Dataset) in-memory:
    sanitize_inplace(nc, include_statics=False): dict report
- Command line copy-cleaner (never modifies input file):
    python ioapi_sanitize.py in.nc out.nc [--include-statics] [--dry-run]
- Rules enforced:
    • Data vars = 4D (TSTEP, LAY, ROW, COL), excluding TFLAG
    • VAR-LIST = concat of 16-char padded data var names
    • NVARS (attr) == #data vars
    • VAR dimension length == NVARS
    • DATE-TIME dimension exists and == 2
    • TFLAG dims = (TSTEP, VAR, DATE-TIME) (created/reshaped if needed)
"""

from __future__ import annotations
import sys
import argparse
from typing import List, Dict, Any
import numpy as np

try:
    from netCDF4 import Dataset  # only needed for CLI path
except Exception:
    Dataset = None  # library still usable without netCDF4 installed

IOAPI_NAME_LEN = 16
TFLAG_NAME = "TFLAG"
TFLAG_NAMES = ('TFLAG', 'ETFLAG')
REQ_DIMS = ("TSTEP", "LAY", "ROW", "COL")

def _pad16(name: str) -> str:
    return name[:IOAPI_NAME_LEN].ljust(IOAPI_NAME_LEN)

def _chunk16(s: str) -> List[str]:
    return [s[i:i+IOAPI_NAME_LEN].strip() for i in range(0, len(s), IOAPI_NAME_LEN)]

def _is_unlimited(nc, dname: str) -> bool:
    dim = nc.dimensions[dname]
    # PNC dims often have .isunlimited; netCDF4 dims have method isunlimited()
    return bool(getattr(dim, "isunlimited", lambda: getattr(dim, "isunlimited", False))())

def _is_data_var(nc, vname: str) -> bool:
    #if vname == TFLAG_NAME:
    if vname in (TFLAG_NAMES):
        return False
    var = nc.variables[vname]
    dims = getattr(var, "dimensions", ())
    return len(dims) >= 4 and dims[:4] == REQ_DIMS

def _list_data_vars(nc, include_statics: bool = False) -> List[str]:
    data = [v for v in nc.variables if _is_data_var(nc, v)]
    if include_statics:
        statics = [v for v in nc.variables
                   #if v not in data and v != TFLAG_NAME
                   if v not in data and v not in  TFLAG_NAMES
                   and getattr(nc.variables[v], "dimensions", ()) == ("ROW", "COL")]
        data.extend(statics)
    return data

def _build_varlist(vnames: List[str]) -> str:
    return "".join(_pad16(n) for n in vnames)

def _get_len(nc, dname: str) -> int | None:
    if dname not in nc.dimensions:
        return None
    dim = nc.dimensions[dname]
    try:
        return len(dim)
    except TypeError:
        return None

def validate(nc, include_statics: bool = False) -> Dict[str, Any]:
    """
    Return a report (dict) describing IOAPI invariants and problems, without modifying nc.
    """
    issues: List[str] = []
    report: Dict[str, Any] = {"issues": issues}

    # required dims present?
    for d in REQ_DIMS:
        if d not in nc.dimensions:
            issues.append(f"missing-dimension:{d}")

    ntime = _get_len(nc, "TSTEP")
    nlay  = _get_len(nc, "LAY")
    nrow  = _get_len(nc, "ROW")
    ncol  = _get_len(nc, "COL")

    report["shape"] = {"TSTEP": ntime, "LAY": nlay, "ROW": nrow, "COL": ncol}

    data_vars = _list_data_vars(nc, include_statics=include_statics)
    report["data_vars"] = data_vars
    nvars_should = len(data_vars)

    # NVARS attribute
    nvars_attr = getattr(nc, "NVARS", None)
    if nvars_attr is None:
        issues.append("missing-attr:NVARS")
    elif nvars_attr != nvars_should:
        issues.append(f"bad-attr:NVARS={nvars_attr} != {nvars_should}")

    # VAR dimension
    var_dim_len = _get_len(nc, "VAR")
    if var_dim_len is None:
        issues.append("missing-dimension:VAR")
    elif var_dim_len != nvars_should:
        issues.append(f"bad-dim:VAR={var_dim_len} != {nvars_should}")

    # VAR-LIST attribute
    expected_varlist = _build_varlist(data_vars)
    varlist_attr = getattr(nc, "VAR-LIST", None)
    if varlist_attr is None:
        issues.append("missing-attr:VAR-LIST")
    else:
        if len(varlist_attr) != len(expected_varlist):
            issues.append(f"bad-attr:VAR-LIST-length={len(varlist_attr)} expected {len(expected_varlist)}")
        else:
            if _chunk16(varlist_attr) != data_vars:
                issues.append("bad-attr:VAR-LIST-names-mismatch")

    # DATE-TIME dimension
    dt = _get_len(nc, "DATE-TIME")
    if dt is None:
        issues.append("missing-dimension:DATE-TIME")
    elif dt != 2:
        issues.append(f"bad-dim:DATE-TIME={dt} != 2")

    # TFLAG variable
    if TFLAG_NAME not in nc.variables:
        issues.append("missing-var:TFLAG")
    else:
        tf = nc.variables[TFLAG_NAME]
        dims = getattr(tf, "dimensions", ())
        if dims[:3] != ("TSTEP", "VAR", "DATE-TIME"):
            issues.append(f"bad-var:TFLAG-dims={dims}")

    return report

def sanitize_inplace(nc, include_statics: bool = False) -> Dict[str, Any]:
    """
    Modify 'nc' in-place to satisfy IOAPI invariants. Returns a report dict.
    Works for PseudoNetCDF or netCDF4-like objects that support createDimension/Variable and setattr on globals.
    """
    report = validate(nc, include_statics=include_statics)
    data_vars = report["data_vars"]
    nvars = len(data_vars)

    # Ensure DATE-TIME dimension = 2
    if "DATE-TIME" not in nc.dimensions:
        nc.createDimension("DATE-TIME", 2)
    elif _get_len(nc, "DATE-TIME") != 2:
        # Some backends can’t resize dims; rebuild would be required.
        # Best-effort: delete+recreate if allowed.
        try:
            del nc.dimensions["DATE-TIME"]  # may fail depending on backend
        except Exception:
            pass
        try:
            nc.createDimension("DATE-TIME", 2)
        except Exception:
            # leave as-is; report already flagged
            pass

    # Ensure VAR dimension size == NVARS
    if "VAR" not in nc.dimensions:
        nc.createDimension("VAR", nvars)
    else:
        if _get_len(nc, "VAR") != nvars:
            # Same caveat as above for resizing
            try:
                del nc.dimensions["VAR"]
            except Exception:
                pass
            try:
                nc.createDimension("VAR", nvars)
            except Exception:
                pass

    # Update NVARS attribute
    try:
        setattr(nc, "NVARS", nvars)
    except Exception:
        pass

    # Update VAR-LIST
    try:
        setattr(nc, "VAR-LIST", _build_varlist(data_vars))
    except Exception:
        pass

    # Ensure TFLAG exists and has correct dims
    if TFLAG_NAME not in nc.variables or getattr(nc.variables[TFLAG_NAME], "dimensions", ())[:3] != ("TSTEP","VAR","DATE-TIME"):
        # Create (or recreate) TFLAG as int32
        try:
            if TFLAG_NAME in nc.variables:
                try:
                    del nc.variables[TFLAG_NAME]
                except Exception:
                    pass
            tf = nc.createVariable(TFLAG_NAME, "i4", ("TSTEP","VAR","DATE-TIME"))
            tf.long_name = "Timestep start date and time"
            tf.units = "YYYYDDD,HHMMSS"
            tf.var_desc = "Timestep start date and time"
            tf[:] = 0
        except Exception:
            # Backend prevents create? leave as-is; report already flagged
            pass

    # Re-validate after fixes
    return validate(nc, include_statics=include_statics)

# ----------------- CLI path: copy-cleaner using netCDF4 -----------------

def _copy_clean(infile: str, outfile: str, include_statics: bool = False) -> Dict[str, Any]:
    if Dataset is None:
        raise RuntimeError("netCDF4 is required for CLI copy-cleaning; pip install netCDF4")

    with Dataset(infile, "r") as src, Dataset(outfile, "w") as dst:
        # Copy all dims except VAR and DATE-TIME (we enforce those below)
        for dname, dim in src.dimensions.items():
            if dname in ("VAR", "DATE-TIME"):
                continue
            dst.createDimension(dname, None if dim.isunlimited() else len(dim))

        # Determine data variables and enforce dims
        data_vars = _list_data_vars(src, include_statics=include_statics)
        nvars_new = len(data_vars)
        dst.createDimension("VAR", nvars_new)
        dst.createDimension("DATE-TIME", 2)

        # Copy global attrs except NVARS/VAR-LIST (we set after)
        for a in src.ncattrs():
            if a in ("NVARS", "VAR-LIST"):
                continue
            setattr(dst, a, getattr(src, a))
        setattr(dst, "NVARS", nvars_new)
        setattr(dst, "VAR-LIST", _build_varlist(data_vars))

        # Copy variables (except TFLAG; handle after to preserve with logic)
        for vname, v in src.variables.items():
            if vname in TFLAG_NAMES:
                continue
            fill = getattr(v, "_FillValue", None)
            kwargs = {}
            if fill is not None:
                kwargs["fill_value"] = fill
            outv = dst.createVariable(vname, v.dtype, v.dimensions, **kwargs)
            # copy attrs
            for a in v.ncattrs():
                if a == "_FillValue":
                    continue
                setattr(outv, a, getattr(v, a))
            outv[:] = v[:]

        # TFLAG: preserve if consistent; else raise; else zeros if absent
        _clone_or_build_tflag(src, dst, nvars_new)

        return validate(dst, include_statics=include_statics)

def _copy_clean_old(infile: str, outfile: str, include_statics: bool = False) -> Dict[str, Any]:
    if Dataset is None:
        raise RuntimeError("netCDF4 is required for CLI copy-cleaning; pip install netCDF4")

    with Dataset(infile, "r") as src, Dataset(outfile, "w") as dst:
        # Copy all dims except VAR and DATE-TIME (we enforce those below)
        for dname, dim in src.dimensions.items():
            if dname in ("VAR", "DATE-TIME"):
                continue
            dst.createDimension(dname, None if dim.isunlimited() else len(dim))
        # enforce dims
        data_vars = _list_data_vars(src, include_statics=include_statics)
        dst.createDimension("VAR", len(data_vars))
        dst.createDimension("DATE-TIME", 2)

        # Copy global attrs except NVARS/VAR-LIST (we set after)
        for a in src.ncattrs():
            if a in ("NVARS", "VAR-LIST"):
                continue
            setattr(dst, a, getattr(src, a))
        setattr(dst, "NVARS", len(data_vars))
        setattr(dst, "VAR-LIST", _build_varlist(data_vars))

        # Copy variables (except TFLAG; handle after)
        for vname, v in src.variables.items():
            #if vname == TFLAG_NAME:
            if vname in TFLAG_NAMES:
                continue
            fill = getattr(v, "_FillValue", None)
            kwargs = {}
            if fill is not None:
                kwargs["fill_value"] = fill
            outv = dst.createVariable(vname, v.dtype, v.dimensions, **kwargs)
            # copy attrs
            for a in v.ncattrs():
                if a == "_FillValue":
                    continue
                setattr(outv, a, getattr(v, a))
            outv[:] = v[:]

        # Ensure TFLAG
        ntime = _get_len(dst, "TSTEP")
        tf = dst.createVariable(TFLAG_NAME, "i4", ("TSTEP","VAR","DATE-TIME"))
        tf.long_name = "Timestep start date and time"
        tf.units = "YYYYDDD,HHMMSS"
        tf.var_desc = "Timestep start date and time"
        tf[:] = 0

        # Return a report built from the destination
        return validate(dst, include_statics=include_statics)

def _clone_or_build_tflag(src, dst, nvars_new: int):
    """
    Create dst TFLAG with dims (TSTEP, VAR, DATE-TIME).
    If src has a TFLAG with consistent slices across VAR, copy/tile it.
    Else if src has inconsistent per-var TFLAG, raise ValueError.
    Else (no TFLAG): create zeros.
    """
    # Ensure dims exist on dst
    if "DATE-TIME" not in dst.dimensions:
        dst.createDimension("DATE-TIME", 2)
    if "VAR" not in dst.dimensions:
        dst.createDimension("VAR", nvars_new)

    tf_dst = dst.createVariable("TFLAG", "i4", ("TSTEP", "VAR", "DATE-TIME"))
    tf_dst.long_name = "Timestep start date and time"
    tf_dst.units = "YYYYDDD,HHMMSS"
    tf_dst.var_desc = "Timestep start date and time"

    if "TFLAG" not in src.variables:
        # No source: zeros
        tf_dst[:] = 0
        return

    tf_src = src.variables["TFLAG"]
    dims_src = getattr(tf_src, "dimensions", ())
    # We only try to “preserve” when source is (TSTEP, VAR, DATE-TIME) or at least 3D in that order
    if dims_src[:3] != ("TSTEP", "VAR", "DATE-TIME"):
        # malformed -> treat as not-preservable; zero-out (or you could raise if you prefer)
        tf_dst[:] = 0
        return

    arr = np.array(tf_src[:], copy=False)
    if arr.ndim < 3 or arr.shape[2] != 2:
        tf_dst[:] = 0
        return

    # Check per-var slice consistency: all equal to the first var’s slice
    # Shape: (TSTEP, VAR_old, 2)
    first = arr[:, 0, :]
    equal_all = np.all(arr == first[:, None, :])  # broadcast compare across VAR axis
    if not equal_all:
        raise ValueError("Inconsistent TFLAG across VAR; refusing to guess. (All per-var TFLAG must match)")

    # Tile to new VAR length
    tiled = np.tile(first[:, None, :], (1, nvars_new, 1))
    tf_dst[:] = tiled

def main():
    ap = argparse.ArgumentParser(description="Validate/fix IOAPI-like NetCDF (VAR-LIST/NVARS/TFLAG).")
    ap.add_argument("infile")
    ap.add_argument("outfile", nargs="?", help="Output path (omit with --dry-run)")
    ap.add_argument("--include-statics", action="store_true", help="Include 2-D (ROW,COL) statics in VAR-LIST")
    ap.add_argument("--dry-run", action="store_true", help="Only print report; do not write output")
    args = ap.parse_args()

    if args.dry_run:
        if Dataset is None:
            print("netCDF4 not installed; cannot open file in CLI mode.", file=sys.stderr)
            sys.exit(2)
        with Dataset(args.infile, "r") as src:
            rep = validate(src, include_statics=args.include_statics)
        print("[REPORT] data_vars:", rep["data_vars"])
        print("[REPORT] issues:", "OK" if not rep["issues"] else rep["issues"])
        sys.exit(0)

    if not args.outfile:
        print("ERROR: outfile is required unless --dry-run is used", file=sys.stderr)
        sys.exit(2)

    rep = _copy_clean(args.infile, args.outfile, include_statics=args.include_statics)
    print("[DONE] Wrote:", args.outfile)
    print("[REPORT] data_vars:", rep["data_vars"])
    print("[REPORT] issues:", "OK" if not rep["issues"] else rep["issues"])

if __name__ == "__main__" and "pytest" not in sys.modules:
    main()

