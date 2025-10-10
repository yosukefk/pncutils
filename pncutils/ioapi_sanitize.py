#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ioapi_sanitize.py — validate & fix IOAPI-like NetCDF metadata.

Features
- In-memory (PNC/netCDF4-like object):
    sanitize_inplace(nc, include_statics=False) -> report dict
- CLI copy-cleaner (never modifies input file):
    python ioapi_sanitize.py in.nc out.nc [--include-statics] [--dry-run] [--strict]

Rules
- Data vars = 4D (TSTEP, LAY, ROW, COL) OR 3D (TSTEP, ROW, COL), excluding TFLAG
- VAR-LIST = concat of 16-char padded data var names (ordered as detected)
- NVARS (attr) == #data vars
- VAR dim length == NVARS   (copy-clean enforces; in-place refuses to resize)
- DATE-TIME dim exists and == 2
- TFLAG dims = (TSTEP, VAR, DATE-TIME)
    * If source TFLAG exists and all per-VAR slices are identical -> preserve by tiling to new VAR
    * If inconsistent across VAR -> raise ValueError
    * If missing/malformed -> create zeros
"""

from __future__ import annotations
from typing import List, Dict, Any
import sys
import argparse
import numpy as np

try:
    from netCDF4 import Dataset  # for CLI
except Exception:
    Dataset = None  # library functions still usable without netCDF4

# ----------------- constants -----------------

IOAPI_NAME_LEN = 16
TFLAG_NAME = "TFLAG"
TFLAG_NAMES = ('TFLAG', 'ETFLAG')
REQ_DIMS_4D = ("TSTEP", "LAY", "ROW", "COL")
REQ_DIMS_3D = ("TSTEP", "ROW", "COL")

# ----------------- tiny utils -----------------

def _pad16(name: str) -> str:
    return name[:IOAPI_NAME_LEN].ljust(IOAPI_NAME_LEN)

def _chunk16(s: str) -> List[str]:
    return [s[i:i+IOAPI_NAME_LEN].strip() for i in range(0, len(s), IOAPI_NAME_LEN)]

def _get_len(nc, dname: str) -> int | None:
    if dname not in nc.dimensions:
        return None
    dim = nc.dimensions[dname]
    try:
        return len(dim)
    except TypeError:
        return None

# ----------------- data-var detection (shared) -----------------

def _is_data_var(nc, vname: str) -> bool:
    #if vname == TFLAG_NAME:
    if vname in (TFLAG_NAMES):
        return False
    var = nc.variables[vname]
    dims = getattr(var, "dimensions", ())
    if len(dims) >= 4 and dims[:4] == REQ_DIMS_4D:
        return True
    if len(dims) >= 3 and dims[:3] == REQ_DIMS_3D:
        return True
    return False

def _detect_data_vars(nc, include_statics: bool = False) -> List[str]:
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

# ----------------- validation (public) -----------------

def validate(nc, include_statics: bool = False) -> Dict[str, Any]:
    """
    Return a report describing IOAPI invariants and problems, without modifying nc.
    """
    issues: List[str] = []
    report: Dict[str, Any] = {"issues": issues}

    report["shape"] = {k: _get_len(nc, k) for k in ("TSTEP", "LAY", "ROW", "COL")}

    data_vars = _detect_data_vars(nc, include_statics=include_statics)
    report["data_vars"] = data_vars
    nvars_should = len(data_vars)

    nva = getattr(nc, "NVARS", None)
    if nva is None:
        issues.append("missing-attr:NVARS")
    elif nva != nvars_should:
        issues.append(f"bad-attr:NVARS={nva} != {nvars_should}")

    var_len = _get_len(nc, "VAR")
    if var_len is None:
        issues.append("missing-dimension:VAR")
    elif var_len != nvars_should:
        issues.append(f"bad-dim:VAR={var_len} != {nvars_should}")

    expected = _build_varlist(data_vars)
    vl = getattr(nc, "VAR-LIST", None)
    if vl is None:
        issues.append("missing-attr:VAR-LIST")
    else:
        if len(vl) != len(expected):
            issues.append(f"bad-attr:VAR-LIST-length={len(vl)} expected {len(expected)}")
        elif _chunk16(vl) != data_vars:
            issues.append("bad-attr:VAR-LIST-names-mismatch")

    dt = _get_len(nc, "DATE-TIME")
    if dt is None:
        issues.append("missing-dimension:DATE-TIME")
    elif dt != 2:
        issues.append(f"bad-dim:DATE-TIME={dt} != 2")

    if TFLAG_NAME not in nc.variables:
        issues.append("missing-var:TFLAG")
    else:
        tf = nc.variables[TFLAG_NAME]
        if getattr(tf, "dimensions", ())[:3] != ("TSTEP", "VAR", "DATE-TIME"):
            issues.append(f"bad-var:TFLAG-dims={getattr(tf, 'dimensions', ())}")

    return report

# ----------------- shared finisher: NVARS/VAR-LIST + TFLAG -----------------

def _finalize_metadata(dst, data_vars: List[str], src_for_tflag=None):
    """
    Set NVARS/VAR-LIST on 'dst' and ensure TFLAG is created/preserved.
    If 'src_for_tflag' provided, attempt to preserve/clone from it; else use 'dst'.
    """
    nvars_new = len(data_vars)

    # Globals
    try: setattr(dst, "NVARS", nvars_new)
    except Exception: pass
    try: setattr(dst, "VAR-LIST", _build_varlist(data_vars))
    except Exception: pass

    # Ensure dims exist on dst
    if "VAR" not in dst.dimensions:
        dst.createDimension("VAR", nvars_new)
    if "DATE-TIME" not in dst.dimensions:
        dst.createDimension("DATE-TIME", 2)

    # Ensure/create TFLAG on dst
    tfs_dst = {}
    for tn in TFLAG_NAMES:
        if tn not in dst.variables:
            tf_dst = dst.createVariable(tn, "i4", ("TSTEP", "VAR", "DATE-TIME"))
            tf_dst.long_name = "Timestep start date and time"
            tf_dst.units     = "YYYYDDD,HHMMSS"
            tf_dst.var_desc  = "Timestep start date and time"
        else:
            tf_dst = dst.variables[tn]
        tfs_dst[tn] = tf_dst

    print(tfs_dst)
    # Try to preserve from source
    src = src_for_tflag or dst
    tf_ok = False
    for tn in TFLAG_NAMES[::-1]:
        print(tn)
        tf_dst = tfs_dst.get(tn, None)
        if not tf_dst: continue
        if tn in getattr(src, "variables", {}):
            tf_src = src.variables[tn]
            print(tf_src)
            dims_src = getattr(tf_src, "dimensions", ())
            arr = None
            if len(dims_src) >= 3 and dims_src[:3] == ("TSTEP", "VAR", "DATE-TIME"):
                arr = np.array(tf_src[:], copy=False)
            if arr is not None and arr.ndim >= 3 and arr.shape[2] == 2:
                # Must match destination TSTEP length to preserve
                tstep_dst = _get_len(dst, "TSTEP")
                if tstep_dst is None or arr.shape[0] != tstep_dst:
                    # different TSTEP length -> cannot safely preserve
                    tf_dst[:] = 0
                    return
                first = arr[:, 0, :]
                equal_all = np.all(arr == first[:, None, :])
                if not equal_all:
                    raise ValueError("Inconsistent TFLAG across VAR; refusing to guess.")
                tiled = np.tile(first[:, None, :], (1, nvars_new, 1))
                tf_dst[:] = tiled
                if tn == TFLAG_NAME:
                    tf_ok = True

    # fallback: zeros
    if not tf_ok:
        tf_dst[:] = 0
    return

# ----------------- in-place sanitize (uses finisher) -----------------

def sanitize_inplace(nc, include_statics: bool = False) -> Dict[str, Any]:
    """
    Modify 'nc' in-place to satisfy IOAPI invariants.
    Refuses to resize VAR or DATE-TIME; use the copy-clean path for those.
    """
    data_vars = _detect_data_vars(nc, include_statics=include_statics)
    nvars_new = len(data_vars)

    var_len = _get_len(nc, "VAR")
    if var_len is not None and var_len != nvars_new:
        raise ValueError(
            f"in-place sanitize would require resizing VAR from {var_len} to {nvars_new}; "
            f"use the copy-clean path instead."
        )

    dt_len = _get_len(nc, "DATE-TIME")
    if dt_len is None:
        nc.createDimension("DATE-TIME", 2)
    elif dt_len != 2:
        raise ValueError("DATE-TIME has wrong length; copy-clean required.")

    if "VAR" not in nc.dimensions:
        nc.createDimension("VAR", nvars_new)

    _finalize_metadata(nc, data_vars, src_for_tflag=nc)
    return validate(nc, include_statics=include_statics)

# ----------------- copy-clean (uses finisher) -----------------

def _copy_clean(infile: str, outfile: str, include_statics: bool = False, strict: bool = False) -> Dict[str, Any]:
    if Dataset is None:
        raise RuntimeError("netCDF4 is required for CLI copy-cleaning; pip install netCDF4")

    with Dataset(infile, "r") as src, Dataset(outfile, "w") as dst:
        data_vars = _detect_data_vars(src, include_statics=include_statics)
        nvars_new = len(data_vars)

        src_var_len = _get_len(src, "VAR")
        if src_var_len is not None:
            if strict and src_var_len != nvars_new:
                raise ValueError(
                    f"Source VAR={src_var_len} but detected #data-vars={nvars_new}. "
                    f"This file mixes VAR usage; refusing to guess (strict)."
                )
            # refuse to EXPAND VAR (source smaller than detected)
            if src_var_len < nvars_new:
                raise ValueError(
                    f"Source VAR={src_var_len} but detected #data-vars={nvars_new} (need to expand). "
                    f"Refusing to fabricate time flags; provide a source with complete TFLAG/VAR."
                )
            # if src_var_len >= nvars_new -> we will shrink in destination (allowed)

        # Copy all dims except VAR/DATE-TIME; enforce those after
        for dname, dim in src.dimensions.items():
            if dname in ("VAR", "DATE-TIME"):
                continue
            dst.createDimension(dname, None if dim.isunlimited() else len(dim))
        dst.createDimension("VAR", nvars_new)
        dst.createDimension("DATE-TIME", 2)

        # Copy globals except NVARS/VAR-LIST (rebuilt)
        for a in src.ncattrs():
            if a in ("NVARS", "VAR-LIST"):
                continue
            setattr(dst, a, getattr(src, a))

        # Copy variables except TFLAG or like (we handle TFLAG via finisher to preserve/validate)
        for vname, v in src.variables.items():
            #if vname == TFLAG_NAME:
            if vname in TFLAG_NAMES:
                continue
            fill = getattr(v, "_FillValue", None)
            kwargs = {}
            if fill is not None:
                kwargs["fill_value"] = fill
            outv = dst.createVariable(vname, v.dtype, v.dimensions, **kwargs)
            for attr in v.ncattrs():
                if attr == "_FillValue":
                    continue
                setattr(outv, attr, getattr(v, attr))
            outv[:] = v[:]

        # Shared finisher (build NVARS/VAR-LIST; preserve/clone TFLAG)
        _finalize_metadata(dst, data_vars, src_for_tflag=src)

        return validate(dst, include_statics=include_statics)

# ----------------- CLI -----------------

def _main():
    ap = argparse.ArgumentParser(description="Validate/fix IOAPI-like NetCDF (VAR-LIST/NVARS/TFLAG).")
    ap.add_argument("infile")
    ap.add_argument("outfile", nargs="?", help="Output path (omit with --dry-run)")
    ap.add_argument("--include-statics", action="store_true",
                    help="Include 2-D (ROW,COL) statics in VAR-LIST")
    ap.add_argument("--dry-run", action="store_true",
                    help="Only print report; do not write output")
    ap.add_argument("--strict", action="store_true",
                    help="Refuse if source VAR != #detected data-vars (default: allow shrinking)")
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

    rep = _copy_clean(args.infile, args.outfile,
                      include_statics=args.include_statics,
                      strict=args.strict)
    print("[DONE] Wrote:", args.outfile)
    print("[REPORT] data_vars:", rep["data_vars"])
    print("[REPORT] issues:", "OK" if not rep["issues"] else rep["issues"])

if __name__ == "__main__" and "pytest" not in sys.modules:
    _main()

