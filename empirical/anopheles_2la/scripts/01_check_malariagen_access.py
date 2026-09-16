#!/usr/bin/env python3
"""Check metadata-only MalariaGEN Ag3 access for Anopheles 2La Stage 1A.

This script deliberately stops at release/sample-set discovery. It must not
touch haplotypes, SNP calls, trees, or other genomic arrays.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import re
import sys
from pathlib import Path
from typing import Any


TARGET_RELEASE = "3.10"
TARGET_SAMPLE_SET = "fontaine-2015-rebuild"

STATUS_AUTH = "AUTHENTICATION_REQUIRED"
STATUS_API = "API_UNAVAILABLE"
STATUS_MISSING = "SAMPLE_SET_NOT_FOUND"
STATUS_SUCCESS = "SUCCESS"


def package_version() -> str:
    try:
        return importlib.metadata.version("malariagen-data")
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def sanitize_error(exc: BaseException) -> str:
    text = str(exc)
    text = re.sub(r"Bearer\s+[A-Za-z0-9._~+/=-]+", "Bearer <redacted>", text, flags=re.I)
    text = re.sub(r"access[_ -]?token[=:]\s*[^,\s]+", "access_token=<redacted>", text, flags=re.I)
    text = re.sub(r"refresh[_ -]?token[=:]\s*[^,\s]+", "refresh_token=<redacted>", text, flags=re.I)
    text = re.sub(r"token=([^&\s]+)", "token=<redacted>", text, flags=re.I)
    text = re.sub(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}", "<redacted-account>", text)
    text = re.sub(r"(/[^\s,;:]+)+/(?:[^/\s,;:]+\.json)", "<redacted-credential-path>", text)
    return text


def is_auth_error(exc: BaseException) -> bool:
    text = f"{type(exc).__name__}: {exc}".lower()
    markers = (
        "401",
        "403",
        "unauthorized",
        "forbidden",
        "permission",
        "anonymous caller",
        "credentials",
        "authentication",
        "storage.objects.get access",
    )
    return any(marker in text for marker in markers)


def find_column(columns: list[str], candidates: list[str]) -> str:
    by_lower = {c.lower(): c for c in columns}
    for candidate in candidates:
        if candidate.lower() in by_lower:
            return by_lower[candidate.lower()]
    return ""


def available_releases(ag3: Any) -> list[str]:
    value = getattr(ag3, "releases", None)
    if value is None:
        return []
    if callable(value):
        value = value()
    try:
        return [str(item) for item in value]
    except TypeError:
        return []


def print_auth_message() -> None:
    print(
        "MalariaGEN GCS authentication is required.\n"
        "Outside Google Colab, configure Google Application Default Credentials.\n"
        "Recommended command:\n"
        "    gcloud auth application-default login",
        file=sys.stderr,
    )


def check(args: argparse.Namespace) -> int:
    print(f"Python: {sys.version.split()[0]}")
    print(f"malariagen-data: {package_version()}")
    try:
        import malariagen_data

        ag3_kwargs = {"check_location": False, "show_progress": False}
        if args.url:
            ag3_kwargs["url"] = args.url
        print("Ag3 initialization: starting")
        ag3 = malariagen_data.Ag3(**ag3_kwargs)
        print("Ag3 initialization: ok")

        releases = available_releases(ag3)
        if releases:
            print(f"Available releases: {', '.join(releases)}")
            print(f"Release {TARGET_RELEASE} visible: {'yes' if TARGET_RELEASE in releases else 'no'}")
        else:
            print("Available releases: not exposed by this package version; checking sample_sets directly")

        sample_sets = ag3.sample_sets(release=TARGET_RELEASE)
        id_col = find_column(list(sample_sets.columns), ["sample_set", "id"])
        if not id_col:
            print(f"Stage 1A access status: {STATUS_API}", file=sys.stderr)
            print(f"Could not identify sample-set column in {list(sample_sets.columns)!r}", file=sys.stderr)
            return 2
        sample_set_ids = [str(value) for value in sample_sets[id_col].tolist()]
        release_visible = bool(sample_set_ids)
        found = TARGET_SAMPLE_SET in sample_set_ids
        if not found:
            found = any(TARGET_SAMPLE_SET.lower() in value.lower() for value in sample_set_ids)
        print(f"Release {TARGET_RELEASE} sample-set metadata visible: {'yes' if release_visible else 'no'}")
        print(f"{TARGET_SAMPLE_SET} located: {'yes' if found else 'no'}")
        if not found:
            print(f"Stage 1A access status: {STATUS_MISSING}", file=sys.stderr)
            return 3
        print(f"Stage 1A access status: {STATUS_SUCCESS}")
        return 0
    except Exception as exc:
        if is_auth_error(exc):
            print(f"Stage 1A access status: {STATUS_AUTH}", file=sys.stderr)
            print_auth_message()
            return 2
        print(f"Stage 1A access status: {STATUS_API}", file=sys.stderr)
        print(f"{type(exc).__name__}: {sanitize_error(exc)}", file=sys.stderr)
        return 2


def main() -> int:
    parser = argparse.ArgumentParser(description="Check authenticated metadata-only MalariaGEN Ag3 access for Stage 1A.")
    parser.add_argument("--url", default="", help="Optional MalariaGEN API storage URL. Leave empty to use the package default.")
    args = parser.parse_args()
    return check(args)


if __name__ == "__main__":
    raise SystemExit(main())
