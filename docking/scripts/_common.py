"""Shared helpers for the numbered pipeline steps.

Deliberately small: path resolution, config loading, logging and a polite HTTP
getter. Anything specific to one step belongs in that step's script.
"""

from __future__ import annotations

import hashlib
import logging
import os
import sys
import time
from pathlib import Path
from typing import Any

import requests
import yaml

# Generic client name on purpose: no personal data goes to external services.
USER_AGENT = "taxane-pgp-docking/0.1 (research pipeline; python-requests)"

_DEFAULT_TIMEOUT = 60
_MAX_RETRIES = 4


def repo_root() -> Path:
    """Project root, i.e. the directory that holds config/ and scripts/."""
    env = os.environ.get("DOCKING_ROOT")
    if env:
        return Path(env).resolve()
    return Path(__file__).resolve().parent.parent


def load_yaml(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    if not path.is_absolute():
        path = repo_root() / path
    with path.open() as handle:
        return yaml.safe_load(handle)


def setup_logging(step_name: str, verbose: bool = False) -> logging.Logger:
    """Log to stderr and to logs/<step_name>.log (appended, so runs accumulate)."""
    log_dir = repo_root() / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger(step_name)
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.handlers.clear()

    fmt = logging.Formatter(
        "%(asctime)s %(levelname)-7s %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    stream = logging.StreamHandler(sys.stderr)
    stream.setFormatter(fmt)
    logger.addHandler(stream)

    file_handler = logging.FileHandler(log_dir / f"{step_name}.log")
    file_handler.setFormatter(fmt)
    logger.addHandler(file_handler)

    logger.propagate = False
    return logger


def http_get(
    url: str,
    *,
    logger: logging.Logger | None = None,
    min_interval: float = 0.25,
    accept: str | None = None,
    allow_404: bool = False,
) -> requests.Response | None:
    """GET with retries and a floor on request spacing.

    `min_interval` keeps us inside PubChem's 5 requests/second policy. Returns
    None only when `allow_404` is set and the server answered 404.
    """
    headers = {"User-Agent": USER_AGENT}
    if accept:
        headers["Accept"] = accept

    last_error: Exception | None = None
    for attempt in range(1, _MAX_RETRIES + 1):
        try:
            response = requests.get(url, headers=headers, timeout=_DEFAULT_TIMEOUT)
        except requests.RequestException as exc:  # network hiccup
            last_error = exc
            if logger:
                logger.warning("attempt %d/%d failed for %s: %s", attempt, _MAX_RETRIES, url, exc)
            time.sleep(min(2 ** attempt, 15))
            continue

        if response.status_code == 404 and allow_404:
            time.sleep(min_interval)
            return None
        if response.status_code in (429, 503):
            wait = min(2 ** attempt, 15)
            if logger:
                logger.warning("throttled (%d) on %s, waiting %ds", response.status_code, url, wait)
            time.sleep(wait)
            continue
        if response.ok:
            time.sleep(min_interval)
            return response

        last_error = RuntimeError(f"HTTP {response.status_code} for {url}")
        if logger:
            logger.warning("%s", last_error)
        time.sleep(min(2 ** attempt, 15))

    raise RuntimeError(f"giving up on {url}") from last_error


def sha256_file(path: str | Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()
