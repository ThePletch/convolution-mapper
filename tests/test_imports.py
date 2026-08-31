"""Smoke tests for the Python package layout."""

from psf_field import SCHEMA_VERSION
import psf_field_core


def test_schema_version() -> None:
    assert SCHEMA_VERSION == "1.0.0"


def test_core_extension_version() -> None:
    assert psf_field_core.__version__ == "0.1.0"
