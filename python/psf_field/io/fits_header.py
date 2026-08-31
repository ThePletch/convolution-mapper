"""Minimal FITS primary-header codec for ImageMeta ingest (C1.4.2 / C1.9.1)."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

FitsValue = str | int | float | bool

# C1.9.1 header names, plus 8-char IDs not named in the table.
FITS_TO_FIELD: dict[str, str] = {
    "SCHEMAV": "schema_version",
    "EXPID": "exposure_id",
    "SESSID": "session_id",
    "NAXIS2": "n_row",
    "NAXIS1": "n_col",
    "LAMBDA": "wavelength_m",
    "PUPILD": "pupil_diameter_m",
    "FOCALLEN": "focal_length_m",
    "PIXSCAL": "pixel_scale_arcsec",
    "OAX": "optical_axis_pixel_x",
    "OAY": "optical_axis_pixel_y",
    "GAIN": "gain_e_per_adu",
    "RDNOISE": "read_noise_e",
    "SATURATE": "saturation_adu",
    "EXPTIME": "exptime_s",
    "KDEFOCUS": "known_defocus_waves",
}


def _card(keyword: str, value: FitsValue) -> bytes:
    key = keyword[:8].ljust(8)
    if isinstance(value, bool):
        rendered = "T" if value else "F"
        body = f"{key}= {rendered:>20}"
    elif isinstance(value, int) and not isinstance(value, bool):
        body = f"{key}= {value:20d}"
    elif isinstance(value, float):
        body = f"{key}= {value:20.10G}"
    else:
        body = f"{key}= '{value}'"
    return body.ljust(80).encode("ascii")


def write_primary_header(path: Path, cards: Mapping[str, FitsValue]) -> None:
    """Write a header-only FITS file (no image data) for ingest tests."""
    records: list[bytes] = [
        _card("SIMPLE", True),
        _card("BITPIX", 8),
        _card("NAXIS", 2),
    ]
    naxis1 = cards.get("NAXIS1", 1)
    naxis2 = cards.get("NAXIS2", 1)
    records.append(_card("NAXIS1", naxis1))
    records.append(_card("NAXIS2", naxis2))
    for key, value in cards.items():
        if key in {"SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2"}:
            continue
        records.append(_card(key, value))
    records.append(b"END".ljust(80))
    block = b"".join(records)
    pad = (2880 - (len(block) % 2880)) % 2880
    path.write_bytes(block + b"\x00" * pad)


def read_primary_header(path: Path) -> dict[str, FitsValue]:
    raw = path.read_bytes()
    if len(raw) < 2880:
        msg = "FITS header is shorter than one 2880-byte block"
        raise ValueError(msg)
    out: dict[str, FitsValue] = {}
    for offset in range(0, len(raw), 80):
        card = raw[offset : offset + 80].decode("ascii", errors="replace")
        if card.startswith("END"):
            break
        if len(card) < 10 or card[8] != "=":
            continue
        key = card[:8].strip()
        payload = card[10:].strip()
        out[key] = _parse_value(payload)
    return out


def _parse_value(payload: str) -> FitsValue:
    comment = payload.split("/", 1)[0].strip()
    if comment in {"T", "F"}:
        return comment == "T"
    if comment.startswith("'"):
        end = comment.find("'", 1)
        return comment[1:end] if end > 0 else comment.strip("'")
    if "." in comment or "E" in comment.upper():
        return float(comment)
    try:
        return int(comment)
    except ValueError:
        return comment


def header_to_image_meta_fields(header: Mapping[str, FitsValue]) -> dict[str, object]:
    """Map FITS keywords onto ImageMeta JSON field names. Does not invent defaults."""
    fields: dict[str, object] = {}
    extras: dict[str, object] = {}
    for key, value in header.items():
        mapped = FITS_TO_FIELD.get(key)
        if mapped is None:
            continue
        extras[mapped] = value
    axis_x = extras.pop("optical_axis_pixel_x", None)
    axis_y = extras.pop("optical_axis_pixel_y", None)
    if axis_x is not None and axis_y is not None:
        fields["optical_axis_pixel"] = [axis_x, axis_y]
    fields.update(extras)
    return fields
