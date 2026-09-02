//! Ingest validators for Contract-1 records, ImageMeta sidecars, and pupil specs.

use serde_json::{Map, Value};

use crate::error::{ErrorModule, PsfFieldError};
use crate::types::{Catalog, ImageMeta, PupilSpec, StarRecord};

const FIELD_TOL_MM: f64 = 1e-6;

/// Deserialize JSON and run C1 ingest checks.
pub fn ingest_star_record(value: &Value) -> Result<StarRecord, PsfFieldError> {
    let rec: StarRecord = serde_json::from_value(value.clone()).map_err(|e| {
        PsfFieldError::input(ErrorModule::Boundary, format!("StarRecord JSON: {e}"))
    })?;
    rec.ingest()
}

pub fn ingest_image_meta(value: &Value) -> Result<ImageMeta, PsfFieldError> {
    let meta: ImageMeta = serde_json::from_value(value.clone())
        .map_err(|e| PsfFieldError::input(ErrorModule::Boundary, format!("ImageMeta JSON: {e}")))?;
    meta.ingest()
}

pub fn ingest_pupil_spec(value: &Value) -> Result<PupilSpec, PsfFieldError> {
    let pupil: PupilSpec = serde_json::from_value(value.clone())
        .map_err(|e| PsfFieldError::input(ErrorModule::Boundary, format!("PupilSpec JSON: {e}")))?;
    pupil.ingest()
}

pub fn ingest_catalog(value: &Value) -> Result<Catalog, PsfFieldError> {
    let catalog: Catalog = serde_json::from_value(value.clone())
        .map_err(|e| PsfFieldError::input(ErrorModule::Boundary, format!("Catalog JSON: {e}")))?;
    catalog.ingest()
}

/// Merge a JSON sidecar over a base ImageMeta object, then ingest (C1.4.2).
///
/// Sidecar keys overlay the base. Missing NOR.13 fields after merge are rejected
/// by `ImageMeta::ingest`. Implementations SHALL NOT invent `optical_axis_pixel`.
pub fn merge_image_meta_sidecar(
    mut base: Map<String, Value>,
    sidecar: Map<String, Value>,
) -> Result<ImageMeta, PsfFieldError> {
    for (k, v) in sidecar {
        base.insert(k, v);
    }
    ingest_image_meta(&Value::Object(base))
}

/// Session-level C1.10 checks: shared S, unique IDs, field coordinates, detector size.
pub fn ingest_session(
    stars: Vec<StarRecord>,
    meta: ImageMeta,
    pupil: PupilSpec,
) -> Result<(Vec<StarRecord>, ImageMeta, PupilSpec), PsfFieldError> {
    let meta = meta.ingest()?;
    let pupil = pupil.ingest()?;
    crate::resample::oversampling_factor(&meta, &pupil)?;
    if stars.is_empty() {
        return Err(PsfFieldError::input(
            ErrorModule::Boundary,
            "session has no stars",
        ));
    }

    let mut star_ids = std::collections::HashSet::new();
    let mut stamp_size: Option<usize> = None;
    let p_mm = meta.mm_per_pixel();
    let mut ingested = Vec::with_capacity(stars.len());

    for star in stars {
        let star = star.ingest()?;
        if star.exposure_id != meta.exposure_id {
            return Err(PsfFieldError::input_star(
                ErrorModule::Boundary,
                "star exposure_id does not match ImageMeta",
                &star.star_id,
            ));
        }
        if star.session_id != meta.session_id {
            return Err(PsfFieldError::input_star(
                ErrorModule::Boundary,
                "star session_id does not match ImageMeta",
                &star.star_id,
            ));
        }
        if !star_ids.insert(star.star_id.clone()) {
            return Err(PsfFieldError::input_star(
                ErrorModule::Boundary,
                "duplicate star_id",
                &star.star_id,
            ));
        }
        let s = star.stamp_size();
        match stamp_size {
            None => stamp_size = Some(s),
            Some(expected) if expected != s => {
                return Err(PsfFieldError::input(
                    ErrorModule::Boundary,
                    "mixed stamp sizes in one exposure_id",
                ));
            }
            Some(_) => {}
        }
        if meta.n_row < s as i64 || meta.n_col < s as i64 {
            return Err(PsfFieldError::input(
                ErrorModule::Boundary,
                "n_row and n_col must be >= S",
            ));
        }

        let recomputed = [
            (star.source_xy_px[0] - meta.optical_axis_pixel[0]) * p_mm,
            (star.source_xy_px[1] - meta.optical_axis_pixel[1]) * p_mm,
        ];
        let dx = (star.field_xy_mm[0] - recomputed[0]).abs();
        let dy = (star.field_xy_mm[1] - recomputed[1]).abs();
        if dx.max(dy) > FIELD_TOL_MM {
            return Err(PsfFieldError::input_star(
                ErrorModule::Boundary,
                format!(
                    "field_xy_mm does not match NOR.8 recomputation (Δ∞={} mm)",
                    dx.max(dy)
                ),
                &star.star_id,
            ));
        }
        ingested.push(star);
    }

    Ok((ingested, meta, pupil))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{ImageMeta, StarRecord};
    use serde_json::json;

    fn valid_meta() -> Value {
        json!({
            "schema_version": "1.0.0",
            "exposure_id": "exp1",
            "session_id": "sess1",
            "n_row": 1024,
            "n_col": 1024,
            "wavelength_m": 5.5e-7,
            "pupil_diameter_m": 0.2,
            "focal_length_m": 1.6,
            "pixel_scale_arcsec": 0.4,
            "optical_axis_pixel": [511.5, 511.5],
            "gain_e_per_adu": 1.5,
            "read_noise_e": 5.0,
            "saturation_adu": 60000.0,
            "exptime_s": 30.0
        })
    }

    fn ones(s: usize, v: f64) -> Vec<Vec<f64>> {
        vec![vec![v; s]; s]
    }

    fn zeros_u8(s: usize) -> Vec<Vec<u8>> {
        vec![vec![0_u8; s]; s]
    }

    fn valid_star(s: usize) -> Value {
        let c = (s as f64 - 1.0) / 2.0;
        json!({
            "schema_version": "1.0.0",
            "star_id": "s1",
            "exposure_id": "exp1",
            "session_id": "sess1",
            "field_xy_mm": [0.0, 0.0],
            "source_xy_px": [511.5, 511.5],
            "stamp": ones(s, 1.0),
            "variance": ones(s, 1.0),
            "centroid_xy_px": [c, c],
            "pixel_mask": zeros_u8(s),
            "flags": ["SELECTED"],
            "flux_sum_adu": (s * s) as f64
        })
    }

    #[test]
    fn rejects_unknown_major() {
        let mut v = valid_star(15);
        v["schema_version"] = json!("2.0.0");
        let err = ingest_star_record(&v).unwrap_err();
        assert!(err.message.contains("unknown schema MAJOR 2"));
    }

    #[test]
    fn rejects_even_stamp_size() {
        let v = valid_star(16);
        let err = ingest_star_record(&v).unwrap_err();
        assert!(err.message.contains("even stamp size"));
    }

    #[test]
    fn sidecar_incomplete_image_meta() {
        let mut base = valid_meta().as_object().cloned().unwrap();
        base.remove("wavelength_m");
        base.remove("optical_axis_pixel");
        let err = merge_image_meta_sidecar(base.clone(), Map::new()).unwrap_err();
        assert_eq!(err.code, crate::error::ErrorCode::Input);

        let mut sidecar = Map::new();
        sidecar.insert("wavelength_m".into(), json!(5.5e-7));
        sidecar.insert("optical_axis_pixel".into(), json!([511.5, 511.5]));
        merge_image_meta_sidecar(base, sidecar).unwrap();
    }

    #[test]
    fn default_catalog_round_trips() {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../docs/contracts/schemas/psf_field_v1_default.catalog.json");
        let text = std::fs::read_to_string(path).unwrap();
        let value: Value = serde_json::from_str(&text).unwrap();
        let catalog = ingest_catalog(&value).unwrap();
        let back = serde_json::to_value(&catalog).unwrap();
        let again = ingest_catalog(&back).unwrap();
        assert_eq!(catalog, again);
        assert_eq!(catalog.catalog_id, "psf_field_v1_default");
    }

    #[test]
    fn ingest_zeros_amplitude_outside_the_mask() {
        let mut pupil = crate::pupil::circular_pupil_spec(128, 512);
        pupil.amplitude = Some(vec![vec![1.0; 128]; 128]);
        let ingested = pupil.ingest().unwrap();
        let amplitude = ingested.amplitude.as_ref().unwrap();
        for (mask_row, amplitude_row) in ingested.mask.iter().zip(amplitude) {
            for (&mask, &value) in mask_row.iter().zip(amplitude_row) {
                if mask == 0.0 {
                    assert_eq!(value, 0.0);
                } else {
                    assert_eq!(value, 1.0);
                }
            }
        }
    }

    #[test]
    fn session_rejects_insane_sampling() {
        let mut meta = ingest_image_meta(&valid_meta()).unwrap();
        meta.pixel_scale_arcsec = 0.01;
        let pupil = crate::pupil::circular_pupil_spec(128, 512)
            .ingest()
            .unwrap();
        let star: StarRecord = serde_json::from_value(valid_star(15)).unwrap();
        let err = ingest_session(vec![star], meta, pupil).unwrap_err();
        assert_eq!(err.code, crate::error::ErrorCode::Input);
        assert!(err.message.contains("sampling insane"));
    }

    #[test]
    fn session_accepts_c10_1_sampling() {
        let mut meta = ImageMeta::c10_1_standard_camera();
        meta.exposure_id = "exp1".to_string();
        meta.session_id = "sess1".to_string();
        let meta = meta.ingest().unwrap();
        let pupil = crate::pupil::circular_pupil_spec(256, 1024)
            .ingest()
            .unwrap();
        let mut star = valid_star(31);
        star["source_xy_px"] = json!([255.5, 255.5]);
        star["centroid_xy_px"] = json!([15.0, 15.0]);
        let star: StarRecord = serde_json::from_value(star).unwrap();
        ingest_session(vec![star], meta, pupil).unwrap();
    }
}
