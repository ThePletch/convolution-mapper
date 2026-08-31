//! Contract-1 input types and C1A extraction config.

use serde::{Deserialize, Serialize};

use crate::error::{ErrorModule, PsfFieldError};
use crate::types::common::{
    check_finite, check_id, check_non_negative, check_positive, check_schema_version,
    check_stamp_size, check_vec2, Flag, SCHEMA_VERSION,
};

/// C1.6 pixel_mask bitfield stored as `u8`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PixelMaskBits;

impl PixelMaskBits {
    pub const INVALID: u8 = 1 << 0;
    pub const SATURATED: u8 = 1 << 1;
    pub const COSMIC: u8 = 1 << 2;
    pub const NEIGHBOR: u8 = 1 << 3;
}

/// C1.1 `StarRecord`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StarRecord {
    pub schema_version: String,
    pub star_id: String,
    pub exposure_id: String,
    pub session_id: String,
    pub field_xy_mm: [f64; 2],
    pub source_xy_px: [f64; 2],
    pub stamp: Vec<Vec<f64>>,
    pub variance: Vec<Vec<f64>>,
    pub centroid_xy_px: [f64; 2],
    pub pixel_mask: Vec<Vec<u8>>,
    pub flags: Vec<Flag>,
    pub flux_sum_adu: f64,
}

impl StarRecord {
    pub fn ingest(self) -> Result<Self, PsfFieldError> {
        check_schema_version(&self.schema_version)?;
        check_id("star_id", &self.star_id)?;
        check_id("exposure_id", &self.exposure_id)?;
        check_id("session_id", &self.session_id)?;
        check_vec2("field_xy_mm", self.field_xy_mm)?;
        check_vec2("source_xy_px", self.source_xy_px)?;
        check_vec2("centroid_xy_px", self.centroid_xy_px)?;
        check_finite("flux_sum_adu", self.flux_sum_adu)?;

        let s = square_len("stamp", &self.stamp)?;
        check_stamp_size(s)?;
        if square_len("variance", &self.variance)? != s {
            return Err(PsfFieldError::input(
                ErrorModule::Boundary,
                "variance shape must match stamp",
            )
            .with_star_id(&self.star_id));
        }
        if square_u8_len("pixel_mask", &self.pixel_mask)? != s {
            return Err(PsfFieldError::input(
                ErrorModule::Boundary,
                "pixel_mask shape must match stamp",
            )
            .with_star_id(&self.star_id));
        }

        let c_star = (s as f64 - 1.0) / 2.0;
        let dx = (self.centroid_xy_px[0] - c_star).abs();
        let dy = (self.centroid_xy_px[1] - c_star).abs();
        if dx.max(dy) > 0.6 {
            return Err(PsfFieldError::input_star(
                ErrorModule::Boundary,
                "centroid_xy_px is more than 0.6 px from the stamp center",
                &self.star_id,
            ));
        }

        for j in 0..s {
            for i in 0..s {
                let var = self.variance[j][i];
                let stamp = self.stamp[j][i];
                let mask = self.pixel_mask[j][i];
                if mask == 0 {
                    if !stamp.is_finite() {
                        return Err(PsfFieldError::input_star(
                            ErrorModule::Boundary,
                            format!("stamp[{j},{i}] is not finite on a valid pixel"),
                            &self.star_id,
                        ));
                    }
                    if !var.is_finite() || var <= 0.0 {
                        return Err(PsfFieldError::input_star(
                            ErrorModule::Boundary,
                            format!("variance[{j},{i}] must be finite and > 0 on a valid pixel"),
                            &self.star_id,
                        ));
                    }
                }
            }
        }

        let mut seen = std::collections::HashSet::new();
        for flag in &self.flags {
            if !seen.insert(*flag) {
                return Err(PsfFieldError::input_star(
                    ErrorModule::Boundary,
                    "flags must be unique",
                    &self.star_id,
                ));
            }
        }
        Ok(self)
    }

    #[must_use]
    pub fn stamp_size(&self) -> usize {
        self.stamp.len()
    }
}

fn square_len(name: &str, m: &[Vec<f64>]) -> Result<usize, PsfFieldError> {
    let n = m.len();
    if n == 0 {
        return Err(PsfFieldError::input(
            ErrorModule::Boundary,
            format!("{name} is empty"),
        ));
    }
    for row in m {
        if row.len() != n {
            return Err(PsfFieldError::input(
                ErrorModule::Boundary,
                format!("{name} must be square"),
            ));
        }
    }
    Ok(n)
}

fn square_u8_len(name: &str, m: &[Vec<u8>]) -> Result<usize, PsfFieldError> {
    let n = m.len();
    if n == 0 {
        return Err(PsfFieldError::input(
            ErrorModule::Boundary,
            format!("{name} is empty"),
        ));
    }
    for row in m {
        if row.len() != n {
            return Err(PsfFieldError::input(
                ErrorModule::Boundary,
                format!("{name} must be square"),
            ));
        }
    }
    Ok(n)
}

/// C1.4 `ImageMeta`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ImageMeta {
    pub schema_version: String,
    pub exposure_id: String,
    pub session_id: String,
    pub n_row: i64,
    pub n_col: i64,
    pub wavelength_m: f64,
    pub pupil_diameter_m: f64,
    pub focal_length_m: f64,
    pub pixel_scale_arcsec: f64,
    pub optical_axis_pixel: [f64; 2],
    pub gain_e_per_adu: f64,
    pub read_noise_e: f64,
    pub saturation_adu: f64,
    pub exptime_s: f64,
    #[serde(default)]
    pub known_defocus_waves: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pixel_size_m: Option<f64>,
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub plate_scale_warning: bool,
}

impl ImageMeta {
    pub fn ingest(mut self) -> Result<Self, PsfFieldError> {
        check_schema_version(&self.schema_version)?;
        check_id("exposure_id", &self.exposure_id)?;
        check_id("session_id", &self.session_id)?;
        if self.n_row < 1 || self.n_col < 1 {
            return Err(PsfFieldError::input(
                ErrorModule::Boundary,
                "n_row and n_col must be >= 1",
            ));
        }
        check_positive("wavelength_m", self.wavelength_m)?;
        check_positive("pupil_diameter_m", self.pupil_diameter_m)?;
        check_positive("focal_length_m", self.focal_length_m)?;
        check_positive("pixel_scale_arcsec", self.pixel_scale_arcsec)?;
        check_vec2("optical_axis_pixel", self.optical_axis_pixel)?;
        check_positive("gain_e_per_adu", self.gain_e_per_adu)?;
        check_non_negative("read_noise_e", self.read_noise_e)?;
        check_positive("saturation_adu", self.saturation_adu)?;
        check_positive("exptime_s", self.exptime_s)?;
        check_finite("known_defocus_waves", self.known_defocus_waves)?;

        let p_mm = mm_per_pixel(self.focal_length_m, self.pixel_scale_arcsec);
        let r_field = 0.5 * ((self.n_col as f64 * p_mm).hypot(self.n_row as f64 * p_mm));
        if r_field == 0.0 {
            return Err(PsfFieldError::input(
                ErrorModule::Boundary,
                "R_field is 0 (NOR.9)",
            ));
        }

        if let Some(pixel_size_m) = self.pixel_size_m {
            check_positive("pixel_size_m", pixel_size_m)?;
            let s_pred = pixel_size_m / self.focal_length_m;
            let s_hdr = pixel_scale_rad(self.pixel_scale_arcsec);
            let delta = (s_pred - s_hdr).abs() / s_hdr;
            if delta > 0.05 {
                return Err(PsfFieldError::input(
                    ErrorModule::Boundary,
                    format!("plate-scale inconsistency δ={delta} exceeds 0.05 (C1.4.1)"),
                ));
            }
            self.plate_scale_warning = delta > 0.01;
        }

        Ok(self)
    }

    #[must_use]
    pub fn mm_per_pixel(&self) -> f64 {
        mm_per_pixel(self.focal_length_m, self.pixel_scale_arcsec)
    }

    #[must_use]
    pub fn r_field_mm(&self) -> f64 {
        let p = self.mm_per_pixel();
        0.5 * ((self.n_col as f64 * p).hypot(self.n_row as f64 * p))
    }
}

#[must_use]
pub fn pixel_scale_rad(pixel_scale_arcsec: f64) -> f64 {
    pixel_scale_arcsec * std::f64::consts::PI / (180.0 * 3600.0)
}

#[must_use]
pub fn mm_per_pixel(focal_length_m: f64, pixel_scale_arcsec: f64) -> f64 {
    focal_length_m * pixel_scale_rad(pixel_scale_arcsec) * 1_000.0
}

/// C1.5 `PupilSpec`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PupilSpec {
    pub schema_version: String,
    pub mask: Vec<Vec<f64>>,
    pub n_pupil: i64,
    pub n_fft: i64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub amplitude: Option<Vec<Vec<f64>>>,
}

impl PupilSpec {
    pub fn ingest(self) -> Result<Self, PsfFieldError> {
        check_schema_version(&self.schema_version)?;
        if !matches!(self.n_pupil, 128 | 256 | 512) {
            return Err(PsfFieldError::input(
                ErrorModule::Boundary,
                format!("n_pupil {} is not in {{128, 256, 512}}", self.n_pupil),
            ));
        }
        if self.n_fft % self.n_pupil != 0 {
            return Err(PsfFieldError::input(
                ErrorModule::Boundary,
                "n_fft must be an integer multiple of n_pupil",
            ));
        }
        let ratio = self.n_fft / self.n_pupil;
        if !matches!(ratio, 2 | 4 | 8) {
            return Err(PsfFieldError::input(
                ErrorModule::Boundary,
                format!("n_fft / n_pupil = {ratio} is not in {{2, 4, 8}}"),
            ));
        }
        let n = square_len("mask", &self.mask)?;
        if n as i64 != self.n_pupil {
            return Err(PsfFieldError::input(
                ErrorModule::Boundary,
                "mask shape must be (n_pupil, n_pupil)",
            ));
        }
        for row in &self.mask {
            for &v in row {
                if v != 0.0 && v != 1.0 {
                    return Err(PsfFieldError::input(
                        ErrorModule::Boundary,
                        "v1 pupil mask values must be in {0, 1}",
                    ));
                }
            }
        }
        if let Some(amp) = &self.amplitude {
            if square_len("amplitude", amp)? != n {
                return Err(PsfFieldError::input(
                    ErrorModule::Boundary,
                    "amplitude shape must match mask",
                ));
            }
        }
        Ok(self)
    }
}

/// C1.8 config stored beside the stars.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Stage1InputConfig {
    pub schema_version: String,
    #[serde(default = "default_stamp_size")]
    pub stamp_size: i64,
}

fn default_stamp_size() -> i64 {
    31
}

impl Default for Stage1InputConfig {
    fn default() -> Self {
        Self {
            schema_version: SCHEMA_VERSION.to_string(),
            stamp_size: 31,
        }
    }
}

impl Stage1InputConfig {
    pub fn ingest(self) -> Result<Self, PsfFieldError> {
        check_schema_version(&self.schema_version)?;
        check_stamp_size(self.stamp_size as usize)?;
        Ok(self)
    }
}

/// C1A extraction config (schema `extraction_config.schema.json`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ExtractionConfig {
    pub schema_version: String,
    #[serde(default = "default_box_size")]
    pub box_size: i64,
    #[serde(default = "default_filter_size")]
    pub filter_size: i64,
    #[serde(default = "default_sigma_clip_sigma")]
    pub sigma_clip_sigma: f64,
    #[serde(default = "default_sigma_clip_maxiters")]
    pub sigma_clip_maxiters: i64,
    pub fwhm: f64,
    #[serde(default = "default_n_sigma")]
    pub n_sigma: f64,
    #[serde(default = "default_sharpness")]
    pub sharpness_range: [f64; 2],
    #[serde(default = "default_roundness")]
    pub roundness_range: [f64; 2],
    #[serde(default = "default_min_sep")]
    pub min_separation_fwhm: f64,
    #[serde(default = "default_snr_min")]
    pub snr_min: f64,
    #[serde(default = "default_max_selected")]
    pub max_selected: i64,
    #[serde(default = "default_selection_mode")]
    pub selection_mode: String,
    #[serde(default)]
    pub holdout_fraction: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub holdout_seed: Option<i64>,
}

fn default_box_size() -> i64 {
    64
}
fn default_filter_size() -> i64 {
    3
}
fn default_sigma_clip_sigma() -> f64 {
    3.0
}
fn default_sigma_clip_maxiters() -> i64 {
    5
}
fn default_n_sigma() -> f64 {
    5.0
}
fn default_sharpness() -> [f64; 2] {
    [0.2, 1.0]
}
fn default_roundness() -> [f64; 2] {
    [-1.0, 1.0]
}
fn default_min_sep() -> f64 {
    1.0
}
fn default_snr_min() -> f64 {
    20.0
}
fn default_max_selected() -> i64 {
    400
}
fn default_selection_mode() -> String {
    "highest_snr".to_string()
}

impl ExtractionConfig {
    pub fn ingest(self) -> Result<Self, PsfFieldError> {
        check_schema_version(&self.schema_version)?;
        if self.fwhm <= 1.0 || !self.fwhm.is_finite() {
            return Err(PsfFieldError::input(
                ErrorModule::Boundary,
                "fwhm must be finite and > 1.0",
            ));
        }
        Ok(self)
    }
}
