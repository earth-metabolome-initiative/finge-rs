//! ECFP-detectable mutation operators for auxiliary negative-sampling tasks.
//!
//! See `docs/plan/` and `docs/research/ecfp-negative-sampling-survey.md` in the
//! repository for the design rationale.

pub mod violation_class;

pub use violation_class::{ViolationClass, ViolationLabel};
