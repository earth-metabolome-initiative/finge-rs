//! Pre-hash field tuple and shared RDKit-style hash for ECFP atom invariants.
//!
//! `AtomInvariantFields` is the canonical decomposition of the RDKit Morgan
//! atom invariant — `[atomic_number, total_degree, total_hydrogens,
//! formal_charge, delta_mass, in_ring]`. The hash that turns this tuple into
//! a 32-bit ECFP invariant is the boost-style `hash_combine` chain that
//! every Morgan-compatible fingerprint implementation uses.
//!
//! Centralising both lets the mutation infrastructure perturb one chemistry
//! field (e.g. set `atomic_number = 200`) and re-hash through the same code
//! path that the standard ECFP pipeline uses on real molecules.

/// Decomposed RDKit-style ECFP atom invariant fields.
///
/// The values mirror RDKit's `getConnectivityInvariants` inputs:
/// `atomic_number` is the periodic-table Z, `total_degree` is heavy degree +
/// total hydrogens, `total_hydrogens` is explicit + implicit hydrogens,
/// `formal_charge` is the (RDKit-normalised) integer charge, `delta_mass` is
/// `isotope_mass - standard_atomic_mass` (zero for the standard isotope),
/// and `in_ring` records ring-atom membership.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct AtomInvariantFields {
    /// Atomic number (`Z`). `0` denotes "no element" / wildcard atoms.
    pub atomic_number: u32,
    /// Total degree counted as `heavy_degree + total_hydrogens`.
    pub total_degree: u32,
    /// Total hydrogen count (explicit + implicit).
    pub total_hydrogens: u32,
    /// Formal charge (signed). Stored as `i32` to accommodate mutations that
    /// push outside the natural `i8` range.
    pub formal_charge: i32,
    /// Isotope mass delta `isotope_mass - standard_atomic_mass`. `0` for the
    /// standard isotope.
    pub delta_mass: i32,
    /// Whether the atom is a member of any ring.
    pub in_ring: bool,
}

impl AtomInvariantFields {
    /// Hashes this field tuple to a 32-bit ECFP atom invariant.
    ///
    /// When `include_ring_membership` is `true` *and* `self.in_ring` is set,
    /// the ring flag (`1`) is appended to the hashed sequence; otherwise it
    /// is omitted. This matches RDKit's behaviour for the `useRingMembership`
    /// flag on the Morgan generator.
    #[inline]
    #[must_use]
    pub fn to_invariant(self, include_ring_membership: bool) -> u32 {
        let charge_bits = self.formal_charge as u32;
        let delta_mass_bits = self.delta_mass as u32;
        if include_ring_membership && self.in_ring {
            morgan_hash_u32_sequence(&[
                self.atomic_number,
                self.total_degree,
                self.total_hydrogens,
                charge_bits,
                delta_mass_bits,
                1,
            ])
        } else {
            morgan_hash_u32_sequence(&[
                self.atomic_number,
                self.total_degree,
                self.total_hydrogens,
                charge_bits,
                delta_mass_bits,
            ])
        }
    }
}

/// Boost-style `hash_combine` used by RDKit's Morgan fingerprint.
///
/// `seed ^= value.wrapping_add(0x9e37_79b9)
///   .wrapping_add(seed << 6)
///   .wrapping_add(seed >> 2)` (all arithmetic wraps modulo 2³²).
#[inline]
pub fn morgan_hash_combine(seed: &mut u32, value: u32) {
    *seed ^= value
        .wrapping_add(0x9e37_79b9)
        .wrapping_add(seed.wrapping_shl(6))
        .wrapping_add(seed.wrapping_shr(2));
}

/// Folds a sequence of `u32` values into a single 32-bit hash with
/// [`morgan_hash_combine`].
#[inline]
#[must_use]
pub fn morgan_hash_u32_sequence(values: &[u32]) -> u32 {
    let mut seed = 0_u32;
    for &value in values {
        morgan_hash_combine(&mut seed, value);
    }
    seed
}

#[cfg(test)]
mod tests {
    use super::{AtomInvariantFields, morgan_hash_combine, morgan_hash_u32_sequence};

    #[test]
    fn hash_combine_matches_boost_formula_from_zero_seed() {
        // For seed = 0, value = 1: result = 0 ^ (1 + 0x9e3779b9 + 0 + 0)
        let mut seed = 0_u32;
        morgan_hash_combine(&mut seed, 1);
        assert_eq!(seed, 0x9e37_79ba);
    }

    #[test]
    fn hash_combine_chains_match_manual_application() {
        let mut chained = 0_u32;
        morgan_hash_combine(&mut chained, 6);
        morgan_hash_combine(&mut chained, 4);
        morgan_hash_combine(&mut chained, 4);
        assert_eq!(morgan_hash_u32_sequence(&[6, 4, 4]), chained);
    }

    #[test]
    fn sequence_depends_on_every_position() {
        let baseline = morgan_hash_u32_sequence(&[6, 4, 4, 0, 0]);
        for i in 0..5 {
            let mut values: [u32; 5] = [6, 4, 4, 0, 0];
            values[i] = values[i].wrapping_add(1);
            let perturbed = morgan_hash_u32_sequence(&values);
            assert_ne!(
                baseline, perturbed,
                "perturbing position {i} changed nothing"
            );
        }
    }

    #[test]
    fn to_invariant_with_ring_differs_from_without() {
        let fields = AtomInvariantFields {
            atomic_number: 6,
            total_degree: 4,
            total_hydrogens: 4,
            formal_charge: 0,
            delta_mass: 0,
            in_ring: true,
        };
        assert_ne!(fields.to_invariant(true), fields.to_invariant(false));
    }

    #[test]
    fn to_invariant_ignores_in_ring_when_flag_off() {
        let with_ring = AtomInvariantFields {
            atomic_number: 6,
            total_degree: 4,
            total_hydrogens: 4,
            formal_charge: 0,
            delta_mass: 0,
            in_ring: true,
        };
        let without_ring = AtomInvariantFields {
            in_ring: false,
            ..with_ring
        };
        assert_eq!(
            with_ring.to_invariant(false),
            without_ring.to_invariant(false),
        );
    }

    #[test]
    fn to_invariant_changes_with_atomic_number() {
        let carbon = AtomInvariantFields {
            atomic_number: 6,
            total_degree: 4,
            total_hydrogens: 4,
            formal_charge: 0,
            delta_mass: 0,
            in_ring: false,
        };
        let alien = AtomInvariantFields {
            atomic_number: 200,
            ..carbon
        };
        assert_ne!(carbon.to_invariant(true), alien.to_invariant(true));
    }

    #[test]
    fn negative_charge_round_trips_via_u32_cast() {
        // Casting a small negative i32 to u32 preserves the bit pattern; the
        // hash must therefore differ from a same-magnitude positive charge.
        let neg = AtomInvariantFields {
            atomic_number: 7,
            total_degree: 3,
            total_hydrogens: 1,
            formal_charge: -1,
            delta_mass: 0,
            in_ring: false,
        };
        let pos = AtomInvariantFields {
            formal_charge: 1,
            ..neg
        };
        assert_ne!(neg.to_invariant(true), pos.to_invariant(true));
    }
}
