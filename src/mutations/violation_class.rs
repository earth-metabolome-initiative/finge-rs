//! Eight ECFP-detectable violation classes, plus a multi-label bitset used as
//! the auxiliary classification target.

/// One ECFP-detectable chemistry violation class.
///
/// Each variant is aligned to a single ECFP invariant channel: any mutation
/// belonging to the variant is, by construction, visible in the ECFP output.
/// The numeric discriminant is the bit index used in [`ViolationLabel`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum ViolationClass {
    /// Atomic number outside the periodic table (`Z = 0` or `Z > 118`).
    ImpossibleAtomicNumber = 0,
    /// Atom degree exceeds the element's known maximum valence.
    Hypervalent = 1,
    /// Hydrogen count incompatible with the element's valence.
    ImpossibleHCount = 2,
    /// Formal charge outside the per-element plausibility window.
    ImpossibleCharge = 3,
    /// Isotope (`Δmass`) not present in the IUPAC table for that element.
    ImpossibleIsotope = 4,
    /// Ring-membership flag inconsistent with the atom's actual topology.
    ImpossibleRingFlag = 5,
    /// Bond invariant outside the RDKit-recognised set `{1, 2, 3, 12}`.
    ImpossibleBondType = 6,
    /// Cross-atom topological inconsistency (e.g. broken aromaticity across
    /// atoms of the same SSSR ring).
    TopologicalPathology = 7,
}

impl ViolationClass {
    /// Total number of violation classes.
    pub const COUNT: usize = 8;

    /// Every variant in declaration order.
    pub const ALL: [Self; Self::COUNT] = [
        Self::ImpossibleAtomicNumber,
        Self::Hypervalent,
        Self::ImpossibleHCount,
        Self::ImpossibleCharge,
        Self::ImpossibleIsotope,
        Self::ImpossibleRingFlag,
        Self::ImpossibleBondType,
        Self::TopologicalPathology,
    ];

    /// Returns the bit index used by [`ViolationLabel`] for this class.
    #[inline]
    #[must_use]
    pub const fn bit_index(self) -> u8 {
        self as u8
    }
}

/// Multi-label classification target — one bit per [`ViolationClass`].
///
/// A graph hits multiple bits when several predicates fire on the same input
/// (e.g. an atom with `q = +7` is both [`ViolationClass::ImpossibleCharge`]
/// and, after the resulting valence mismatch, [`ViolationClass::Hypervalent`]).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct ViolationLabel(u8);

impl ViolationLabel {
    /// Returns the all-zero label (no violations).
    #[inline]
    #[must_use]
    pub const fn empty() -> Self {
        Self(0)
    }

    /// Returns a label with every bit set.
    #[inline]
    #[must_use]
    pub const fn all() -> Self {
        // `ViolationClass::COUNT == 8` ⇒ every bit of the backing `u8` is set.
        Self(u8::MAX)
    }

    /// Sets the bit corresponding to `class`.
    #[inline]
    pub fn set(&mut self, class: ViolationClass) {
        self.0 |= 1u8 << class.bit_index();
    }

    /// Clears the bit corresponding to `class`.
    #[inline]
    pub fn clear(&mut self, class: ViolationClass) {
        self.0 &= !(1u8 << class.bit_index());
    }

    /// Returns whether the bit for `class` is set.
    #[inline]
    #[must_use]
    pub const fn is_set(self, class: ViolationClass) -> bool {
        (self.0 >> class.bit_index()) & 1 == 1
    }

    /// Returns whether no bits are set.
    #[inline]
    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.0 == 0
    }

    /// Returns the number of set bits.
    #[inline]
    #[must_use]
    pub const fn count(self) -> u32 {
        self.0.count_ones()
    }

    /// Returns the raw byte backing the label (bit `i` represents class `i`).
    #[inline]
    #[must_use]
    pub const fn as_u8(self) -> u8 {
        self.0
    }

    /// Returns an iterator over every set class, in declaration order.
    #[inline]
    pub fn iter(self) -> impl Iterator<Item = ViolationClass> {
        ViolationClass::ALL
            .into_iter()
            .filter(move |class| self.is_set(*class))
    }
}

impl core::ops::BitOr for ViolationLabel {
    type Output = Self;

    #[inline]
    fn bitor(self, rhs: Self) -> Self {
        Self(self.0 | rhs.0)
    }
}

impl core::ops::BitOrAssign for ViolationLabel {
    #[inline]
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}

impl From<ViolationClass> for ViolationLabel {
    #[inline]
    fn from(class: ViolationClass) -> Self {
        let mut label = Self::empty();
        label.set(class);
        label
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec::Vec;

    use super::{ViolationClass, ViolationLabel};

    #[test]
    fn bit_index_matches_discriminant() {
        for (i, class) in ViolationClass::ALL.iter().enumerate() {
            assert_eq!(usize::from(class.bit_index()), i);
        }
    }

    #[test]
    fn empty_label_has_no_bits() {
        let label = ViolationLabel::empty();
        assert!(label.is_empty());
        assert_eq!(label.count(), 0);
        for class in ViolationClass::ALL {
            assert!(!label.is_set(class));
        }
    }

    #[test]
    fn default_label_is_empty() {
        assert_eq!(ViolationLabel::default(), ViolationLabel::empty());
    }

    #[test]
    fn set_and_clear_round_trip() {
        let mut label = ViolationLabel::empty();
        label.set(ViolationClass::Hypervalent);
        label.set(ViolationClass::ImpossibleIsotope);
        assert!(label.is_set(ViolationClass::Hypervalent));
        assert!(label.is_set(ViolationClass::ImpossibleIsotope));
        assert!(!label.is_set(ViolationClass::ImpossibleCharge));
        assert_eq!(label.count(), 2);

        label.clear(ViolationClass::Hypervalent);
        assert!(!label.is_set(ViolationClass::Hypervalent));
        assert!(label.is_set(ViolationClass::ImpossibleIsotope));
        assert_eq!(label.count(), 1);
    }

    #[test]
    fn all_label_sets_every_class() {
        let label = ViolationLabel::all();
        assert_eq!(label.count() as usize, ViolationClass::COUNT);
        for class in ViolationClass::ALL {
            assert!(label.is_set(class));
        }
    }

    #[test]
    fn iter_yields_set_classes_in_declaration_order() {
        let mut label = ViolationLabel::empty();
        label.set(ViolationClass::TopologicalPathology);
        label.set(ViolationClass::ImpossibleAtomicNumber);
        label.set(ViolationClass::ImpossibleCharge);

        let collected: Vec<ViolationClass> = label.iter().collect();
        assert_eq!(
            collected,
            alloc::vec![
                ViolationClass::ImpossibleAtomicNumber,
                ViolationClass::ImpossibleCharge,
                ViolationClass::TopologicalPathology,
            ]
        );
    }

    #[test]
    fn bitor_unions_labels() {
        let a = ViolationLabel::from(ViolationClass::Hypervalent);
        let b = ViolationLabel::from(ViolationClass::ImpossibleBondType);
        let union = a | b;
        assert!(union.is_set(ViolationClass::Hypervalent));
        assert!(union.is_set(ViolationClass::ImpossibleBondType));
        assert_eq!(union.count(), 2);

        let mut acc = ViolationLabel::empty();
        acc |= a;
        acc |= b;
        assert_eq!(acc, union);
    }

    #[test]
    fn from_class_sets_only_that_bit() {
        let label = ViolationLabel::from(ViolationClass::ImpossibleRingFlag);
        assert_eq!(label.count(), 1);
        assert!(label.is_set(ViolationClass::ImpossibleRingFlag));
    }

    #[test]
    fn as_u8_round_trips_bit_pattern() {
        let mut label = ViolationLabel::empty();
        label.set(ViolationClass::ImpossibleAtomicNumber);
        label.set(ViolationClass::ImpossibleBondType);
        assert_eq!(label.as_u8(), 0b0100_0001);
    }
}
