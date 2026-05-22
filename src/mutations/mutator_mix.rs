//! [`MutatorMix`] — weighted-random sampler over a collection of mutators
//! that also computes the multi-label [`ViolationLabel`] for each generated
//! negative.

use alloc::{boxed::Box, vec::Vec};

use rand_core::RngCore;

use crate::{
    mutations::{
        invalidated_graph::InvalidatedGraph,
        mutator::{Mutator, MutatorError},
        mutators::{
            HypervalentMutator, ImpossibleAtomicNumberMutator, ImpossibleBondTypeMutator,
            ImpossibleChargeMutator, ImpossibleHCountMutator, ImpossibleIsotopeMutator,
            ImpossibleRingFlagMutator, TopologicalPathologyMutator,
        },
        predicate::{
            HypervalentPredicate, ImpossibleAtomicNumberPredicate, ImpossibleBondTypePredicate,
            ImpossibleChargePredicate, ImpossibleHCountPredicate, ImpossibleIsotopePredicate,
            ImpossibleRingFlagPredicate, TopologicalPathologyPredicate, ViolationPredicate,
        },
        violation_class::ViolationLabel,
    },
    traits::{EcfpGraph, MolecularAtom},
};

/// One entry in the mix: a mutator and its sampling weight.
type WeightedMutator<'mix, G> = (Box<dyn Mutator<G> + 'mix>, f32);

/// Weighted-random sampler over a configurable set of mutators plus a fixed
/// list of predicates that label each generated negative.
///
/// Build with [`MutatorMix::new`] (empty) or
/// [`MutatorMix::with_default_mutators_and_predicates`] (the eight built-in
/// pairs, each with weight `1.0`). Append more via [`MutatorMix::add_mutator`]
/// and [`MutatorMix::add_predicate`].
///
/// Sampling:
/// 1. Pick a mutator by weight.
/// 2. Apply it to the input graph.
/// 3. Run every predicate on the resulting wrapper and OR the hits into a
///    [`ViolationLabel`].
///
/// `sample()` returns `Err` when the chosen mutator declines the input. The
/// caller is expected to retry with another input (or, if iterating a fixed
/// corpus, simply skip and move on).
///
/// The `'mix` lifetime parameter lets the mix carry mutators / predicates
/// that borrow data shorter than `'static` (e.g. when the graph type `G`
/// itself contains a borrow, as with `SmilesRdkitGraph<'_>`).
pub struct MutatorMix<'mix, G: EcfpGraph>
where
    G::NodeSymbol: MolecularAtom,
{
    mutators: Vec<WeightedMutator<'mix, G>>,
    predicates: Vec<Box<dyn ViolationPredicate<InvalidatedGraph<G>> + 'mix>>,
}

impl<G> Default for MutatorMix<'_, G>
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    #[inline]
    fn default() -> Self {
        Self::new()
    }
}

impl<'mix, G> MutatorMix<'mix, G>
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    /// Creates an empty mix.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        Self {
            mutators: Vec::new(),
            predicates: Vec::new(),
        }
    }

    /// Creates a mix populated with the eight built-in mutators (each at
    /// weight `1.0`) and the matching eight built-in predicates.
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(feature = "smiles-support")]
    /// # fn run() {
    /// use finge_rs::{
    ///     CountEcfpFingerprint, Fingerprint, MutatorMix, SmilesRdkitGraph,
    ///     SmilesRdkitScratch,
    /// };
    /// use rand_chacha::{ChaCha8Rng, rand_core::SeedableRng};
    ///
    /// let smiles: smiles_parser::smiles::Smiles =
    ///     "c1ccccc1O".parse().expect("phenol is valid SMILES");
    /// let mut scratch = SmilesRdkitScratch::default();
    /// let inner = scratch.prepare(&smiles);
    ///
    /// let mix = MutatorMix::<SmilesRdkitGraph<'_>>::with_default_mutators_and_predicates();
    /// let mut rng = ChaCha8Rng::seed_from_u64(0xC0FFEE);
    /// let (negative, label) = mix
    ///     .sample(inner, &mut rng)
    ///     .expect("phenol accepts at least one of the eight mutators");
    ///
    /// assert!(!label.is_empty());
    /// let _fingerprint = CountEcfpFingerprint::new(2, 2048).compute(&negative);
    /// # }
    /// # #[cfg(not(feature = "smiles-support"))]
    /// # fn run() {}
    /// # run();
    /// ```
    #[inline]
    #[must_use]
    pub fn with_default_mutators_and_predicates() -> Self {
        Self::new()
            .add_mutator(Box::new(ImpossibleAtomicNumberMutator), 1.0)
            .add_mutator(Box::new(HypervalentMutator), 1.0)
            .add_mutator(Box::new(ImpossibleHCountMutator), 1.0)
            .add_mutator(Box::new(ImpossibleChargeMutator), 1.0)
            .add_mutator(Box::new(ImpossibleIsotopeMutator), 1.0)
            .add_mutator(Box::new(ImpossibleRingFlagMutator), 1.0)
            .add_mutator(Box::new(ImpossibleBondTypeMutator), 1.0)
            .add_mutator(Box::new(TopologicalPathologyMutator), 1.0)
            .add_predicate(Box::new(ImpossibleAtomicNumberPredicate))
            .add_predicate(Box::new(HypervalentPredicate))
            .add_predicate(Box::new(ImpossibleHCountPredicate))
            .add_predicate(Box::new(ImpossibleChargePredicate))
            .add_predicate(Box::new(ImpossibleIsotopePredicate))
            .add_predicate(Box::new(ImpossibleRingFlagPredicate))
            .add_predicate(Box::new(ImpossibleBondTypePredicate))
            .add_predicate(Box::new(TopologicalPathologyPredicate))
    }

    /// Appends a mutator with the given non-negative `weight`.
    ///
    /// Non-positive weights are clamped to `0.0` (the mutator stays in the
    /// list but is never picked unless every other entry also has zero
    /// weight).
    #[inline]
    #[must_use]
    pub fn add_mutator(mut self, mutator: Box<dyn Mutator<G> + 'mix>, weight: f32) -> Self {
        self.mutators.push((mutator, weight.max(0.0)));
        self
    }

    /// Appends a predicate to the labelling set.
    #[inline]
    #[must_use]
    pub fn add_predicate(
        mut self,
        predicate: Box<dyn ViolationPredicate<InvalidatedGraph<G>> + 'mix>,
    ) -> Self {
        self.predicates.push(predicate);
        self
    }

    /// Returns the number of registered mutators.
    #[inline]
    #[must_use]
    pub fn mutator_count(&self) -> usize {
        self.mutators.len()
    }

    /// Returns the number of registered predicates.
    #[inline]
    #[must_use]
    pub fn predicate_count(&self) -> usize {
        self.predicates.len()
    }

    /// Draws one labelled negative from the mix.
    ///
    /// # Errors
    ///
    /// Returns [`MutatorError`] from the chosen mutator. Common reasons:
    /// [`MutatorError::GraphTooSmall`] for empty input;
    /// [`MutatorError::NoEligibleAtom`] / [`MutatorError::NoEligibleBond`]
    /// when the mutator can't find a target (e.g. the topological-pathology
    /// mutator on an aliphatic molecule).
    pub fn sample(
        &self,
        graph: G,
        rng: &mut dyn RngCore,
    ) -> Result<(InvalidatedGraph<G>, ViolationLabel), MutatorError> {
        if self.mutators.is_empty() {
            return Err(MutatorError::NoEligibleAtom);
        }
        let index = self.pick_mutator_index(rng);
        let mutator = self
            .mutators
            .get(index)
            .map(|(m, _)| m.as_ref())
            .ok_or(MutatorError::NoEligibleAtom)?;
        let mut wrapper = InvalidatedGraph::new(graph);
        mutator.mutate_in_place(&mut wrapper, rng)?;

        let mut label = ViolationLabel::empty();
        for predicate in &self.predicates {
            if predicate.check(&wrapper) {
                label.set(predicate.class());
            }
        }
        Ok((wrapper, label))
    }

    /// Weighted random index in `0..mutators.len()`. Falls back to a uniform
    /// random pick when every weight is zero.
    fn pick_mutator_index(&self, rng: &mut dyn RngCore) -> usize {
        let total: f32 = self.mutators.iter().map(|(_, w)| *w).sum();
        if total <= 0.0 || !total.is_finite() {
            return (rng.next_u32() as usize) % self.mutators.len();
        }
        // Uniform `[0, total)` from a 24-bit random float (matches f32
        // mantissa width — avoids losing precision near `total`).
        let r = ((rng.next_u32() >> 8) as f32) * (total / (1u32 << 24) as f32);
        let mut cumulative = 0.0;
        for (idx, (_, weight)) in self.mutators.iter().enumerate() {
            cumulative += *weight;
            if r < cumulative {
                return idx;
            }
        }
        self.mutators.len() - 1
    }
}

#[cfg(test)]
mod tests {
    use rand_chacha::{ChaCha8Rng, rand_core::SeedableRng};

    use super::MutatorMix;
    use crate::{
        MutatorError, ViolationClass, ViolationLabel,
        smiles_support_impl::{SmilesRdkitGraph, SmilesRdkitScratch},
    };

    fn prepared(smiles: &str) -> (SmilesRdkitScratch, smiles_parser::smiles::Smiles) {
        let parsed: smiles_parser::smiles::Smiles =
            smiles.parse().expect("fixture SMILES should parse");
        (SmilesRdkitScratch::default(), parsed)
    }

    #[test]
    fn default_mix_has_eight_mutators_and_eight_predicates() {
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::with_default_mutators_and_predicates();
        assert_eq!(mix.mutator_count(), 8);
        assert_eq!(mix.predicate_count(), 8);
    }

    #[test]
    fn empty_mix_sample_returns_err() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::new();
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        assert!(mix.sample(inner, &mut rng).is_err());
    }

    #[test]
    fn sampled_label_always_has_at_least_one_bit_set() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::with_default_mutators_and_predicates();

        for seed in 0..32u64 {
            let mut rng = ChaCha8Rng::seed_from_u64(seed);
            if let Ok((_, label)) = mix.sample(inner, &mut rng) {
                assert!(
                    !label.is_empty(),
                    "seed {seed} produced an empty ViolationLabel"
                );
            }
        }
    }

    #[test]
    fn samples_from_ethanol_cover_at_least_seven_classes() {
        // Integration sanity check from the plan: 100 samples from "CCO"
        // with the default mix cover the seven classes whose mutators can
        // operate on aliphatic input. `TopologicalPathology` requires an
        // aromatic bond, so its bit does not appear here — see the next
        // test for the aromatic case.
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::with_default_mutators_and_predicates();
        let mut rng = ChaCha8Rng::seed_from_u64(0xBEEF_CAFE);

        let mut union = ViolationLabel::empty();
        let mut successes = 0;
        for _ in 0..100 {
            if let Ok((_, label)) = mix.sample(inner, &mut rng) {
                union |= label;
                successes += 1;
            }
        }
        assert!(successes >= 50, "too few successes: {successes}");
        for class in [
            ViolationClass::ImpossibleAtomicNumber,
            ViolationClass::Hypervalent,
            ViolationClass::ImpossibleHCount,
            ViolationClass::ImpossibleCharge,
            ViolationClass::ImpossibleIsotope,
            ViolationClass::ImpossibleRingFlag,
            ViolationClass::ImpossibleBondType,
        ] {
            assert!(
                union.is_set(class),
                "class {class:?} never appeared in 100 samples"
            );
        }
    }

    #[test]
    fn topology_class_fires_on_aromatic_input() {
        let (mut scratch, parsed) = prepared("c1ccccc1");
        let inner = scratch.prepare(&parsed);
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::with_default_mutators_and_predicates();
        let mut rng = ChaCha8Rng::seed_from_u64(0x1234);

        let mut union = ViolationLabel::empty();
        for _ in 0..100 {
            if let Ok((_, label)) = mix.sample(inner, &mut rng) {
                union |= label;
            }
        }
        assert!(
            union.is_set(ViolationClass::TopologicalPathology),
            "TopologicalPathology never fired on benzene over 100 samples",
        );
    }

    #[test]
    fn zero_weight_mutator_is_never_picked() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::new()
            .add_mutator(
                alloc::boxed::Box::new(crate::ImpossibleAtomicNumberMutator),
                1.0,
            )
            .add_mutator(
                alloc::boxed::Box::new(crate::TopologicalPathologyMutator),
                0.0,
            );
        for seed in 0..64u64 {
            let mut rng = ChaCha8Rng::seed_from_u64(seed);
            let result = mix.sample(inner, &mut rng);
            assert!(
                matches!(result, Ok(_) | Err(MutatorError::GraphTooSmall)),
                "seed {seed} unexpectedly hit the zero-weight mutator: {result:?}",
            );
        }
    }

    #[test]
    fn default_constructs_empty_mix() {
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::default();
        assert_eq!(mix.mutator_count(), 0);
        assert_eq!(mix.predicate_count(), 0);
    }

    #[test]
    fn pick_mutator_index_falls_back_when_all_weights_zero() {
        // Both mutators registered with weight 0.0 — the cumulative-sum
        // path is unreachable (total == 0.0), so `pick_mutator_index`
        // takes the uniform-fallback branch. The atomic-number mutator
        // always succeeds on CCO, so a successful sample proves the
        // fallback path was taken.
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::new()
            .add_mutator(
                alloc::boxed::Box::new(crate::ImpossibleAtomicNumberMutator),
                0.0,
            )
            .add_mutator(alloc::boxed::Box::new(crate::HypervalentMutator), 0.0);
        let mut rng = ChaCha8Rng::seed_from_u64(0xA5A5);
        let result = mix.sample(inner, &mut rng);
        assert!(
            result.is_ok(),
            "all-zero-weight fallback must still pick a mutator: {result:?}",
        );
    }
}
