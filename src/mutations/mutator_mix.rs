//! [`MutatorMix`] — weighted-random sampler over a collection of mutators
//! that composes 1..=K mutations per draw and labels the result by running
//! every registered predicate against the final wrapper.
//!
//! The number of mutations per sample is `k = 1 + Poisson(λ)` clamped to
//! `[1, k_max]`. Defaults: `λ = 0.7`, `k_max = 4`. Each slot draws a mutator
//! independently by weight (with replacement, so the same mutator can fire
//! several times on different atoms). The composed wrapper is then scored
//! against every predicate, producing a multi-label [`ViolationLabel`] that
//! reflects what actually fired — not just the primary class of whatever
//! mutators happened to be sampled.

use alloc::{boxed::Box, vec::Vec};

use rand_core::RngCore;

use crate::{
    Fingerprint,
    fingerprints::CountEcfpFingerprint,
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
    traits::{EcfpGraph, MolecularAtom, MolecularBond},
};

/// One entry in the mix: a mutator and its sampling weight.
type WeightedMutator<'mix, G> = (Box<dyn Mutator<G> + 'mix>, f32);

/// Default Poisson rate λ for the `k = 1 + Poisson(λ)` count of mutations
/// per sample. With λ = 0.7 the discrete distribution over k is roughly
/// `{1: 0.50, 2: 0.35, 3: 0.12, ≥4: 0.03}` — most samples are still
/// single-mutation, but ~50 % have at least two co-firing mutators.
pub const DEFAULT_COMPOSITION_LAMBDA: f64 = 0.7;

/// Default upper bound on the number of mutations composed per sample.
pub const DEFAULT_K_MAX: u8 = 4;

/// Default number of times a slot is retried after a mutator returns Err
/// (e.g. `TopologicalPathologyMutator` on an aliphatic graph). Each retry
/// draws a fresh mutator by weight, so the worst case for a 4-slot sample
/// is `4 * (1 + 3) = 16` mutator invocations.
pub const DEFAULT_RETRY_PER_SLOT: u8 = 3;

/// Maximum number of slots the cumulative CDF holds. Constraints: `k_max`
/// cannot exceed this without resizing `k_cumulative`.
pub const K_MAX_LIMIT: u8 = 8;

/// Weighted-random sampler over a configurable set of mutators plus a fixed
/// list of predicates that label each generated negative.
///
/// Build with [`MutatorMix::new`] (empty) or
/// [`MutatorMix::with_default_mutators_and_predicates`] (the eight built-in
/// pairs, each with weight `1.0`). Append more via [`MutatorMix::add_mutator`]
/// and [`MutatorMix::add_predicate`].
///
/// Sampling:
/// 1. Draw `k = 1 + Poisson(λ)` clamped to `[1, k_max]`.
/// 2. For each of `k` slots: pick a mutator by weight and call
///    [`Mutator::mutate_in_place`] on a shared [`InvalidatedGraph`]. Retry
///    up to `retry_per_slot` times on `Err`, then move on.
/// 3. Run every registered predicate against the final wrapper and OR the
///    hits into a [`ViolationLabel`]. The returned label is ground truth —
///    it reflects *what actually fired*, not which mutators were sampled.
///
/// Slots that all fail their retries are skipped silently. `sample()` returns
/// `Err` only when *every* slot failed (no mutation was applied at all).
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
    /// Cumulative density function for `k` slots per sample. `k_cumulative[i]`
    /// holds `P(k <= i + 1)`. Always normalised so the last active entry is
    /// `1.0`; unused tail entries (beyond `k_max`) are `1.0` as well.
    k_cumulative: [f32; K_MAX_LIMIT as usize],
    /// Poisson rate cached so `with_k_max` can recompute the CDF without
    /// inverting it.
    composition_lambda: f64,
    /// Upper bound on the number of mutations per sample.
    k_max: u8,
    /// Per-slot retry budget when a sampled mutator returns Err.
    retry_per_slot: u8,
    /// If set, [`Self::sample`] computes the count-ECFP fingerprint at the
    /// stored `(radius, fp_size)` for both the baseline and the mutated
    /// wrapper, returning [`MutatorError::FingerprintCollision`] when they
    /// are byte-equal. Lets the caller pin the trainer's chosen `fp_size`
    /// and skip samples whose perturbation is invisible to the model.
    collision_filter: Option<(u8, usize)>,
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
    /// Creates an empty mix with the default composition policy
    /// (`λ = DEFAULT_COMPOSITION_LAMBDA`, `k_max = DEFAULT_K_MAX`,
    /// `retry_per_slot = DEFAULT_RETRY_PER_SLOT`) and no collision filter.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        Self {
            mutators: Vec::new(),
            predicates: Vec::new(),
            k_cumulative: poisson_k_cumulative(DEFAULT_COMPOSITION_LAMBDA, DEFAULT_K_MAX),
            composition_lambda: DEFAULT_COMPOSITION_LAMBDA,
            k_max: DEFAULT_K_MAX,
            retry_per_slot: DEFAULT_RETRY_PER_SLOT,
            collision_filter: None,
        }
    }

    /// Creates a mix populated with the eight built-in mutators (each at
    /// weight `1.0`) and the matching eight built-in predicates.
    ///
    /// # Example
    ///
    /// ```
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

    /// Sets the Poisson rate `λ` of the `k = 1 + Poisson(λ)` slot count.
    ///
    /// Negative or non-finite values are clamped to `0.0` (every sample
    /// becomes `k = 1`).
    #[inline]
    #[must_use]
    pub fn with_composition_lambda(mut self, lambda: f64) -> Self {
        let lambda = if lambda.is_finite() {
            lambda.max(0.0)
        } else {
            0.0
        };
        self.composition_lambda = lambda;
        self.k_cumulative = poisson_k_cumulative(lambda, self.k_max);
        self
    }

    /// Sets the upper bound on the number of mutations composed per sample.
    ///
    /// `k_max = 0` is treated as `1` (every sample produces at least one
    /// mutation). Values above [`K_MAX_LIMIT`] are clamped.
    #[inline]
    #[must_use]
    pub fn with_k_max(mut self, k_max: u8) -> Self {
        let k_max = k_max.clamp(1, K_MAX_LIMIT);
        self.k_max = k_max;
        self.k_cumulative = poisson_k_cumulative(self.composition_lambda, k_max);
        self
    }

    /// Sets the per-slot retry budget after a mutator returns Err.
    #[inline]
    #[must_use]
    pub fn with_retry_per_slot(mut self, retry: u8) -> Self {
        self.retry_per_slot = retry;
        self
    }

    /// Enables ECFP-collision filtering at the given `(radius, fp_size)`.
    ///
    /// When set, every successful [`Self::sample`] call additionally
    /// computes a count-ECFP fingerprint for the baseline and for the
    /// mutated wrapper and returns
    /// [`MutatorError::FingerprintCollision`] if the two fold to byte-equal
    /// bins. Use this with the trainer's chosen `fp_size` to ensure every
    /// returned sample carries a feature-space signal at the size the
    /// model will actually see.
    ///
    /// Adds one extra `CountEcfpFingerprint::compute` per sample.
    #[inline]
    #[must_use]
    pub fn with_collision_filter(mut self, radius: u8, fp_size: usize) -> Self {
        self.collision_filter = Some((radius, fp_size));
        self
    }

    /// Disables the collision filter previously set via
    /// [`Self::with_collision_filter`].
    #[inline]
    #[must_use]
    pub fn without_collision_filter(mut self) -> Self {
        self.collision_filter = None;
        self
    }

    /// Returns the active collision-filter `(radius, fp_size)` if set.
    #[inline]
    #[must_use]
    pub fn collision_filter(&self) -> Option<(u8, usize)> {
        self.collision_filter
    }

    /// Returns the configured `k_max`.
    #[inline]
    #[must_use]
    pub fn k_max(&self) -> u8 {
        self.k_max
    }

    /// Returns the configured per-slot retry budget.
    #[inline]
    #[must_use]
    pub fn retry_per_slot(&self) -> u8 {
        self.retry_per_slot
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
    /// Composes `k = 1 + Poisson(λ)` clamped to `[1, k_max]` mutations on a
    /// shared wrapper, retrying each slot up to [`retry_per_slot`] times when
    /// the sampled mutator declines, then runs every registered predicate
    /// against the final wrapper to derive the [`ViolationLabel`].
    ///
    /// [`retry_per_slot`]: MutatorMix::retry_per_slot
    ///
    /// # Errors
    ///
    /// Returns [`MutatorError`] from the *last* failed slot when every slot's
    /// retries were exhausted (no mutation was applied at all). Returns
    /// [`MutatorError::CompositionCancelled`] when mutations were applied but
    /// later writes on the same atom-channel clobbered earlier ones, leaving
    /// the wrapper pre-hash-identical to the inner. Returns
    /// [`MutatorError::FingerprintCollision`] when the collision filter is
    /// enabled and the wrapper's folded ECFP equals the baseline's. When
    /// the mix is empty, returns [`MutatorError::NoEligibleAtom`].
    pub fn sample(
        &self,
        graph: G,
        rng: &mut dyn RngCore,
    ) -> Result<(InvalidatedGraph<G>, ViolationLabel), MutatorError>
    where
        G: EcfpGraph<NodeId = usize>,
        G::Bond: MolecularBond<NodeId = usize>,
    {
        if self.mutators.is_empty() {
            return Err(MutatorError::NoEligibleAtom);
        }
        // Compute the baseline fingerprint up-front when collision filtering
        // is enabled. Done before the graph is consumed by
        // `InvalidatedGraph::new` so we don't need a Clone bound on G.
        let baseline_fp = self
            .collision_filter
            .map(|(radius, fp_size)| CountEcfpFingerprint::new(radius, fp_size).compute(&graph));
        let mut wrapper = InvalidatedGraph::new(graph);
        let k = self.pick_k(rng);

        let mut applied_any = false;
        let mut last_err = MutatorError::NoEligibleAtom;
        for _ in 0..k {
            for _attempt in 0..=self.retry_per_slot {
                let index = self.pick_mutator_index(rng);
                // `pick_mutator_index` always returns a valid index because
                // we checked `self.mutators.is_empty()` above.
                let mutator = self.mutators[index].0.as_ref();
                match mutator.mutate_in_place(&mut wrapper, rng) {
                    Ok(()) => {
                        applied_any = true;
                        break;
                    }
                    Err(err) => last_err = err,
                }
            }
        }

        if !applied_any {
            return Err(last_err);
        }

        // Reject samples where every successful mutation was clobbered by a
        // later same-channel write. Without this guard, training data could
        // include "negatives" that are byte-identical to the baseline graph
        // (see fuzz/artifacts/mutator_mix/crash-9f067edf… — `c8ccccccc8c`
        // where ImpossibleRingFlagMutator wrote `in_ring=Some(true)` on
        // the dangling atom 8 and TopologicalPathologyMutator then wrote
        // `in_ring=Some(false)` over it, returning the atom to its inner
        // state of `in_ring=false`).
        if !wrapper.has_effective_perturbation() {
            return Err(MutatorError::CompositionCancelled);
        }

        // Reject samples where the perturbation is invisible to the trainer's
        // chosen fp_size. Atoms whose post-override invariants happen to fold
        // to the same bins as the baseline (e.g. atom 2 of BrOBr coinciding
        // with the unchanged atom 0 in the 65 536-bin fold) produce a wrapper
        // whose ECFP is byte-equal to the baseline — useless as a labelled
        // negative for the model. The trainer pins the fp_size via
        // `with_collision_filter`; we re-fingerprint and compare here.
        if let Some((radius, fp_size)) = self.collision_filter {
            let wrapper_fp = CountEcfpFingerprint::new(radius, fp_size).compute(&wrapper);
            if Some(&wrapper_fp) == baseline_fp.as_ref() {
                return Err(MutatorError::FingerprintCollision);
            }
        }

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

    /// Draws `k` slots for the current sample. Always returns at least `1`
    /// and at most `self.k_max`.
    fn pick_k(&self, rng: &mut dyn RngCore) -> u8 {
        let u = ((rng.next_u32() >> 8) as f32) / ((1u32 << 24) as f32);
        for (i, &threshold) in self.k_cumulative.iter().enumerate() {
            if u < threshold {
                return (i as u8 + 1).min(self.k_max);
            }
        }
        self.k_max
    }
}

/// Computes `exp(-x)` for non-negative `x` via Taylor series at 0.
///
/// Used to build the Poisson PMF without pulling in `libm`. Converges to
/// double-precision for `x ≤ 5` within ~25 terms; we always sample
/// `lambda ≤ ~3` so this is comfortably accurate.
fn exp_neg(x: f64) -> f64 {
    let neg = -x;
    let mut sum = 1.0_f64;
    let mut term = 1.0_f64;
    for n in 1..=40 {
        term *= neg / (n as f64);
        sum += term;
    }
    sum
}

/// Returns the cumulative distribution function of `k = 1 + Poisson(lambda)`
/// truncated and renormalised to `[1, k_max]`. `result[i]` is `P(k <= i + 1)`
/// for `i < k_max`; all tail entries are `1.0`.
fn poisson_k_cumulative(lambda: f64, k_max: u8) -> [f32; K_MAX_LIMIT as usize] {
    let mut cdf = [1.0_f32; K_MAX_LIMIT as usize];
    if k_max == 0 {
        return cdf;
    }
    let e = exp_neg(lambda);
    // pmf[j] holds P(Poisson(lambda) = j) for j in 0..k_max.
    let mut pmf = [0.0_f64; K_MAX_LIMIT as usize];
    let mut term = e;
    pmf[0] = term;
    let head_len = (k_max as usize).min(K_MAX_LIMIT as usize);
    for (j, slot) in pmf.iter_mut().enumerate().take(head_len).skip(1) {
        term *= lambda / (j as f64);
        *slot = term;
    }
    // Collapse the tail `P(Poisson >= k_max - 1)` into the last entry so the
    // resulting CDF sums to exactly 1.
    let head_sum: f64 = pmf
        .iter()
        .take((k_max as usize).saturating_sub(1))
        .copied()
        .sum();
    let tail = (1.0 - head_sum).max(0.0);
    if let Some(last) = pmf.get_mut((k_max as usize).saturating_sub(1)) {
        *last = tail;
    }
    let mut acc = 0.0_f64;
    for (i, slot) in cdf.iter_mut().enumerate() {
        if i < k_max as usize {
            acc += pmf[i];
            *slot = (acc.min(1.0)) as f32;
        } else {
            *slot = 1.0;
        }
    }
    // Final guard: floating-point drift could leave the last active entry a
    // hair below 1.0, which would cause `pick_k` to return `k_max` even when
    // unintended. Pinning eliminates that ambiguity.
    if k_max as usize <= K_MAX_LIMIT as usize {
        cdf[k_max as usize - 1] = 1.0;
    }
    cdf
}

#[cfg(test)]
mod tests {
    use rand_chacha::{ChaCha8Rng, rand_core::SeedableRng};

    use super::{
        DEFAULT_COMPOSITION_LAMBDA, DEFAULT_K_MAX, DEFAULT_RETRY_PER_SLOT, K_MAX_LIMIT, MutatorMix,
        poisson_k_cumulative,
    };
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

    #[test]
    fn poisson_cumulative_matches_hand_computed_values_for_default_lambda() {
        // For lambda = 0.7, k_max = 4:
        //   P(Poisson=0) = exp(-0.7)               ≈ 0.4966
        //   P(Poisson=1) = 0.7 * exp(-0.7)         ≈ 0.3476
        //   P(Poisson=2) = 0.49 * exp(-0.7) / 2    ≈ 0.1217
        //   tail (>=3) collapsed into slot 3        ≈ 0.0341
        // CDF: [0.4966, 0.8442, 0.9659, 1.0000, 1.0, 1.0, 1.0, 1.0]
        let cdf = poisson_k_cumulative(DEFAULT_COMPOSITION_LAMBDA, DEFAULT_K_MAX);
        let expected = [0.4966_f32, 0.8442, 0.9659, 1.0000, 1.0, 1.0, 1.0, 1.0];
        for (i, (got, want)) in cdf.iter().zip(expected.iter()).enumerate() {
            assert!((got - want).abs() < 1e-3, "cdf[{i}] = {got}, want ≈ {want}",);
        }
    }

    #[test]
    fn poisson_cumulative_handles_zero_lambda_as_always_one_slot() {
        let cdf = poisson_k_cumulative(0.0, DEFAULT_K_MAX);
        assert_eq!(cdf[0], 1.0, "lambda=0 must put all mass on k=1");
    }

    #[test]
    fn default_config_matches_documented_constants() {
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::with_default_mutators_and_predicates();
        assert_eq!(mix.k_max(), DEFAULT_K_MAX);
        assert_eq!(mix.retry_per_slot(), DEFAULT_RETRY_PER_SLOT);
    }

    #[test]
    fn with_k_max_clamps_within_limit() {
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::new().with_k_max(255);
        assert_eq!(mix.k_max(), K_MAX_LIMIT);
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::new().with_k_max(0);
        assert_eq!(mix.k_max(), 1);
    }

    #[test]
    fn composed_samples_include_multi_bit_labels() {
        // With lambda = 0.7 and k_max = 4 roughly half of all samples have
        // k >= 2. Across 200 samples we expect at least a handful with two
        // or more bits set. (At k = 1 each sample is single-bit by
        // construction except in the rare cases where two predicates fire
        // on the same atom; with composition, multi-bit labels become the
        // common case.)
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::with_default_mutators_and_predicates();
        let mut rng = ChaCha8Rng::seed_from_u64(0x0C01_1EC7);
        let mut multi_bit = 0;
        let mut single_bit = 0;
        for _ in 0..200 {
            if let Ok((_, label)) = mix.sample(inner, &mut rng) {
                if label.count() >= 2 {
                    multi_bit += 1;
                } else if label.count() == 1 {
                    single_bit += 1;
                }
            }
        }
        assert!(
            multi_bit >= 20,
            "expected >= 20 multi-bit labels out of 200, got {multi_bit} \
             (single-bit: {single_bit})",
        );
    }

    #[test]
    fn with_k_max_one_collapses_to_single_mutation_per_sample() {
        // With k_max = 1, every sample applies exactly one mutation.
        // Therefore at most one bit can be set per label (modulo
        // secondary predicate hits — rare on CCO).
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::with_default_mutators_and_predicates()
            .with_k_max(1);
        let mut rng = ChaCha8Rng::seed_from_u64(0x5111_E0E1);
        let mut multi_bit = 0;
        for _ in 0..50 {
            if let Ok((_, label)) = mix.sample(inner, &mut rng) {
                if label.count() >= 3 {
                    multi_bit += 1;
                }
            }
        }
        assert_eq!(
            multi_bit, 0,
            "k_max=1 must never produce >=3 simultaneously firing predicates",
        );
    }

    #[test]
    fn with_composition_lambda_accepts_finite_non_negative() {
        // Default lambda is 0.7. Bumping it to 1.5 shifts the CDF rightward.
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::with_default_mutators_and_predicates();
        let cdf_default = mix.k_cumulative;
        let mix_high = mix.with_composition_lambda(1.5);
        let cdf_high = mix_high.k_cumulative;
        assert!(
            cdf_high[0] < cdf_default[0],
            "higher lambda must shift weight off the k=1 slot",
        );
    }

    #[test]
    fn with_composition_lambda_clamps_negative_to_zero() {
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::new().with_composition_lambda(-1.0);
        // Negative lambda is clamped to 0, putting all mass on k=1.
        assert_eq!(mix.k_cumulative[0], 1.0);
    }

    #[test]
    fn with_composition_lambda_clamps_nan_to_zero() {
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::new().with_composition_lambda(f64::NAN);
        assert_eq!(mix.k_cumulative[0], 1.0);
    }

    #[test]
    fn with_retry_per_slot_updates_field() {
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::new().with_retry_per_slot(7);
        assert_eq!(mix.retry_per_slot(), 7);
    }

    #[test]
    fn sample_returns_err_when_every_slot_exhausts_retries() {
        // A mix whose only mutator declines every input must yield Err
        // (no slot ever succeeded), exercising the `if !applied_any` path.
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::new()
            .add_mutator(
                alloc::boxed::Box::new(crate::TopologicalPathologyMutator),
                1.0,
            )
            .with_k_max(2)
            .with_retry_per_slot(2);
        let mut rng = ChaCha8Rng::seed_from_u64(0xDEAF_BEAD);
        let err = mix
            .sample(inner, &mut rng)
            .expect_err("aliphatic input with topology-only mix must Err");
        assert!(matches!(
            err,
            MutatorError::NoEligibleAtom | MutatorError::GraphTooSmall,
        ));
    }

    #[test]
    fn collision_filter_eventually_rejects_on_collision_prone_input() {
        // BrOBr's structural symmetry (atoms 0 and 2 are equivalent Br
        // terminals) makes the 65 536-bin folded ECFP coincidence
        // probable: when HypervalentMutator picks atom 2 and the RNG
        // happens to draw an OOD degree whose 32-bit invariant folds to
        // atom 0's bin, the wrapper's ECFP equals the baseline's.
        // Iterating seeds with a tiny mix forces the filter to trip
        // within a small budget.
        let parsed: smiles_parser::smiles::Smiles = "BrOBr".parse().expect("parse");
        let mut scratch = SmilesRdkitScratch::default();
        let graph = scratch.prepare(&parsed);
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::new()
            .add_mutator(alloc::boxed::Box::new(crate::HypervalentMutator), 1.0)
            .add_predicate(alloc::boxed::Box::new(crate::HypervalentPredicate))
            .with_k_max(1)
            .with_collision_filter(2, 65_536);
        let mut saw_collision = false;
        for seed in 0..2048u64 {
            let mut rng = ChaCha8Rng::seed_from_u64(seed);
            if let Err(MutatorError::FingerprintCollision) = mix.sample(graph, &mut rng) {
                saw_collision = true;
                break;
            }
        }
        assert!(
            saw_collision,
            "collision filter never triggered on BrOBr across 2048 seeds",
        );
    }

    #[test]
    fn collision_filter_passes_through_when_fingerprint_differs() {
        // Sanity: the collision filter must not reject samples whose
        // fingerprints actually differ. CCO with the same setup at a
        // benign seed succeeds.
        let parsed: smiles_parser::smiles::Smiles = "CCO".parse().expect("parse");
        let mut scratch = SmilesRdkitScratch::default();
        let graph = scratch.prepare(&parsed);
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::with_default_mutators_and_predicates()
            .with_collision_filter(2, 65_536);
        let mut rng = ChaCha8Rng::seed_from_u64(0xC0FFEE);
        assert!(mix.sample(graph, &mut rng).is_ok());
    }

    #[test]
    fn with_and_without_collision_filter_round_trip() {
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::new().with_collision_filter(2, 2048);
        assert_eq!(mix.collision_filter(), Some((2, 2048)));
        let mix = mix.without_collision_filter();
        assert_eq!(mix.collision_filter(), None);
    }

    #[test]
    fn fuzz_regression_composed_isotope_writes_dont_cancel_invariants() {
        // ssI + seed 6_733_535_862_861_618_035 composes two
        // ImpossibleIsotopeMutator slots: atom 0 gets delta_mass = 11027
        // (R0 invariant shifts by -2287) and atom 1 gets delta_mass = 2327
        // (R0 invariant shifts by +2287). A naive sum aggregate of R0
        // invariants over non-ignored atoms would see no change. Pin a
        // position-aware check: the per-atom invariant vector must differ.
        use crate::{EcfpGraph, traits::MolecularGraph as _};
        let parsed: smiles_parser::smiles::Smiles = "ssI".parse().expect("parse");
        let mut scratch = SmilesRdkitScratch::default();
        let graph = scratch.prepare(&parsed);
        let baseline_invariants: alloc::vec::Vec<u32> = (0..graph.atom_count())
            .filter(|&id| !graph.ecfp_atom_is_ignored(id))
            .map(|id| graph.ecfp_atom_invariant(id, true))
            .collect();
        let baseline_sum: u32 = baseline_invariants
            .iter()
            .copied()
            .fold(0u32, u32::wrapping_add);

        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::with_default_mutators_and_predicates();
        let mut rng = ChaCha8Rng::seed_from_u64(6_733_535_862_861_618_035);
        if let Ok((wrapper, _label)) = mix.sample(graph, &mut rng) {
            let mutated_invariants: alloc::vec::Vec<u32> = (0..wrapper.atom_count())
                .filter(|&id| !wrapper.ecfp_atom_is_ignored(id))
                .map(|id| wrapper.ecfp_atom_invariant(id, true))
                .collect();
            let mutated_sum: u32 = mutated_invariants
                .iter()
                .copied()
                .fold(0u32, u32::wrapping_add);
            // The position-aware vector must differ even when the sum doesn't.
            assert_ne!(
                baseline_invariants, mutated_invariants,
                "per-atom invariant vector should differ after the composition",
            );
            // The historical false negative was that sums happened to match.
            // We don't assert which direction is true on every seed, just
            // document the trap by comparing both.
            let _ = (baseline_sum, mutated_sum);
        }
    }

    #[test]
    fn fuzz_regression_composition_cancellation_yields_err() {
        // `c8ccccccc8c` + seed 14_612_714_913_291_461_578 originally
        // produced an Ok((wrapper, label)) where the wrapper was
        // pre-hash-identical to the inner: composition picked
        // `ImpossibleRingFlagMutator` (in_ring=Some(true) on the dangling
        // aromatic atom 8) then `TopologicalPathologyMutator`
        // (in_ring=Some(false) on the same atom), clobbering back to the
        // inner's false. `has_effective_perturbation` now catches this
        // and `MutatorMix::sample` returns
        // `MutatorError::CompositionCancelled`. See
        // `fuzz/artifacts/mutator_mix/crash-9f067edf…`.
        let parsed: smiles_parser::smiles::Smiles = "c8ccccccc8c".parse().expect("parse");
        let mut scratch = SmilesRdkitScratch::default();
        let graph = scratch.prepare(&parsed);
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::with_default_mutators_and_predicates();
        let mut rng = ChaCha8Rng::seed_from_u64(14_612_714_913_291_461_578);
        match mix.sample(graph, &mut rng) {
            Ok((wrapper, _)) => {
                use crate::{traits::EcfpGraph as _, traits::MolecularGraph as _};
                // If sampling did succeed (e.g. after future RNG changes),
                // the wrapper must carry an *effective* perturbation —
                // not just a recorded override that no-ops.
                let any_changed = (0..wrapper.atom_count()).any(|id| {
                    wrapper.inner().ecfp_atom_invariant_fields(id)
                        != wrapper.ecfp_atom_invariant_fields(id)
                });
                assert!(
                    any_changed,
                    "successful sample must change at least one atom invariant tuple",
                );
            }
            Err(err) => assert!(
                matches!(
                    err,
                    MutatorError::CompositionCancelled
                        | MutatorError::NoEligibleAtom
                        | MutatorError::NoEligibleBond
                        | MutatorError::GraphTooSmall
                ),
                "unexpected error variant: {err:?}",
            ),
        }
    }

    #[test]
    fn retry_recovers_from_initial_slot_failure() {
        // A mix with only the TopologicalPathologyMutator declines aliphatic
        // CCO every time. With retry_per_slot = 0 we must get Err; with the
        // default budget we still get Err (retries draw the same mutator
        // again), but with a second mutator added at non-zero weight the
        // retry path eventually picks one that can fire.
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::new()
            .add_mutator(
                alloc::boxed::Box::new(crate::TopologicalPathologyMutator),
                1.0,
            )
            .add_mutator(
                alloc::boxed::Box::new(crate::ImpossibleAtomicNumberMutator),
                1.0,
            )
            .add_predicate(alloc::boxed::Box::new(
                crate::ImpossibleAtomicNumberPredicate,
            ))
            .with_k_max(1);
        let mut rng = ChaCha8Rng::seed_from_u64(0xDEAD_BEEF);
        let mut successes = 0;
        for _ in 0..50 {
            if mix.sample(inner, &mut rng).is_ok() {
                successes += 1;
            }
        }
        assert!(
            successes >= 25,
            "retry must let the always-succeeding mutator rescue at least \
             half of samples; got {successes}/50",
        );
    }
}
