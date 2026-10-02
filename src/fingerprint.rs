/// Generic interface for computing a fingerprint from a molecular graph.
pub trait Fingerprint<G> {
    /// The produced fingerprint representation.
    type Output;

    /// Computes a fingerprint for the provided molecular graph.
    fn compute(&self, graph: &G) -> Self::Output;
}
