//! Inverted-index backend for sparse Tanimoto queries.
//!
//! A posting list per feature id accelerates the sparse scan. Landed only if
//! Phase 5 profiling shows the pairwise merge is the bottleneck at target n.
