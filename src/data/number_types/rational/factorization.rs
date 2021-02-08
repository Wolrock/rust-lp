use std::num::NonZeroI8;
use crate::data::number_types::traits::factorization::NonzeroFactorizable;
use crate::data::number_types::rational::Rational;
use crate::data::number_types::nonzero::{NonzeroSigned, NonzeroSign};
use crate::algorithm::utilities::merge_sparse_indices;
use std::ops::Add;

/// Prime factorization representation of a nonzero rational number.
///
/// Includes a sign.
pub struct Factorization {
    /// Whether the number is negative.
    sign: NonzeroSign,
    /// `(prime factor, power)` tuples.
    ///
    /// The factors should all be smaller than 64 bits and can have negative powers; that is, appear
    /// in the denominator. The powers can't be zero, as this is a sparse representation.
    ///
    /// When this field is empty, the value `1` or `-1` is represented, depending on `sign`.
    factors: Vec<(u64, NonZeroI8)>,
}

impl<R> NonzeroFactorizable for R
where
    R: Rational<
        Numerator: NonzeroFactorizable,
        Denominator: NonzeroFactorizable,
    > + NonzeroSigned,
{
    type Factor = <R::Numerator as NonzeroFactorizable>::Factor;
    type Power = <R::Numerator as NonzeroFactorizable>::Factor;

    fn small_factors(&self) -> Vec<(Self::Factor, Self::Power)> {
        let numerator_factors = self.numerator().small_factors().into_iter();
        let denominator_factors = self.denominator().small_factors().into_iter();

        // TODO: Neg the denom factors
        merge_sparse_indices(numerator_factors, denominator_factors, Add::add).collect()
    }
}
