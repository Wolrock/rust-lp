use crate::data::number_types::nonzero::{NonzeroSigned, Nonzero};

pub trait NonzeroFactorizable: NonzeroSigned + Nonzero {
    /// Some prime greater than 1.
    type Factor: Nonzero + Ord;
    /// How often the factor appears in the number.
    type Power: Nonzero;

    fn small_factors(&self) -> Vec<(Self::Factor, Self::Power)>;
}
