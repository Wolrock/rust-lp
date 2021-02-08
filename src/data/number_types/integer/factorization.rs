use std::num::{NonZeroI8, NonZeroU8};
use crate::data::number_types::traits::factorization::NonzeroFactorizable;

macro_rules! impl_factorization_non_zero {
    ($t:ident, $f:ident) => {
        impl NonzeroFactorizable for $t {
            type Factor = $f;
            type Power = u8;
            // type Power = NonZeroU8;

            fn small_factors(&self) -> Vec<(Self::Factor, Self::Power)> {
                unimplemented!()
            }
        }
    }
}

// impl_factorization_non_zero!(NonZeroI8, NonZeroU8);
impl_factorization_non_zero!(i8, u8);
