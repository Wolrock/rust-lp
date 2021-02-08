use crate::data::linear_algebra::traits::Element;
use crate::data::number_types::traits::Field;

pub trait Rhs =
    Field +
    Element +
;
