//! # Interactions with fixed size integers
use std::convert::From;
use std::ops::{Add, AddAssign, Div, Mul};

use num::FromPrimitive;
use num::One;
use num::Zero;

use super::Big;

macro_rules! define_interations {
    ($t:ident) => {
        mod $t {
            use super::*;

            mod creation {
                use super::*;

                impl From<$t> for Big {
                    fn from(value: $t) -> Self {
                        Self(num::BigRational::new(value.into(), One::one()))
                    }
                }

                impl From<&$t> for Big {
                    fn from(value: &$t) -> Self {
                        Self::from(*value)
                    }
                }
            }

            mod compare {
                use super::*;

                impl PartialEq<$t> for Big {
                    fn eq(&self, other: &$t) -> bool {
                            // TODO(PERFORMANCE): Make sure that this is just as efficient as a native algorithm.
                        self.0.denom().is_one() && self.0.numer() == &num::BigInt::from(*other)
                    }
                }
            }

            mod field {
                use super::*;

                mod add {
                    use super::*;

                    impl Add<&$t> for Big {
                        type Output = Self;

                        fn add(self, rhs: &$t) -> Self::Output {
                            // TODO(PERFORMANCE): Make sure that this is just as efficient as a native algorithm.
                            let (mut numer, denom) = self.0.into();
                            numer += &denom * rhs;
                            Self((numer, denom).into())
                        }
                    }

                    impl Add<Option<&$t>> for Big {
                        type Output = Self;

                        fn add(self, rhs: Option<&$t>) -> Self::Output {
                            // TODO(PERFORMANCE): Make sure that this is just as efficient as a native algorithm.
                            match rhs {
                                None => self,
                                Some(rhs) => Add::add(self, rhs),
                            }
                        }
                    }

                    impl Add<&$t> for &Big {
                        type Output = Big;

                        fn add(self, rhs: &$t) -> Self::Output {
                            // TODO(PERFORMANCE): Make sure that this is just as efficient as a native algorithm.
                            self.clone().add(rhs)
                        }
                    }

                    impl Add<Option<&$t>> for &Big {
                        type Output = Big;

                        fn add(self, rhs: Option<&$t>) -> Self::Output {
                            // TODO(PERFORMANCE): Make sure that this is just as efficient as a native algorithm.
                            let copy = self.clone();
                            match rhs {
                                None => copy,
                                Some(rhs) => Add::add(copy, rhs),
                            }
                        }
                    }

                    impl AddAssign<$t> for Big {
                        fn add_assign(&mut self, rhs: $t) {
                            // TODO(PERFORMANCE): Make sure that this is just as efficient as a native algorithm.
                            let uptyped = num::BigRational::new(rhs.into(), One::one());
                            self.0 += uptyped;
                        }
                    }

                    impl AddAssign<&$t> for Big {
                        fn add_assign(&mut self, rhs: &$t) {
                            // TODO(PERFORMANCE): Make sure that this is just as efficient as a native algorithm.
                            self.add_assign(*rhs)
                        }
                    }
                }

                mod mul {
                    use super::*;

                    impl Mul<&$t> for Big {
                        type Output = Big;

                        fn mul(self, rhs: &$t) -> Self::Output {
                            // TODO(PERFORMANCE): Make sure that this is just as efficient as a native algorithm.
                            let (mut numer, denom) = self.0.into();
                            numer *= *rhs;
                            Self((numer, denom).into())
                        }
                    }

                    impl Mul<&$t> for &Big {
                        type Output = Big;

                        fn mul(self, rhs: &$t) -> Self::Output {
                            // TODO(PERFORMANCE): Make sure that this is just as efficient as a native algorithm.
                            let (mut numer, denom) = (self.0.numer().clone(), self.0.denom().clone());
                            numer *= *rhs;
                            Big((numer, denom).into())
                        }
                    }

                    impl Mul<Option<&$t>> for Big {
                        type Output = Big;

                        fn mul(self, rhs: Option<&$t>) -> Self::Output {
                            // TODO(PERFORMANCE): Make sure that this is just as efficient as a native algorithm.
                            match rhs {
                                None => Big::zero(),
                                Some(rhs) => Mul::mul(self, rhs),
                            }
                        }
                    }

                    impl Mul<Option<&$t>> for &Big {
                        type Output = Big;

                        fn mul(self, rhs: Option<&$t>) -> Self::Output {
                            // TODO(PERFORMANCE): Make sure that this is just as efficient as a native algorithm.
                            match rhs {
                                None => Big::zero(),
                                Some(rhs) => Mul::mul(self, rhs),
                            }
                        }
                    }
                }

                mod div {
                    use super::*;

                    impl Div<&$t> for Big {
                        type Output = Big;

                        fn div(self, rhs: &$t) -> Self::Output {
                            // TODO(PERFORMANCE): Make sure that this is just as efficient as a native algorithm.
                            let (numer, mut denom) = self.0.into();
                            denom *= *rhs;
                            Self((numer, denom).into())
                        }
                    }
                }
            }
        }
    }
}

define_interations!(i32);
define_interations!(i64);
define_interations!(i128);
define_interations!(isize);
define_interations!(u32);
define_interations!(u64);
define_interations!(u128);
define_interations!(usize);
