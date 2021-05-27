use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::lower_upper_remultiply_factor::permutation::Permutation;

// TODO(Performance): Maybe generic type for col major and row major
#[derive(Eq, PartialEq, Clone, Debug)]
pub struct Rotate {
    pub start: usize,
    pub end: usize,
}

impl Rotate {
    fn new(start: usize, end: usize) -> Self {
        debug_assert!(start <= end);
        Self {
            start: start,
            end: end,
        }
    }
}

impl Permutation for Rotate {
    // Indices are exclusive
    fn forward(&self, i: &mut usize) {
        debug_assert!(*i < self.end);
        if *i < self.start {
        } else if *i < self.end - 1 {
            *i += 1;
        } else if *i < self.end {
            *i = self.start;
        }
    }

    //interpretation of end/start?
    fn backward(&self, i: &mut usize) {
        debug_assert!(*i < self.end);
        if *i < self.start {
        } else if *i > self.start {
            *i -= 1;
        } else if *i < self.start + 1 {
            *i = self.end;
        }
    }

    fn len(&self) -> usize {
        todo!()
    }

    // fn forward_sorted<T>(&self, items: &mut [(usize, T)]) {
    //     todo!();
    // }

    // fn backward_sorted<T>(&self, items: &mut [(usize, T)]) {
    //     todo!();
    // }
}

#[cfg(test)]
mod test {
    use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::lower_upper_remultiply_factor::permutation::Permutation;
    use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::lower_upper_remultiply_factor::permutation::rotate::Rotate;
    #[test]
    fn no_change() {
        let mut items = vec![(0, 0), (1, 1), (2, 2)];
        // TODO(Debug): constructor
        Rotate::new(2, 3).forward_sorted(&mut items);
        let expected = vec![(0, 0), (1, 1), (2, 2)];
        assert_eq!(items, expected);
    }

    fn one_forward() {}
}
