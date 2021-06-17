//! # LU decomposition
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fmt;
use std::fmt::Display;
use std::iter;

use std::thread::sleep;
use std::time::{Duration, SystemTime};

use num::Zero;

use crate::algorithm::two_phase::matrix_provider::column::{Column, OrderedColumn};
use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::lower_upper_remultiply_factor::permutation::{
    FullPermutation, Permutation, Rotate,
};
use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::BasisInverse;
use crate::algorithm::two_phase::tableau::inverse_maintenance::{ops, ColumnComputationInfo};
use crate::data::linear_algebra::traits::SparseElement;
use crate::data::linear_algebra::vector::{SparseVector, Vector};

use crate::RB;
use std::array::IntoIter;
use std::borrow::Borrow;
use std::iter::FromIterator;

mod decomposition;
mod permutation;

// time_basis_change: SystemTime = SystemTime::now();

/// Decompose a matrix `B` into `PBQ = LU` where
///
/// * `P` is a row permutation
/// * `Q` is a column permutation
/// * `L` is lower triangular with `1`'s on the diagonal
/// * `U` is upper triangular
///
/// Note that permutations `P` and `Q` have the transpose equal to their inverse as they are
/// orthogonal matrices.
///
/// `P` and `Q` are "full" permutations, not to be confused with the simpler "rotating" permutations
/// in the `updates` field.
#[derive(Eq, PartialEq, Clone, Debug)]
pub struct LUDecomposition<F> {
    /// Row permutation `P`.
    ///
    /// The `forward` application of the permutation to rows of `M` corresponds to `PM`.
    /// FROM INITAL LU
    row_permutation: FullPermutation,
    /// Column permutation `Q`.
    ///
    /// The `backward`
    /// FROM INITAL LU
    column_permutation: FullPermutation,
    /// Lower triangular matrix `L`.
    ///
    /// Column major, one's on diagonal implied. So the length of the vector is `m - 1`.
    lower_triangular: Vec<Vec<(usize, F)>>,
    /// Upper triangular matrix `U`.
    ///
    /// Column major.
    // TODO(PERFORMANCE): Consider storing the diagonal separately.
    upper_triangular: Vec<Vec<(usize, F)>>,

    updates: Vec<(FullPermutation, FullPermutation)>, // Store row and column permutations of one update (row_update, column_update)
}

impl<F> BasisInverse for LUDecomposition<F>
where
    F: ops::Internal + ops::InternalHR,
{
    type F = F;
    type ColumnComputationInfo = ColumnAndSpike<Self::F>;

    fn identity(m: usize) -> Self {
        Self {
            row_permutation: FullPermutation::identity(m),
            column_permutation: FullPermutation::identity(m),
            lower_triangular: vec![vec![]; m - 1],
            upper_triangular: (0..m).map(|i| vec![(i, F::one())]).collect(),
            updates: vec![],
        }
    }

    fn invert<C: Column + OrderedColumn>(columns: Vec<C>) -> Self
    where
        Self::F: ops::Column<C::F>,
    {
        let m = columns.len();
        let mut rows = vec![Vec::new(); m];
        for (j, column) in columns.into_iter().enumerate() {
            for (i, value) in column.iter() {
                rows[*i].push((j, value.into()));
            }
        }
        debug_assert!(rows.iter().all(|row| row.is_sorted_by_key(|&(j, _)| j)));

        Self::rows(rows)
    }

    fn change_basis(&mut self, pivot_row_index: usize, column: Self::ColumnComputationInfo) {
        let m = self.m();
        // Index in entering column

        let pivot_column_index = {
            // Simplex iteration pivots on `pivot_row_index` p
            // -> Replace column p of Basis by column a_q (column of pivot in row p) of A
            // -> `pivot_column_index` represents this column
            //
            // Column with a pivot in `pivot_row_index` is leaving
            let mut pivot_column_index = pivot_row_index;
            // Compute and store the column permutations applied through Q
            // -> Q places spike in column p in position m and moves other columns to the left
            self.column_permutation.forward(&mut pivot_column_index);
            for (p, q) in &self.updates {
                Permutation::forward(q, &mut pivot_column_index);
                Permutation::forward(p, &mut pivot_column_index);
            }
            pivot_column_index
        };

        //********************************
        // Initial active block boundaries:
        //********************************
        // Assume index starting from 0 and being inclusive
        let Self::ColumnComputationInfo {
            column: _,
            mut spike,
        } = column;

        self.upper_triangular[pivot_column_index] = spike;

        let mut active_block_row = self.upper_triangular[pivot_column_index].last().unwrap().0;
        let mut active_block_column = pivot_column_index;
        let l = &mut self.lower_triangular;
        let u = &mut self.upper_triangular;
        if active_block_column == active_block_row {
            self.updates.push((
                FullPermutation::identity(u.len()),
                FullPermutation::identity(u.len()),
            ));
            return;
        }
        // Make L into square matrix
        l.push(vec![]);

        // lower triangular to unit lower triangular
        for col in 0..l.len() {
            l[col].insert(0, (col, F::one()));
        }

        // Get blocks needed for update
        let mut u_12: (usize, Vec<Vec<(usize, F)>>) =
            get_block(&u, 1, 2, active_block_column, active_block_row);

        let mut u_23: (usize, Vec<Vec<(usize, F)>>) =
            get_block(&u, 2, 3, active_block_column, active_block_row);

        let mut u_22: (usize, Vec<Vec<(usize, F)>>) =
            get_block(&u, 2, 2, active_block_column, active_block_row);

        let mut l_21: (usize, Vec<Vec<(usize, F)>>) =
            get_block(&l, 2, 1, active_block_column, active_block_row);

        let l_22: (usize, Vec<Vec<(usize, F)>>) =
            get_block(&l, 2, 2, active_block_column, active_block_row);

        let mut l_32: (usize, Vec<Vec<(usize, F)>>) =
            get_block(&l, 3, 2, active_block_column, active_block_row);

        let active_block_product = multiply_matrices(&l_22, &u_22);

        // Bring product into row major for LU decomposition
        let mut active_block_product_T =
            (active_block_product.0, vec![vec![]; active_block_product.0]);
        for j in 0..active_block_product.0 {
            for (i, x) in &active_block_product.1[j] {
                active_block_product_T.1[*i].push((j, x.clone()));
            }
        }

        // Compute LU of active block
        let mut active_block_lu = LUDecomposition::<F>::rows(active_block_product_T.1);
        let mut l_bar = (
            active_block_product.0,
            active_block_lu.lower_triangular.clone(),
        );
        let mut u_bar = (
            active_block_product.0,
            active_block_lu.upper_triangular.clone(),
        );
        // make l_bar square
        l_bar.1.push(vec![]);

        let mut id = vec![vec![]; active_block_product.0];
        // l_bar to unit lower triangular
        for col in 0..l_bar.1.len() {
            l_bar.1[col].insert(0, (col, F::one()));
            id[col].push((col, F::one()));
        }
        let mut p_bar = generate_permutation_matrix(&active_block_lu.row_permutation);
        let mut q_bar = generate_permutation_matrix(&active_block_lu.column_permutation);

        // Compute required inversions for updates
        &active_block_lu.column_permutation.invert();
        &active_block_lu.row_permutation.invert();

        let mut l_bar_inv = Vec::new();
        let mut l_22_inv = Vec::new();

        for i in 0..active_block_product.0 {
            let unit_column =
                BTreeMap::from_iter(IntoIter::new([(id[i][0].0, id[i][0].1.clone())]));

            l_bar_inv.push(active_block_lu.invert_lower_right(unit_column));
        }
        let l_bar_inv = (active_block_product.0, l_bar_inv);

        // remove unit diagonal and last column for inversion in LU decomposition
        for col in 0..l_22.1.len() - 1 {
            active_block_lu.lower_triangular[col] = l_22.1[col].clone();
            active_block_lu.lower_triangular[col].remove(0);
        }

        for i in 0..active_block_product.0 - 1 {
            let mut unit_column = BTreeMap::new();
            unit_column.insert(i, F::one());

            l_22_inv.push(active_block_lu.invert_lower_right(unit_column));
        }
        // make l_22 square
        l_22_inv.push(vec![]);
        // Make l_22 inv unit triangular
        for col in 0..l_22.1.len() {
            l_22_inv[col].insert(0, (col, F::one()));
        }

        let l_22_inv = (active_block_product.0, l_22_inv);

        let p_bar_inv = generate_permutation_matrix::<F>(&active_block_lu.row_permutation);
        let q_bar_inv = generate_permutation_matrix::<F>(&active_block_lu.column_permutation);

        // Compute new entries for L matrix
        l_21 = multiply_matrices(&p_bar_inv, &l_21);
        l_32 = multiply_matrices(
            &multiply_matrices(&l_32, &l_22_inv),
            &multiply_matrices(&p_bar, &l_bar),
        );

        // Compute new entries for U matrix
        // TODO(Debug): u_23 is not sorted? might be a bug. insert block deals with this by using binary search and insert
        u_12 = multiply_matrices(&u_12, &q_bar);
        u_23 = multiply_matrices(
            &multiply_matrices(&l_bar_inv, &p_bar_inv),
            &multiply_matrices(&l_22, &u_23),
        );

        // Restore L and U to upper triangular
        insert_block(l, &l_bar, active_block_column, active_block_column);
        insert_block(l, &l_32, active_block_row + 1, active_block_column);
        insert_block(l, &l_21, active_block_column, 0);
        insert_block(u, &u_12, 0, active_block_column);
        insert_block(u, &u_bar, active_block_column, active_block_column);
        insert_block(u, &u_23, active_block_column, active_block_row + 1);

        let row_update = FullPermutation::new(
            (0..l.len())
                .map(|j| {
                    if j >= active_block_column && j <= active_block_row {
                        p_bar.1[j - active_block_column][0].0 + active_block_column
                    } else {
                        j
                    }
                })
                .collect(),
        );
        let column_update = FullPermutation::new(
            (0..l.len())
                .map(|j| {
                    if j >= active_block_column && j <= active_block_row {
                        q_bar.1[j - active_block_column][0].0 + active_block_column
                    } else {
                        j
                    }
                })
                .collect(),
        );
        // Restore l into non square and non unit lower
        // TODO(Debug): might be problematic if L is not unit lower_triangular but only lower triangular
        l.pop();
        for i in (0..l.len()) {
            if !l[i].is_empty() {
                l[i].remove(0);
            }
        }
        self.updates.push((row_update, column_update));
    }

    fn generate_column<C: Column>(&self, original_column: C) -> Self::ColumnComputationInfo
    where
        Self::F: ops::Column<C::F>,
    {
        // PBQ = LU

        // goal: y = B^{-1} c

        // apply P and Q permutations of initial LU to rhs
        let rhs = original_column
            .iter()
            .map(|(mut i, v)| {
                self.row_permutation.forward(&mut i);
                (i, v.into())
            })
            // Also sorts after the row permutation
            .collect::<BTreeMap<_, _>>();

        // TODO(PERFORMANCE): refactor such that this is not required
        let mut rhs: Vec<(usize, F)> = Vec::from_iter(rhs.into_iter());

        // apply row updates to rhs
        for (p, _) in self.updates.iter() {
            // apply last inverse row permutation of updates vector
            rhs = rhs
                .iter()
                .map(|(mut i, v)| {
                    p.backward(&mut i);
                    (i, v.clone())
                })
                .collect();
            // apply last inverse column permutation of updates vector
            // q.backward_unsorted(&mut rhs[..])
        }

        let rhs = BTreeMap::from_iter(rhs.into_iter());

        // Compute L^-1 c
        let mut w = self.invert_lower_right(rhs);

        let spike = w.clone();

        // Compute U^-1 c
        let mut column = self.invert_upper_right(w.into_iter().collect());

        // Apply column permutation updates
        for (_, q) in self.updates.iter().rev() {
            q.backward_sorted(&mut column);
        }
        // Apply initial column permutation
        self.column_permutation.forward_sorted(&mut column);
        column.sort_unstable_by_key(|&(i, _)| i);

        Self::ColumnComputationInfo {
            column: SparseVector::new(column, self.m()),
            spike,
        }
    }

    fn generate_element<C: Column + OrderedColumn>(
        &self,
        i: usize,
        original_column: C,
    ) -> Option<Self::F>
    where
        Self::F: ops::Column<C::F>,
    {
        self.generate_column(original_column)
            .into_column()
            .get(i)
            .cloned()
    }

    fn should_refactor(&self) -> bool {
        // TODO(ENHANCEMENT): What would be a good decision rule?
        self.updates.len() > 10
    }

    fn basis_inverse_row(&self, mut row: usize) -> SparseVector<Self::F, Self::F> {
        self.column_permutation.forward(&mut row);

        for (_, q) in &self.updates {
            q.backward(&mut row);
        }

        // unit vector where v[row]=1
        let initial_rhs = iter::once((row, Self::F::one())).collect();

        let mut w = self.invert_upper_left(initial_rhs);
        let mut tuples = self.invert_lower_left(w.into_iter().collect());

        for (p, _) in self.updates.iter().rev() {
            p.backward_sorted(&mut tuples);
        }

        self.row_permutation.backward_sorted(&mut tuples);

        SparseVector::new(tuples, self.m())
    }

    fn m(&self) -> usize {
        self.row_permutation.len()
        // == self.column_permutation.len()
        // == self.upper_triangular.len()
        // == self.lower_triangular.len() + 1
    }
}

impl<F> LUDecomposition<F>
where
    F: ops::Internal + ops::InternalHR,
{
    fn invert_lower_right(&self, mut rhs: BTreeMap<usize, F>) -> Vec<(usize, F)> {
        let mut result = Vec::new();

        while let Some((row, rhs_value)) = rhs.pop_first() {
            let column = row;

            let result_value = rhs_value; // Implicit 1 on the diagonal.

            if column != self.m() - 1 {
                // Introduce new nonzeros
                for &(row, ref value) in &self.lower_triangular[column] {
                    insert_or_shift_maybe_remove(row, &result_value * value, &mut rhs);
                }
            }

            result.push((row, result_value));
        }

        result
    }

    fn invert_upper_right(&self, mut rhs: BTreeMap<usize, F>) -> Vec<(usize, F)> {
        let mut result = Vec::new();

        while let Some((row, rhs_value)) = rhs.pop_last() {
            let column = row;
            let result_value = self.compute_result(column, rhs_value);
            self.update_rhs(column, &result_value, &mut rhs);
            result.push((row, result_value));
        }

        result.reverse();
        debug_assert!(result.is_sorted_by_key(|&(i, _)| i));

        result
    }

    fn invert_upper_right_row(&self, mut rhs: BTreeMap<usize, F>, target_row: usize) -> Option<F> {
        while let Some((row, rhs_value)) = rhs.pop_last() {
            let column = row;
            match row.cmp(&target_row) {
                Ordering::Less => return None,
                Ordering::Equal => return Some(self.compute_result(column, rhs_value)),
                Ordering::Greater => {
                    let result_value = self.compute_result(column, rhs_value);
                    self.update_rhs(column, &result_value, &mut rhs);
                }
            }
        }

        None
    }

    fn compute_result(&self, column: usize, rhs_value: F) -> F {
        debug_assert_eq!(
            self.upper_triangular[column].last().unwrap().0,
            column,
            "Needs to have a diagonal element",
        );

        let diagonal_item = &self.upper_triangular[column].last().unwrap().1;
        rhs_value / diagonal_item
    }

    fn update_rhs(&self, column: usize, result_value: &F, rhs: &mut BTreeMap<usize, F>) {
        let nr_column_items = self.upper_triangular[column].len();
        for &(row, ref value) in &self.upper_triangular[column][..(nr_column_items - 1)] {
            insert_or_shift_maybe_remove(row, result_value * value, rhs);
        }
    }

    fn invert_lower_left(&self, mut rhs: BTreeMap<usize, F>) -> Vec<(usize, F)> {
        let mut result = Vec::new();

        while let Some((column, rhs_value)) = rhs.pop_last() {
            let row = column;
            let result_value = rhs_value;

            // Introduce new nonzeros, by rows
            for j in 0..column {
                // TODO(PERFORMANCE): Avoid scanning all columns for row values
                let has_row = self.lower_triangular[j].binary_search_by_key(&row, |&(i, _)| i);
                if let Ok(data_index) = has_row {
                    let value = &self.lower_triangular[j][data_index].1;
                    insert_or_shift_maybe_remove(j, &result_value * value, &mut rhs);
                }
            }

            result.push((row, result_value));
        }

        result.reverse();
        debug_assert!(result.is_sorted_by_key(|&(i, _)| i));

        result
    }

    fn invert_upper_left(&self, mut rhs: BTreeMap<usize, F>) -> Vec<(usize, F)> {
        let mut result = Vec::new();

        while let Some((column, rhs_value)) = rhs.pop_first() {
            let row = column;
            let diagonal_item = &self.upper_triangular[column].last().unwrap().1;
            let result_value = rhs_value / diagonal_item;

            // Introduce new nonzeros, by rows
            for j in (column + 1)..self.m() {
                // TODO(PERFORMANCE): Avoid scanning all columns for row values
                let has_row = self.upper_triangular[j].binary_search_by_key(&row, |&(i, _)| i);
                if let Ok(data_index) = has_row {
                    let value = &self.upper_triangular[j][data_index].1;
                    insert_or_shift_maybe_remove(j, &result_value * value, &mut rhs);
                }
            }

            result.push((column, result_value));
        }

        debug_assert!(result.windows(2).all(|w| w[0].0 < w[1].0));

        result
    }
}

fn generate_permutation_matrix<F: num::One>(p: &FullPermutation) -> (usize, Vec<Vec<(usize, F)>>) {
    (
        p.len(),
        (0..p.len())
            .map(|mut j| {
                p.forward(&mut j);
                vec![(j, F::one())]
            })
            .collect(),
    )
}

fn multiply_matrices<F: Clone + ops::Internal + ops::InternalHR>(
    a: &(usize, Vec<Vec<(usize, F)>>),
    b: &(usize, Vec<Vec<(usize, F)>>),
) -> (usize, Vec<Vec<(usize, F)>>) {
    assert_eq!(a.1.len(), b.0);
    let a_size = (a.0, a.1.len());
    let b_size = (b.0, b.1.len());
    // TODO(Debug): might need fix in b_size.1
    let mut product = vec![vec![]; b_size.1];
    for column in 0..b_size.1 {
        let mut column_product: Vec<Vec<(usize, F)>> = Vec::new();
        let mut row_indices: Vec<usize> = vec![];
        for (row, val) in &b.1[column] {
            column_product.push(a.1[*row].iter().map(|(i, x)| (*i, x * val)).collect());
            let mut temp: Vec<usize> = column_product
                .last()
                .unwrap()
                .iter()
                .map(|&(i, _)| i)
                .filter(|i| !row_indices.contains(i))
                .collect();
            row_indices.append(&mut temp);
        }
        for row in row_indices {
            let mut row_sum = F::zero();
            let column_indices = (0..column_product.len()).filter_map(|j| {
                column_product[j]
                    .binary_search_by_key(&row, |&(i, _)| i)
                    .ok()
                    .map(|i| (i, j))
            });
            for (i, j) in column_indices {
                row_sum += &column_product[j][i].1;
            }
            if row_sum != F::zero() {
                product[column].push((row, row_sum));
            }
        }
    }
    (a_size.0, product)
}

fn get_block<F: Clone>(
    m: &Vec<Vec<(usize, F)>>,
    row: usize,    // Row of block in block representation of matrix M (exclusive)
    column: usize, // Column of block in block representation of matrix M (exclusive)
    active_block_column: usize, // Inclusive
    active_block_row: usize, // Inclusive
) -> (usize, Vec<Vec<(usize, F)>>) {
    // columns in block
    let (c1, c2) = match column {
        1 => (0, active_block_column),
        2 => (active_block_column, active_block_row + 1),
        _ => (active_block_row + 1, m.len()),
    };

    // rows in block
    let (r1, r2) = match row {
        1 => (0, active_block_column),
        2 => (active_block_column, active_block_row + 1),
        _ => (active_block_row + 1, m.len()),
    };

    let mut block: (usize, Vec<Vec<(usize, F)>>) = (
        r2 - r1,
        m[c1..c2]
            .iter()
            .map(|vec| {
                vec.iter()
                    .filter(|&(i, _)| r1 <= *i && *i < r2)
                    .map(|(i, x)| (*i - r1, x.clone()))
                    .collect()
            })
            .collect(),
    );
    block
}

// Inserts block into matrix inplace
// Assume block has exclusive column length in first arg of tuple
fn insert_block<F: Clone>(
    m: &mut Vec<Vec<(usize, F)>>,          // Matrix M to insert into
    block: &(usize, Vec<Vec<(usize, F)>>), // Block to insert
    row: usize,                            // Row where the block begins in M (inclusive)
    column: usize,                         // Column          "
) {
    // Restore correct row indices in relation to M
    let aligned_block: (usize, Vec<Vec<(usize, F)>>) = (
        block.0,
        (0..block.1.len())
            .map(|j| {
                block.1[j]
                    .iter()
                    .map(|(i, x)| (*i + row, x.clone()))
                    .collect()
            })
            .collect(),
    );
    // Update the row entries belonging to the block in M
    for j in (column..column + aligned_block.1.len()) {
        // zero out existing entries
        m[j] = m[j]
            .iter()
            .filter(|&(i, _)| *i < row || *i >= (row + aligned_block.0))
            .cloned()
            .collect();
        for i in 0..aligned_block.1[j - column].len() {
            let new_entry = &aligned_block.1[j - column][i].0;
            let search = m[j].binary_search_by_key(new_entry, |&(i, _)| i);
            match search {
                Ok(x) => m[j][x] = aligned_block.1[j - column][i].clone(),
                Err(x) => m[j].insert(x, aligned_block.1[j - column][i].clone()),
            }
        }
    }
}

fn insert_or_shift_maybe_remove<F>(index: usize, change: F, rhs: &mut BTreeMap<usize, F>)
where
    F: ops::Internal + ops::InternalHR + Zero,
{
    match rhs.get_mut(&index) {
        None => {
            rhs.insert(index, -change);
        }
        Some(existing) => {
            *existing -= change;
            if existing.is_zero() {
                rhs.remove(&index);
            }
        }
    }
}

#[derive(Debug)]
pub struct ColumnAndSpike<F> {
    column: SparseVector<F, F>,
    spike: Vec<(usize, F)>,
}

impl<F: SparseElement<F>> ColumnComputationInfo<F> for ColumnAndSpike<F> {
    fn column(&self) -> &SparseVector<F, F> {
        &self.column
    }

    fn into_column(self) -> SparseVector<F, F> {
        self.column
    }
}

impl<F> Display for LUDecomposition<F>
where
    F: ops::Internal + ops::InternalHR,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let width = 10;
        let column_width = 3;

        writeln!(f, "Lower:")?;
        write!(f, "{:>width$} |", "", width = column_width)?;
        for j in 0..self.m() {
            write!(f, "{0:^width$}", j, width = width)?;
        }
        writeln!(f)?;
        let total_width = column_width + 1 + 1 + self.m() * width;
        writeln!(f, "{}", "-".repeat(total_width))?;

        for i in 0..self.m() {
            write!(f, "{0:>width$} |", i, width = column_width)?;
            for j in 0..self.m() {
                let value = match j.cmp(&i) {
                    Ordering::Greater => "".to_string(),
                    Ordering::Equal => "1".to_string(),
                    Ordering::Less => {
                        match self.lower_triangular[j].binary_search_by_key(&i, |&(i, _)| i) {
                            Ok(index) => self.lower_triangular[j][index].1.to_string(),
                            Err(_) => "0".to_string(),
                        }
                    }
                };
                write!(f, "{0:^width$}", value, width = width)?;
            }
            writeln!(f)?;
        }
        writeln!(f)?;

        writeln!(f, "Upper:")?;
        write!(f, "{:>width$} |", "", width = column_width)?;
        for j in 0..self.m() {
            write!(f, "{0:^width$}", j, width = width)?;
        }
        writeln!(f)?;
        writeln!(f, "{}", "-".repeat(total_width))?;

        for i in 0..self.m() {
            write!(f, "{0:>width$} |", i, width = column_width)?;
            for j in 0..self.m() {
                let value = match j.cmp(&i) {
                    Ordering::Equal | Ordering::Greater => {
                        match self.upper_triangular[j].binary_search_by_key(&i, |&(i, _)| i) {
                            Ok(index) => self.upper_triangular[j][index].1.to_string(),
                            Err(_) => "0".to_string(),
                        }
                    }
                    Ordering::Less => "".to_string(),
                };
                write!(f, "{0:^width$}", value, width = width)?;
            }
            writeln!(f)?;
        }
        writeln!(f)?;

        writeln!(f, "Row permutation:")?;
        writeln!(f, "{:?}", self.row_permutation)?;
        writeln!(f, "Column permutation:")?;
        writeln!(f, "{:?}", self.column_permutation)?;

        writeln!(f, "Updates:")?;
        for (i, (eta, t)) in self.updates.iter().enumerate() {
            writeln!(f, "Update {}: ", i)?;
            writeln!(f, "R: {:?}", eta)?;
            writeln!(f, "pivot index: {}", t)?;
        }
        writeln!(f)
    }
}
#[cfg(test)]
mod test {
    use std::collections::BTreeMap;

    use num::FromPrimitive;

    use crate::{R64, RB};
    use crate::algorithm::two_phase::matrix_provider::column::identity::{IdentityColumnStruct, One};
    use crate::algorithm::two_phase::matrix_provider::matrix_data;
    use crate::algorithm::two_phase::matrix_provider::matrix_data::Column;
    use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::BasisInverse;
    use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::lower_upper_remultiply_factor::ColumnAndSpike;
    use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::lower_upper_remultiply_factor::permutation::FullPermutation;
    use crate::algorithm::two_phase::tableau::inverse_maintenance::ColumnComputationInfo;
    use crate::data::linear_algebra::vector::{SparseVector, Vector};
    use crate::data::number_types::rational::{Rational64, RationalBig};
    use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::lower_upper_remultiply_factor::LUDecomposition;

    mod matmul {
        use super::*;

        #[test]
        fn identity_empty() {
            let identity = LUDecomposition::<RationalBig>::identity(2);

            let column = BTreeMap::new();
            let result = identity.invert_upper_right(column);
            assert!(result.is_empty());
            let column = BTreeMap::new();
            let result = identity.invert_upper_left(column);
            assert!(result.is_empty());
            let column = BTreeMap::new();
            let result = identity.invert_lower_right(column);
            assert!(result.is_empty());
            let column = BTreeMap::new();
            let result = identity.invert_lower_left(column);
            assert!(result.is_empty());
        }

        #[test]
        fn identity_single() {
            let identity = LUDecomposition::<RationalBig>::identity(2);

            let column = vec![(0, RB!(1))];
            let result = identity.invert_upper_right(column.clone().into_iter().collect());
            assert_eq!(result, column);
            let column = vec![(0, RB!(1))];
            let result = identity.invert_upper_left(column.clone().into_iter().collect());
            assert_eq!(result, column);
            let column = vec![(0, RB!(1))];
            let result = identity.invert_lower_right(column.clone().into_iter().collect());
            assert_eq!(result, column);
            let column = vec![(0, RB!(1))];
            let result = identity.invert_lower_left(column.clone().into_iter().collect());
            assert_eq!(result, column);

            let column = vec![(1, RB!(1))];
            let result = identity.invert_upper_right(column.clone().into_iter().collect());
            assert_eq!(result, column);
            let column = vec![(1, RB!(1))];
            let result = identity.invert_upper_left(column.clone().into_iter().collect());
            assert_eq!(result, column);
            let column = vec![(1, RB!(1))];
            let result = identity.invert_lower_right(column.clone().into_iter().collect());
            assert_eq!(result, column);
            let column = vec![(1, RB!(1))];
            let result = identity.invert_lower_left(column.clone().into_iter().collect());
            assert_eq!(result, column);
        }

        #[test]
        fn identity_double() {
            let identity = LUDecomposition::<RationalBig>::identity(2);

            let column = vec![(0, RB!(1)), (1, RB!(1))];
            let result = identity.invert_upper_right(column.clone().into_iter().collect());
            assert_eq!(result, column);
            let column = vec![(0, RB!(1)), (1, RB!(1))];
            let result = identity.invert_upper_left(column.clone().into_iter().collect());
            assert_eq!(result, column);
            let column = vec![(0, RB!(1)), (1, RB!(1))];
            let result = identity.invert_lower_right(column.clone().into_iter().collect());
            assert_eq!(result, column);
            let column = vec![(0, RB!(1)), (1, RB!(1))];
            let result = identity.invert_lower_left(column.clone().into_iter().collect());
            assert_eq!(result, column);
        }

        #[test]
        fn offdiagonal_empty() {
            let offdiag = LUDecomposition {
                row_permutation: FullPermutation::identity(2),
                column_permutation: FullPermutation::identity(2),
                lower_triangular: vec![vec![(1, RB!(1))]],
                upper_triangular: vec![vec![(0, RB!(1))], vec![(1, RB!(1))]],
                updates: vec![],
            };

            let column: Column<Rational64> = Column::Sparse {
                constraint_values: vec![],
                slack: None,
                mock_array: [],
            };
            let result = offdiag.generate_column(column);
            assert_eq!(result.column, SparseVector::new(vec![], 2));
        }

        #[test]
        fn offdiagonal_single() {
            let offdiag = LUDecomposition {
                row_permutation: FullPermutation::identity(2),
                column_permutation: FullPermutation::identity(2),
                lower_triangular: vec![vec![(1, RB!(1))]],
                upper_triangular: vec![vec![(0, RB!(1))], vec![(1, RB!(1))]],
                updates: vec![],
            };

            let column = Column::Sparse {
                constraint_values: vec![(0, R64!(1))],
                slack: None,
                mock_array: [],
            };
            let result = offdiag.generate_column(column);
            assert_eq!(
                result.column,
                SparseVector::new(vec![(0, RB!(1)), (1, RB!(-1))], 2)
            );

            let column = Column::Sparse {
                constraint_values: vec![(1, R64!(1))],
                slack: None,
                mock_array: [],
            };
            let result = offdiag.generate_column(column);
            assert_eq!(result.column, SparseVector::new(vec![(1, RB!(1))], 2));
        }
    }

    mod change_basis {

        use super::*;

        /// Spike is the column which is already there.
        #[test]
        fn no_change() {
            let mut initial = LUDecomposition::<RationalBig>::identity(3);

            let spike = vec![(1, RB!(1))];
            let column_computation_info = ColumnAndSpike {
                column: SparseVector::new(spike.clone(), 3),
                spike,
            };
            initial.change_basis(1, column_computation_info);
            let modified = initial;

            let mut expected = LUDecomposition::<RationalBig>::identity(3);
            expected
                .updates
                .push((FullPermutation::identity(3), FullPermutation::identity(3)));
            assert_eq!(modified, expected);
        }
    }

    #[test]
    fn from_identity_2() {
        let mut identity = LUDecomposition::<RationalBig>::identity(2);

        let spike = vec![(0, RB!(1)), (1, RB!(1))];
        let column_computation_info = ColumnAndSpike {
            column: SparseVector::new(spike.clone(), 2),
            spike,
        };
        identity.change_basis(0, column_computation_info);
        let modified = identity;

        let expected = LUDecomposition {
            row_permutation: FullPermutation::identity(2),
            column_permutation: FullPermutation::identity(2),
            lower_triangular: vec![vec![(1, RB!(1))]],
            upper_triangular: vec![vec![(0, RB!(1))], vec![(1, RB!(1))]],
            updates: vec![(
                FullPermutation::new(vec![0, 1]),
                FullPermutation::new(vec![0, 1]),
            )],
        };
        assert_eq!(modified, expected);
    }

    /// Doesn't require any `r`, permutations only are sufficient
    /// TODO: small example
    #[test]
    fn from_5x5_identity_no_r() {
        let m = 5;
        let mut initial = LUDecomposition::<RationalBig>::identity(5);

        let spike = vec![(0, RB!(2)), (1, RB!(3)), (2, RB!(5)), (3, RB!(7))];
        let column_computation_info = ColumnAndSpike {
            column: SparseVector::new(spike.clone(), m),
            spike,
        };
        initial.change_basis(1, column_computation_info);
        let modified = initial;

        let expected = LUDecomposition {
            row_permutation: FullPermutation::identity(m),
            column_permutation: FullPermutation::identity(m),
            lower_triangular: vec![vec![], vec![(2, RB!(5, 3)), (3, RB!(7, 3))], vec![], vec![]],
            upper_triangular: vec![
                vec![(0, RB!(1))],
                vec![(0, RB!(2)), (1, RB!(3))],
                vec![(2, RB!(1))],
                vec![(3, RB!(1))],
                vec![(4, RB!(1))],
            ],
            updates: vec![(FullPermutation::identity(m), FullPermutation::identity(m))],
        };
        assert_eq!(modified, expected);
    }

    /// Does require an `r`, permutations are not sufficient.
    #[test]
    fn from_4x4_identity() {
        let m = 4;
        let mut initial = LUDecomposition {
            row_permutation: FullPermutation::identity(m),
            column_permutation: FullPermutation::identity(m),
            lower_triangular: vec![vec![]; m - 1],
            upper_triangular: vec![
                vec![(0, RB!(1))],
                vec![(1, RB!(1))],
                vec![(2, RB!(4))],
                vec![(1, RB!(5)), (3, RB!(6))],
            ],
            updates: vec![],
        };

        let spike = vec![(1, RB!(2)), (2, RB!(3)), (3, RB!(4))];
        let column_computation_info = ColumnAndSpike {
            // They are the same in this case, row permutation and lower are identity
            column: SparseVector::new(spike.clone(), m),
            spike,
        };
        initial.change_basis(1, column_computation_info);
        let modified = initial;

        let expected = LUDecomposition {
            row_permutation: FullPermutation::identity(m),
            column_permutation: FullPermutation::identity(m),
            lower_triangular: vec![vec![], vec![], vec![(3, RB!(2))]],
            upper_triangular: vec![
                vec![(0, RB!(1, 1))],
                vec![(1, RB!(4, 1))],
                vec![(1, RB!(3, 1)), (2, RB!(2, 1))],
                vec![(2, RB!(5, 1)), (3, RB!(-4, 1))],
            ],
            updates: vec![(
                FullPermutation::new(vec![0, 2, 1, 3]),
                FullPermutation::new(vec![0, 2, 1, 3]),
            )],
        };
        assert_eq!(modified, expected);

        // Columns
        assert_eq!(
            modified
                .generate_column(IdentityColumnStruct((0, One)))
                .into_column(),
            SparseVector::standard_basis_vector(0, m),
        );
        assert_eq!(
            modified
                .generate_column(IdentityColumnStruct((1, One)))
                .into_column(),
            SparseVector::new(vec![(1, RB!(-3, 4)), (2, RB!(9, 16)), (3, RB!(1, 2))], m),
        );
        assert_eq!(
            modified
                .generate_column(IdentityColumnStruct((2, One)))
                .into_column(),
            SparseVector::new(vec![(2, RB!(1, 4))], m),
        );
        assert_eq!(
            modified
                .generate_column(IdentityColumnStruct((3, One)))
                .into_column(),
            SparseVector::new(vec![(1, RB!(5, 8)), (2, RB!(-15, 32)), (3, RB!(-1, 4))], m),
        );

        // Rows
        assert_eq!(
            modified.basis_inverse_row(0),
            SparseVector::standard_basis_vector(0, m),
        );
        assert_eq!(
            modified.basis_inverse_row(1),
            SparseVector::new(vec![(1, RB!(-3, 4)), (3, RB!(5, 8))], m),
        );
        assert_eq!(
            modified.basis_inverse_row(2),
            SparseVector::new(vec![(1, RB!(9, 16)), (2, RB!(1, 4)), (3, RB!(-15, 32))], m),
        );
        assert_eq!(
            modified.basis_inverse_row(3),
            SparseVector::new(vec![(1, RB!(1, 2)), (3, RB!(-1, 4))], m),
        );
    }

    ///From "A review of the LU update in the simplex algorithm" by Joseph M. Elble and
    ///Nikolaes V. Sahinidis, Int. J. Mathematics in Operational Research, Vol. 4, No. 4, 2012.
    #[test]
    fn from_5x5_identity() {
        let m = 5;
        let mut initial = LUDecomposition {
            row_permutation: FullPermutation::identity(m),
            column_permutation: FullPermutation::identity(m),
            lower_triangular: vec![vec![]; m - 1],
            upper_triangular: vec![
                vec![(0, RB!(11))],
                vec![(0, RB!(12)), (1, RB!(22))],
                vec![(0, RB!(13)), (1, RB!(23)), (2, RB!(33))],
                vec![(0, RB!(14)), (1, RB!(24)), (2, RB!(34)), (3, RB!(44))],
                vec![
                    (0, RB!(15)),
                    (1, RB!(25)),
                    (2, RB!(35)),
                    (3, RB!(45)),
                    (4, RB!(55)),
                ],
            ],
            updates: vec![],
        };

        let spike = vec![(0, RB!(12)), (1, RB!(22)), (2, RB!(32)), (3, RB!(42))];
        let column_computation_info = ColumnAndSpike {
            column: SparseVector::new(spike.clone(), m),
            spike,
        };
        initial.change_basis(1, column_computation_info);
        let modified = initial;

        let expected = LUDecomposition {
            row_permutation: FullPermutation::identity(m),
            column_permutation: FullPermutation::identity(m),
            lower_triangular: vec![
                vec![],
                vec![(2, RB!(16, 21)), (3, RB!(11, 21))],
                vec![(3, RB!(23, 33))],
                vec![],
            ],
            upper_triangular: vec![
                vec![(0, RB!(11, 1))],
                vec![(0, RB!(12, 1)), (1, RB!(42, 1))],
                vec![(0, RB!(13, 1)), (2, RB!(33, 1))],
                vec![
                    (0, RB!(14, 1)),
                    (1, RB!(44, 1)),
                    (2, RB!(10, 21)),
                    (3, RB!(430, 693)),
                ],
                vec![
                    (0, RB!(15, 1)),
                    (1, RB!(45, 1)),
                    (2, RB!(5, 7)),
                    (3, RB!(215, 231)),
                    (4, RB!(55, 1)),
                ],
            ],
            updates: vec![(
                FullPermutation::new(vec![0, 3, 2, 1, 4]),
                FullPermutation::identity(m),
            )],
        };
        assert_eq!(modified, expected);

        // Columns
        assert_eq!(
            modified
                .generate_column(IdentityColumnStruct((0, One)))
                .into_column(),
            SparseVector::new(vec![(0, RB!(1, 11))], m),
        );
        assert_eq!(
            modified
                .generate_column(IdentityColumnStruct((1, One)))
                .into_column(),
            SparseVector::new(
                vec![
                    (0, RB!(-2, 11)),
                    (1, RB!(-363, 215)),
                    (2, RB!(-1, 43)),
                    (3, RB!(693, 430))
                ],
                m
            ),
        );
        assert_eq!(
            modified
                .generate_column(IdentityColumnStruct((2, One)))
                .into_column(),
            SparseVector::new(
                vec![
                    (0, RB!(1, 11)),
                    (1, RB!(253, 215)),
                    (2, RB!(2, 43)),
                    (3, RB!(-483, 430))
                ],
                m
            ),
        );
        assert_eq!(
            modified
                .generate_column(IdentityColumnStruct((3, One)))
                .into_column(),
            SparseVector::new(vec![(1, RB!(1, 86)), (2, RB!(-1, 43)), (3, RB!(1, 86))], m),
        );
        assert_eq!(
            modified
                .generate_column(IdentityColumnStruct((4, One)))
                .into_column(),
            SparseVector::new(
                vec![(1, RB!(1, 110)), (3, RB!(-3, 110)), (4, RB!(1, 55))],
                m
            ),
        );
        // Sum of two
        assert_eq!(
            modified
                .generate_column(matrix_data::Column::TwoSlack(
                    [(0, RB!(1)), (1, RB!(1))],
                    []
                ))
                .into_column(),
            SparseVector::new(
                vec![
                    (0, RB!(-1, 11)),
                    (1, RB!(-363, 215)),
                    (2, RB!(-1, 43)),
                    (3, RB!(693, 430))
                ],
                m
            ),
        );

        // Rows
        assert_eq!(
            modified.basis_inverse_row(0),
            SparseVector::new(vec![(0, RB!(1, 11)), (1, RB!(-2, 11)), (2, RB!(1, 11))], m),
        );
        assert_eq!(
            modified.basis_inverse_row(1),
            SparseVector::new(
                vec![
                    (1, RB!(-363, 215)),
                    (2, RB!(253, 215)),
                    (3, RB!(1, 86)),
                    (4, RB!(1, 110))
                ],
                m
            ),
        );
        assert_eq!(
            modified.basis_inverse_row(2),
            SparseVector::new(vec![(1, RB!(-1, 43)), (2, RB!(2, 43)), (3, RB!(-1, 43))], m),
        );
        assert_eq!(
            modified.basis_inverse_row(3),
            SparseVector::new(
                vec![
                    (1, RB!(693, 430)),
                    (2, RB!(-483, 430)),
                    (3, RB!(1, 86)),
                    (4, RB!(-3, 110))
                ],
                m
            ),
        );
        assert_eq!(
            modified.basis_inverse_row(4),
            SparseVector::new(vec![(4, RB!(1, 55))], m),
        );
    }
}
