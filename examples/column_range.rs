use num::FromPrimitive;

use rust_lp::algorithm::two_phase::matrix_provider::matrix_data::MatrixData;
use rust_lp::algorithm::two_phase::tableau::inverse_maintenance::carry::basis_inverse_rows::BasisInverseRows;
use rust_lp::algorithm::two_phase::tableau::inverse_maintenance::carry::Carry;
use rust_lp::data::linear_algebra::matrix::{ColumnMajor, Order};
use rust_lp::data::linear_algebra::vector::{DenseVector, Vector};
use rust_lp::data::linear_program::elements::VariableType;
use rust_lp::data::linear_program::general_form::Variable;
use rust_lp::data::number_types::rational::RationalBig;
use rust_lp::RB;
use rust_lp::algorithm::two_phase::tableau::inverse_maintenance::InverseMaintener;
use rust_lp::algorithm::two_phase::tableau::kind::non_artificial::NonArtificial;
use rust_lp::algorithm::two_phase::phase_two;
use rust_lp::algorithm::two_phase::strategy::pivot_rule::FirstProfitable;
use rust_lp::algorithm::two_phase::tableau::Tableau;

fn main() {
    let input_matrix = [
        [3, 3, 3],
        [2, 3, 3],
        [1, 2, 3],
    ];

    let m = input_matrix.len();
    let n = input_matrix[0].len();

    let column_extreme = |f: fn(_) -> Option<i32>| (0..n)
        .map(|j| f((0..m).map(move |i| input_matrix[i][j])).unwrap())
        .collect::<Vec<_>>();
    let column_min = column_extreme(Iterator::min);
    let column_max = column_extreme(Iterator::max);

    // Variables
    let x = (0..m).map(|_| Variable {
            variable_type: VariableType::Continuous,
            cost: RB!(0),
            lower_bound: Some(RB!(0)),
            upper_bound: Some(RB!(1)),
            shift: RB!(0),
            flipped: false,
        })
        .collect::<Vec<_>>();
    let m_lower = (0..n).map(|_| Variable {
        variable_type: VariableType::Continuous,
        cost: RB!(-1),
        lower_bound: Some(RB!(0)),
        upper_bound: None,
        shift: RB!(0),
        flipped: false,
    })
        .collect::<Vec<_>>();
    let m_upper = (0..n).map(|_| Variable {
        variable_type: VariableType::Continuous,
        cost: RB!(1),
        lower_bound: Some(RB!(0)),
        upper_bound: None,
        shift: RB!(0),
        flipped: false,
    })
        .collect::<Vec<_>>();
    let variables = [x, m_lower, m_upper].concat();

    // Constraints
    let mut row_major_constraints = Vec::new();
    let mut b = Vec::new();
    let mut basis_columns = Vec::new();
    let mut nr_upper_bounded_constraints = 0;
    let mut had_a_min = vec![false; n];
    for j in 0..n {
        for i in 0..m {
            if input_matrix[i][j] == column_min[j] {
                row_major_constraints.push(vec![(i, RB!(1)), (m + j, RB!(1))]);
                b.push(RB!(input_matrix[i][j]));
                basis_columns.push(if had_a_min[j] {
                    m + 2 * n + nr_upper_bounded_constraints
                } else {
                    had_a_min[j] = true;
                    m + j
                });
                nr_upper_bounded_constraints += 1;
            }
        }
    }
    let mut nr_lower_bounded_constraints = 0;
    let mut had_a_max = vec![false; n];
    for j in 0..n {
        for i in 0..m {
            if input_matrix[i][j] == column_max[j] {
                row_major_constraints.push(vec![(i, RB!(1)), (m + n + j, RB!(1))]);
                b.push(RB!(input_matrix[i][j]));
                basis_columns.push(if had_a_max[j] {
                    m + 2 * n + nr_upper_bounded_constraints + nr_lower_bounded_constraints
                } else {
                    had_a_max[j] = true;
                    m + n + j
                });
                nr_lower_bounded_constraints += 1;
            }
        }
    }

    for i in 0..m {
        basis_columns.push(m + 2 * n + nr_upper_bounded_constraints + nr_lower_bounded_constraints + i);
    }

    let mut constraints = vec![vec![]; m + 2 * n];
    for (row_index, row) in row_major_constraints.into_iter().enumerate() {
        for (column_index, value) in row {
            constraints[column_index].push((row_index, value));
        }
    }
    let constraints = ColumnMajor::new(constraints, nr_upper_bounded_constraints + nr_lower_bounded_constraints, m + 2 * n);
    let b = DenseVector::new(b, nr_upper_bounded_constraints + nr_lower_bounded_constraints);

    let matrix = MatrixData::new(
        &constraints,
        &b,
        Vec::with_capacity(0),
        0, 0, nr_upper_bounded_constraints, nr_lower_bounded_constraints,
        &variables,
    );

    type IM = Carry<RationalBig, BasisInverseRows<RationalBig>>;
    let inverse_maintainer = IM::from_basis(&basis_columns, &matrix);

    let mut tableau = Tableau::<_, NonArtificial<_>>::new_with_inverse_maintainer(
        &matrix, inverse_maintainer, basis_columns.into_iter().collect(),
    );
    phase_two::primal::<_, _, FirstProfitable>(&mut tableau);
}
