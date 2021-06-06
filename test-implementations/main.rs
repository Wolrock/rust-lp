fn main() {
    let m1: (usize, Vec<Vec<(usize, f64)>>) = (2, vec![vec![(0, 1.)], vec![(1, 1.)]]);
    let m2: (usize, Vec<Vec<(usize, f64)>>) = (2, vec![vec![(0, 1.), (1, 1.)]]);

    println!("{:?}", multiply_matrices(&m1, &m2));
}
// TODO: Tests
fn change_basis(
    u: &mut Vec<Vec<(usize, f64)>>,
    l: &mut Vec<Vec<(usize, f64)>>,
    u_bar: &mut Vec<Vec<(usize, f64)>>,
    l_bar: &Vec<Vec<(usize, f64)>>,
    q_bar: &Vec<Vec<(usize, f64)>>,
    p_bar: &Vec<Vec<(usize, f64)>>,
    active_block_column: usize,
    active_block_row: usize,
) {
    // TODO: Can probably be done in one loop or inplace
    //  Blocks needed for update
    let mut u_12: (usize, Vec<Vec<(usize, f64)>>) =
        get_block(&u, 1, 2, active_block_row, active_block_column);

    let mut u_23: (usize, Vec<Vec<(usize, f64)>>) =
        get_block(&u, 2, 3, active_block_row, active_block_column);

    let mut u_22: (usize, Vec<Vec<(usize, f64)>>) =
        get_block(&u, 2, 2, active_block_row, active_block_column);

    // Empty vectors dont get inserted
    let mut l_21: (usize, Vec<Vec<(usize, f64)>>) =
        get_block(&l, 2, 1, active_block_row, active_block_column);

    let l_22: (usize, Vec<Vec<(usize, f64)>>) =
        get_block(&l, 2, 2, active_block_row, active_block_column);

    let mut l_32: (usize, Vec<Vec<(usize, f64)>>) =
        get_block(&l, 3, 2, active_block_row, active_block_column);

    let l_22_inv = invert(&l_22);
    let l_bar_inv = invert(&l_bar);
    let p_bar_inv = invert(&p_bar);
    let q_bar_inv = invert(&q_bar);

    // Compute new entries for L matrix
    // TODO: Consider associativity + can probably be done inplace
    l_21 = multiply_matrices(&p_bar_inv, &l_21);

    l_32 = multiply_matrices(
        &multiply_matrices(&l_32, &l_22_inv),
        &multiply_matrices(&p_bar, l_bar),
    );

    // Compute new entries for U matrix
    // TODO: Consider associativity + can probably be done inplace
    u_12 = multiply_matrices(&u_12, &q_bar);
    u_23 = multiply_matrices(
        &multiply_matrices(&l_bar_inv, &p_bar_inv),
        &multiply_matrices(&l_22, &u_23),
    );

    // Restore L and U to upper triangular
    // Could be avoided if blocks use inplace update
    l[active_block_column..=active_block_row].clone_from_slice(l_bar);
    insert_block(&mut l, l_32, active_block_column, active_block_column);
    insert_block(&mut l, &l_32, active_block_column, active_block_column);

    for column in 0..active_block_column {
        l[column].append(&mut l_21[column]);
        l[column].sort_by_key(|&(i, _)| i);
    }

    u[active_block_column..=active_block_row].clone_from_slice(&u_12);
    for column in active_block_column..=active_block_row {
        u[column].append(&mut u_bar[column]);
    }
    for column in active_block_row..u.len() {
        u[column].append(&mut u_23[column]);
        u[column].sort_by_key(|&(i, _)| i);
    }
}

fn multiply_matrices(
    a: &(usize, Vec<Vec<(usize, f64)>>),
    b: &(usize, Vec<Vec<(usize, f64)>>),
) -> (usize, Vec<Vec<(usize, f64)>>) {
    assert_eq!(a.1.len(), b.0);
    let a_size = (a.0, a.1.len());
    let b_size = (b.0, b.1.len());
    let mut product: Vec<Vec<(usize, f64)>> = vec![Vec::new(); b_size.1];
    for column in 0..b_size.1 {
        let mut column_product: Vec<Vec<(usize, f64)>> = Vec::new();
        let mut row_indices: Vec<usize> = vec![];
        for (row, val) in &b.1[column] {
            column_product.push(a.1[*row].iter().map(|&(i, x)| (i, x * val)).collect());
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
            let mut row_sum = 0.;
            let column_indices = (0..column_product.len()).filter_map(|j| {
                column_product[j]
                    .binary_search_by_key(&row, |&(i, _)| i)
                    .ok()
                    .map(|i| (i, j))
            });
            for (i, j) in column_indices {
                row_sum += column_product[j][i].1;
            }
            if row_sum != 0. {
                product[column].push((row, row_sum));
            }
        }
    }
    (a_size.0, product)
}

fn get_block(
    m: &Vec<Vec<(usize, f64)>>,
    row: usize,    // Row of block in block representation of matrix M (exclusive)
    column: usize, // Column of block in block representation of matrix M (exclusive)
    active_block_column: usize, // Inclusive
    active_block_row: usize, // Inclusive
) -> (usize, Vec<Vec<(usize, f64)>>) {
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

    let mut block: (usize, Vec<Vec<(usize, f64)>>) = (
        r2 - r1,
        m[c1..c2]
            .iter()
            .map(|vec| {
                vec.iter()
                    .filter(|&(i, _)| r1 <= *i && *i < r2)
                    .map(|&(i, x)| (i - r1, x))
                    .collect()
            })
            .collect(),
    );
    block
}

// Inserts block into matrix inplace
// Assume block has exclusive column length in first arg of tuple
fn insert_block(
    m: &mut Vec<Vec<(usize, f64)>>,          // Matrix M to insert into
    block: &(usize, Vec<Vec<(usize, f64)>>), // Block to insert
    row: usize,                              // Row where the block begins in M (inclusive)
    column: usize,                           // Column          "
) {
    // Restore correct row indices in relation to M
    let aligned_block: (usize, Vec<Vec<(usize, f64)>>) = (
        block.0,
        (0..block.1.len())
            .map(|j| block.1[j].iter().map(|&(i, x)| (i + row, x)).collect())
            .collect(),
    );

    // Update the row entries belonging to the block in M
    for j in (column..column + aligned_block.1.len()) {
        for i in 0..aligned_block.0 {
            let search = m[j].binary_search_by_key(&aligned_block.1[j - column][i].0, |&(i, _)| i);
            match search {
                Ok(x) => m[j][x] = aligned_block.1[j - column][i],
                Err(x) => m[j].insert(x, aligned_block.1[j - column][i]),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::multiply_matrices;
    use crate::{get_block, insert_block};
    #[test]
    fn multiply_identity() {
        let id: (usize, Vec<Vec<(usize, f64)>>) =
            (3, vec![vec![(0, 1.)], vec![(1, 1.)], vec![(2, 1.)]]);
        let m1: (usize, Vec<Vec<(usize, f64)>>) = (
            3,
            vec![
                vec![(0, 1.), (1, 3.)],
                vec![(0, 2.)],
                vec![(1, 4.), (2, 1.)],
            ],
        );
        assert_eq!(multiply_matrices(&m1, &id), m1);
        assert_eq!(multiply_matrices(&id, &m1), m1);
        assert_eq!(multiply_matrices(&id, &id), id);
    }

    #[test]
    fn multiply_zero() {
        let zero: (usize, Vec<Vec<(usize, f64)>>) = (3, vec![vec![], vec![], vec![]]);
        let m1: (usize, Vec<Vec<(usize, f64)>>) = (
            3,
            vec![
                vec![(0, 1.), (1, 3.)],
                vec![(0, 2.)],
                vec![(1, 4.), (2, 1.)],
            ],
        );
        assert_eq!(multiply_matrices(&m1, &zero), zero);
        assert_eq!(multiply_matrices(&zero, &m1), zero);
        assert_eq!(multiply_matrices(&zero, &zero), zero);
    }

    #[test]
    fn multiply_non_square() {
        let m1: (usize, Vec<Vec<(usize, f64)>>) = (2, vec![vec![(0, 1.)], vec![(1, 1.)]]);
        let m2: (usize, Vec<Vec<(usize, f64)>>) = (2, vec![vec![(0, 1.), (1, 1.)]]);
        let m3: (usize, Vec<Vec<(usize, f64)>>) = (2, vec![vec![(0, 1.), (1, 2.)], vec![(1, 3.)]]);
        assert_eq!(multiply_matrices(&m1, &m2), m2);
        assert_eq!(
            multiply_matrices(&m3, &m2),
            (2, vec![vec![(0, 1.0), (1, 5.0)]])
        );
    }

    #[test]
    fn insert_column() {
        let mut id: Vec<Vec<(usize, f64)>> = vec![vec![(0, 1.)], vec![(1, 1.)], vec![(2, 1.)]];
        let block1: (usize, Vec<Vec<(usize, f64)>>) = (3, vec![vec![(0, 1.), (1, 1.), (2, 1.)]]);

        let result = vec![
            vec![(0, 1.), (1, 1.), (2, 1.)],
            vec![(1, 1.)],
            vec![(2, 1.)],
        ];
        insert_block(&mut id, &block1, 0, 0);
        assert_eq!(result, id);
    }

    #[test]
    fn insert_non_square() {
        let mut id: Vec<Vec<(usize, f64)>> = vec![vec![(0, 1.)], vec![(1, 1.)], vec![(2, 1.)]];
        let block1: (usize, Vec<Vec<(usize, f64)>>) = (
            3,
            vec![
                vec![(0, 1.), (1, 1.), (2, 1.)],
                vec![(0, 1.), (1, 1.), (2, 1.)],
            ],
        );

        let result = vec![
            vec![(0, 1.)],
            vec![(0, 1.), (1, 1.), (2, 1.)],
            vec![(0, 1.), (1, 1.), (2, 1.)],
        ];
        insert_block(&mut id, &block1, 0, 1);
        assert_eq!(result, id);
    }
    #[test]
    fn insert_square() {
        let mut id: Vec<Vec<(usize, f64)>> = vec![vec![(0, 1.)], vec![(1, 1.)], vec![(2, 1.)]];
        let block1: (usize, Vec<Vec<(usize, f64)>>) =
            (2, vec![vec![(0, 2.), (1, 2.)], vec![(0, 2.), (1, 2.)]]);

        let result = vec![
            vec![(0, 1.)],
            vec![(1, 2.), (2, 2.)],
            vec![(1, 2.), (2, 2.)],
        ];
        insert_block(&mut id, &block1, 1, 1);
        assert_eq!(result, id);
    }
    #[test]
    fn insert_element() {
        let mut id: Vec<Vec<(usize, f64)>> = vec![vec![(0, 1.)], vec![(1, 1.)], vec![(2, 1.)]];
        let block1: (usize, Vec<Vec<(usize, f64)>>) = (1, vec![vec![(0, 2.)]]);

        let result = vec![vec![(0, 1.)], vec![(1, 2.)], vec![(2, 1.)]];
        insert_block(&mut id, &block1, 1, 1);
        assert_eq!(result, id);
    }

    #[test]
    fn insert_full_matrix() {
        let mut ones: Vec<Vec<(usize, f64)>> = vec![
            vec![(0, 1.), (1, 1.), (2, 1.)],
            vec![(0, 1.), (1, 1.), (2, 1.)],
            vec![(0, 1.), (1, 1.), (2, 1.)],
        ];
        let block1: (usize, Vec<Vec<(usize, f64)>>) = (
            3,
            vec![
                vec![(0, 2.), (1, 2.), (2, 2.)],
                vec![(0, 2.), (1, 2.), (2, 2.)],
                vec![(0, 2.), (1, 2.), (2, 2.)],
            ],
        );

        let result = vec![
            vec![(0, 2.), (1, 2.), (2, 2.)],
            vec![(0, 2.), (1, 2.), (2, 2.)],
            vec![(0, 2.), (1, 2.), (2, 2.)],
        ];
        insert_block(&mut ones, &block1, 0, 0);
        assert_eq!(result, ones);
    }
    // TODO(Debug): Test for get_block

    #[test]
    fn get_matrix() {
        let id: Vec<Vec<(usize, f64)>> = vec![vec![(0, 1.)], vec![(1, 1.)], vec![(2, 1.)]];
        assert_eq!(get_block(&id, 2, 2, 0, 2), (3, id));
    }

    #[test]
    fn get_active_block() {
        let id: Vec<Vec<(usize, f64)>> = vec![vec![(0, 1.)], vec![(1, 1.)], vec![(2, 1.)]];
        assert_eq!(get_block(&id, 2, 2, 1, 1), (1, vec![vec![(0, 1.)]]));
    }
}
