/// v is col major
/// col major -> row major
pub fn convert_vec2d_to_row_major(v: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let row_n = v[0].len();
    let mut vals_col = vec![Vec::new(); row_n];
    for v in v.iter() {
        for (x, v_row) in v.iter().zip(vals_col.iter_mut()) {
            v_row.push(*x);
        }
    }
    vals_col
}

/// row major -> col major
pub fn convert_vec2_row_to_col_major(v: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    convert_vec2d_to_row_major(v)
}

/// 2d vec to 1dvec
/// for col and row major
///
/// if v is col major, shape is (col, row)
/// if v is row major, shape is (row, col)
///
/// input: v
/// output: (1d vec, shape of v)
/// ex. v: [[1,10],[2,20],[3,30]] -> ([1,10,2,20,3,30], (3, 2))
pub fn convert_vec2d_to_vec1(v: Vec<Vec<f64>>) -> (Vec<f64>, (usize, usize)) {
    let v1: Vec<f64> = v.iter().flat_map(|x| x.iter()).cloned().collect();
    (v1, (v.len(), v[0].len()))
    //(v1, v[0].len(), v.len())
}

/// col major
pub fn mean_v2(v: &Vec<Vec<f64>>) -> Vec<f64> {
    let mut v_mean = Vec::new();
    for v in v.iter() {
        v_mean.push(mean_vec(v));
        //v_mean.push(v.iter().sum::<f64>() / v.len() as f64);
    }
    v_mean
}

pub fn mean_vec(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len() as f64
}

/// col major
pub fn std_v2(v: &Vec<Vec<f64>>) -> Vec<f64> {
    let mut v_std = Vec::new();

    let means = mean_v2(v);

    for (v, mean) in v.iter().zip(means.iter()) {
        let v_std_i = std_vec(v, *mean);
        //let diff2 = v.iter().map(|x| (x - mean).powi(2)).sum::<f64>();
        //let v_std_i = (diff2 / v.len() as f64).sqrt();
        v_std.push(v_std_i);
    }

    v_std
}

pub fn std_vec(v: &[f64], mean: f64) -> f64 {
    let diff2 = v.iter().map(|x| (x - mean).powi(2)).sum::<f64>();
    let v_std_i = (diff2 / v.len() as f64).sqrt();
    v_std_i
}

pub fn std_vector(v: &[f64]) -> f64 {
    let mean = mean_vec(v);
    std_vec(v, mean)
}

// TODO: use norm_vec()
/// col major
pub fn norm_vec2d(v: Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let means = mean_v2(&v);
    let stds = std_v2(&v);

    let mut v_norm = v;
    //let mut v_norm = v.clone();

    for (row_v, (mean, std)) in v_norm.iter_mut().zip(means.iter().zip(stds.iter())) {
        for x in row_v.iter_mut() {
            *x = (*x - mean) / std;
        }
    }

    v_norm
}

//pub fn norm_vec(v: &[f64], mean: f64, std: f64) -> Vec<f64> {
//    v.iter().map(|x| (x - mean) / std).collect::<Vec<f64>>()
//}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_convert_vec2_col_to_row_major() {
        let v = vec![vec![1.0f64, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let v_ans = vec![vec![1.0f64, 4.0], vec![2.0, 5.0], vec![3.0, 6.0]];
        assert_eq!(convert_vec2d_to_row_major(&v), v_ans);
    }

    #[test]
    fn test_convert_vec2_to_vec() {
        // v: row major
        // ex. v: [[1,10],[2,20],[3,30]] -> ([1,10,2,20,3,30], 3, 2)
        let v = vec![vec![1.0f64, 10.0], vec![2.0, 20.0], vec![3.0, 30.0]];
        let v_ans = vec![1.0f64, 10.0, 2.0, 20.0, 3.0, 30.0];
        assert_eq!(convert_vec2d_to_vec1(v), (v_ans, (3, 2)));
        // col=2, row=3
        //let v = vec![vec![1.0f64, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        //let v_ans = vec![1.0f64, 2.0, 3.0, 4.0, 5.0, 6.0];
        //assert_eq!(convert_vec2d_to_vec1(v), (v_ans, (2, 3)));
    }

    #[test]
    fn test_mean_v2() {
        let v = vec![vec![1.0f64, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let v_ans = vec![2.0f64, 5.0];
        assert_eq!(mean_v2(&v), v_ans);
    }

    #[test]
    fn test_std_v2() {
        let v = vec![vec![1.0f64, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let v_ans = vec![(2.0f64 / 3.0).sqrt(), (2.0f64 / 3.0).sqrt()];
        assert_eq!(std_v2(&v), v_ans);
    }

    #[test]
    fn test_norm_vec2() {
        let v = vec![vec![1.0f64, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let std = (2.0f64 / 3.0).sqrt();
        let v_ans = vec![
            vec![-1.0f64 / std, 0.0, 1.0 / std],
            vec![-1.0 / std, 0.0, 1.0 / std],
        ];
        let normed_v = norm_vec2d(v);
        assert_eq!(normed_v, v_ans);

        let means = mean_v2(&normed_v);
        let stds = std_v2(&normed_v);
        assert_float_absolute_eq!(means[0], 0.0);
        assert_float_absolute_eq!(means[1], 0.0);
        assert_float_absolute_eq!(stds[0], 1.0);
        assert_float_absolute_eq!(stds[1], 1.0);
    }
}
