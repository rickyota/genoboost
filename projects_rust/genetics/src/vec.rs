pub mod vec2;

pub use vec2::*;

use std::collections::HashSet;
use std::hash::Hash;

pub fn fill<T: Copy>(vec: &mut [T], val: T) {
    vec.iter_mut().for_each(|v| *v = val);
}

pub fn sort_float<T: PartialOrd>(vec: &mut [T]) {
    vec.sort_by(|a, b| a.partial_cmp(b).unwrap());
}

pub fn sort_and_argsort<T: PartialOrd + Copy>(vec: &[T]) -> (Vec<T>, Vec<usize>) {
    let mut vec_pair: Vec<(usize, T)> = vec.iter().enumerate().map(|(i, v)| (i, *v)).collect();
    vec_pair.sort_by(|(_, vi), (_, vj)| vi.partial_cmp(vj).unwrap());

    let (vec_index, vec_sorted) = vec_pair.iter().cloned().unzip();

    (vec_sorted, vec_index)
}

// is_sorted() is nightly
// https://stackoverflow.com/questions/51272571/how-do-i-check-if-a-slice-is-sorted
pub fn is_sorted<T: Ord>(vec: &[T]) -> bool {
    vec.windows(2).all(|w| w[0] <= w[1])
}

#[inline]
pub fn count_true(vec: &[bool]) -> usize {
    vec.iter().filter(|&v| *v).count()
}

pub fn is_all_same_float<T: PartialOrd>(vec: &[T]) -> bool {
    vec.windows(2).all(|w| w[0] == w[1])
}

// not using now
/// create HashMap<index in, index>: n_in to n
/// return hashmap and n
//pub fn create_hashmap_index_from_use(use_vec: &[bool]) -> (HashMap<usize, usize>, usize) {
//    let n_in = use_vec.len();
//    let mut n_in_to_n = HashMap::with_capacity(n_in);
//    let mut ni: usize = 0;
//    for (vi, v) in use_vec.iter().enumerate() {
//        if *v {
//            n_in_to_n.insert(vi, ni);
//            ni += 1;
//        }
//    }
//    let n = ni;
//    (n_in_to_n, n)
//}

pub fn or_bool_vec_mut(v1: &mut [bool], v2: &[bool]) {
    v1.iter_mut()
        .zip(v2.iter())
        .for_each(|(b1, b2)| *b1 = *b1 || *b2);
}

pub fn or_bool_vec(v1: &[bool], v2: &[bool]) -> Vec<bool> {
    v1.iter()
        .zip(v2.iter())
        .map(|(&b1, &b2)| b1 || b2)
        .collect()
}

pub fn and_bool_vec_mut(v1: &mut [bool], v2: &[bool]) {
    v1.iter_mut()
        .zip(v2.iter())
        .for_each(|(b1, b2)| *b1 = *b1 && *b2);
}

pub fn and_bool_vec(v1: &[bool], v2: &[bool]) -> Vec<bool> {
    v1.iter()
        .zip(v2.iter())
        .map(|(&b1, &b2)| b1 && b2)
        .collect()
}

/// v2.len() = v1.count_true()
pub fn and_in_bool_vec(v1: &[bool], v2: &[bool]) -> Vec<bool> {
    assert_eq!(count_true(&v1), v2.len());
    //v1.iter().zip(v2.iter()).map(|(b1, b2)| b1 && b2).collect()

    let mut v = vec![false; v1.len()];

    let mut i2 = 0usize;
    for (i1, &x1) in v1.iter().enumerate() {
        if x1 {
            v[i1] = x1 && v2[i2];
            i2 += 1;
        }
    }

    //for (i1, x) in v.iter_mut().enumerate() {
    //    let x1 = v1[i1];
    //    if x1 {
    //        *x = x1 && v2[i2];
    //        i2 += 1;
    //    }
    //}

    v
}

// not tested
// pub fn is_unique<T: Eq + Hash>(v: &[T]) -> bool {
//     let mut uniq = HashSet::new();
//     v.iter().all(move |x| uniq.insert(x))
// }

pub fn uniq_clone<T: Clone + Eq + Hash>(v: &[T]) -> HashSet<T> {
    v.into_iter()
        .map(|x| x.clone())
        .collect::<HashSet<T>>()
        .into_iter()
        .collect()
}

// pub fn uniq_clone(v: &[String]) -> HashSet<String> {
//     v.into_iter()
//         .map(|x| x.clone())
//         .collect::<HashSet<String>>()
//         .into_iter()
//         .collect()
// }

//pub fn convert_type_string<T: std::str::FromStr>(v: Vec<String>) -> Vec<T>
pub fn convert_type_string<T: std::str::FromStr>(v: &[String]) -> Vec<T>
where
    <T as std::str::FromStr>::Err: std::fmt::Debug,
{
    v.into_iter()
        .map(|x| x.parse::<T>().unwrap())
        .collect::<Vec<T>>()
}

// similar to Vec::extract_if() deinfed in nightly
pub fn extract_if_iter<'a, T, I>(v: I, f: &[bool]) -> Vec<T>
where
    I: Iterator<Item = &'a T>,
    T: Clone + 'a,
{
    v.zip(f.iter())
        .filter(|(_, &b)| b)
        .map(|(x, _)| x.clone())
        .collect::<Vec<T>>()
}

pub fn extract_if_into_iter<'a, T, I>(v: I, f: &[bool]) -> Vec<T>
where
    I: IntoIterator<Item = T>,
{
    v.into_iter()
        .zip(f.iter())
        .filter(|(_, &b)| b)
        .map(|(x, _)| x)
        .collect::<Vec<T>>()
}

//fn check_if_unique(v: &[String]){
pub fn has_unique_elements<T>(iter: T) -> bool
where
    T: IntoIterator,
    T::Item: Eq + Hash,
{
    let mut uniq = HashSet::new();
    iter.into_iter().all(move |x| uniq.insert(x))
}

pub fn histogram(v: &[usize], bin_n: usize) -> (Vec<usize>, Vec<usize>) {
    // create rough histogram

    let mut hist = vec![0; bin_n];
    // ex. max=103, n=5
    // bin_width = 20
    let bin_width = v.iter().max().unwrap() / bin_n;
    // (20,40,60,80); len=bin_n-1
    let bins = (1..bin_n).map(|i| i * bin_width).collect::<Vec<usize>>();
    for &x in v.iter() {
        if x < bins[0] {
            hist[0] += 1;
        } else if x >= bins[bin_n - 2] {
            hist[bin_n - 1] += 1;
        } else {
            for i in 0..(bin_n - 1) {
                if (x >= bins[i]) && (x < bins[i + 1]) {
                    hist[i + 1] += 1;
                    break;
                }
            }
        }
    }
    (hist, bins)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sort_and_argsort() {
        let v = [1.0, 5.0, 2.0, 4.0, 3.0];
        let v_sorted_exp = [1.0, 2.0, 3.0, 4.0, 5.0];
        let v_index_exp = [0, 2, 4, 3, 1];
        let (v_sorted, v_index) = sort_and_argsort(&v);
        assert_eq!(v_sorted, v_sorted_exp);
        assert_eq!(v_index, v_index_exp);
    }

    #[test]
    fn test_is_sorted() {
        let v = [0usize, 1, 2];
        assert_eq!(is_sorted(&v), true);
    }

    #[test]
    fn test_is_sorted_false() {
        let v = [1usize, 0, 2];
        assert_eq!(is_sorted(&v), false);
    }

    #[test]
    fn test_count_true() {
        let v = [true, true, false];
        assert_eq!(count_true(&v), 2);
    }

    //#[test]
    //fn test_create_hashmap_index_from_use(){}

    #[test]
    fn test_or_bool_vec_mut() {
        let mut v1 = [true, true, false, false];
        let v2 = [true, false, true, false];
        let v_ans = [true, true, true, false];
        or_bool_vec_mut(&mut v1, &v2);
        assert_eq!(v1, v_ans);
    }

    #[test]
    fn test_and_bool_vec_mut() {
        let mut v1 = [true, true, false, false];
        let v2 = [true, false, true, false];
        let v_ans = [true, false, false, false];
        and_bool_vec_mut(&mut v1, &v2);
        assert_eq!(v1, v_ans);
    }

    #[test]
    fn test_and_bool_vec() {
        let v1 = [true, true, false, false];
        let v2 = [true, false, true, false];
        let v_ans = [true, false, false, false];
        let v = and_bool_vec(&v1, &v2);
        assert_eq!(v, v_ans);
    }

    #[test]
    fn test_and_in_bool_vec() {
        let v1 = [false, true, true, false];
        let v2 = [false, true];
        let v_ans = [false, false, true, false];
        let v = and_in_bool_vec(&v1, &v2);
        assert_eq!(v, v_ans);
    }

    #[test]
    fn test_vals_convert_type() {
        let vals_str: Vec<String> = vec!["1.0".to_string(), "1.5".to_string(), "-0.5".to_string()];

        let vals_ans: Vec<f64> = vec![1.0f64, 1.5, -0.5];

        let vals = convert_type_string::<f64>(&vals_str);

        assert_eq!(vals, vals_ans);
    }

    #[test]
    fn test_extract_if_iter() {
        let v1 = [1, 2, 3, 4, 5];
        let v2 = [true, false, true, false, true];
        let v_ans = [1, 3, 5];
        let v = extract_if_iter(v1.iter(), &v2);
        assert_eq!(v, v_ans);
    }

    #[test]
    fn test_extract_if_into_iter() {
        let v1 = [
            "1".to_string(),
            "2".to_string(),
            "3".to_string(),
            "4".to_string(),
            "5".to_string(),
        ]; //["1", "2", "3", "4", "5"];
        let v2 = [true, false, true, false, true];
        let v_ans = ["1".to_string(), "3".to_string(), "5".to_string()];
        let v = extract_if_iter(v1.iter(), &v2);
        assert_eq!(v, v_ans);
    }

    #[test]
    fn test_historam() {
        let v = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let bin_n = 5;
        let (hist, bins) = histogram(&v, bin_n);
        assert_eq!(bins, [2, 4, 6, 8]);
        let hist_exp = [2, 2, 2, 2, 3];
        assert_eq!(hist, hist_exp);
    }
}
