use std::collections::HashMap;

/* pub fn push<T: Copy>(vec: &mut Vec<T>, val: T, len: usize) {
    for _ in 0..len {
        vec.push(val);
    }
} */

pub fn fill<T: Copy>(vec: &mut [T], val: T) {
    vec.iter_mut().for_each(|v| *v = val);
}

pub fn sort_float<T: PartialOrd>(vec: &mut [T]) {
    vec.sort_by(|a, b| a.partial_cmp(b).unwrap());
}

pub fn sort_and_argsort<T>(vec: &[T]) -> (Vec<T>, Vec<usize>)
where
    T: PartialOrd + Copy,
{
    //let mut v_s = v.to_vec();
    let mut vec_pair: Vec<(usize, T)> = vec.iter().enumerate().map(|(i, v)| (i, *v)).collect();
    vec_pair.sort_by(|(_, vi), (_, vj)| vi.partial_cmp(vj).unwrap());

    let (vec_index, vec_sorted) = vec_pair.iter().cloned().unzip();

    (vec_sorted, vec_index)
}

/*
pub fn sort_and_argsort(vec: &[f64]) -> (Vec<f64>, Vec<usize>) {
    //let mut v_s = v.to_vec();
    let mut vec_pair: Vec<(usize, f64)> = vec.iter().enumerate().map(|(i, v)| (i, *v)).collect();
    vec_pair.sort_by(|(_, vi), (_, vj)| vi.partial_cmp(vj).unwrap());

    let (vec_index, vec_sorted) = vec_pair.iter().cloned().unzip();

    (vec_sorted, vec_index)
}
 */

// https://stackoverflow.com/questions/51272571/how-do-i-check-if-a-slice-is-sorted
pub fn is_sorted<T>(vec: &[T]) -> bool
where
    T: Ord,
{
    vec.windows(2).all(|w| w[0] <= w[1])
}

#[inline]
pub fn count_true(vec: &[bool]) -> usize {
    vec.iter().filter(|&v| *v).count()
}

/// create HashMap<index in, index>: n_in to n
/// return hashmap and n
pub fn create_hashmap_index_from_use(use_vec: &[bool]) -> (HashMap<usize, usize>, usize) {
    let n_in = use_vec.len();
    let mut n_in_to_n = HashMap::with_capacity(n_in);
    let mut ni: usize = 0;
    for (vi, v) in use_vec.iter().enumerate() {
        if *v {
            n_in_to_n.insert(vi, ni);
            ni += 1;
        }
    }
    let n = ni;
    (n_in_to_n, n)
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

    /*
    #[test]
    fn test_with_capacity_align_u8() {
        let len: usize = 5;
        let v = with_capacity_align_u8(len);
        assert_eq!(v.as_ptr() as usize % 256, 0);
        assert_eq!(v.capacity(), len);
        assert_eq!(v.len(), 0);
    }

    #[test]
    #[should_panic]
    fn test_with_capacity_align_u8_panic() {
        let len: usize = 0;
        let v = with_capacity_align_u8(len);
        assert_eq!(v.as_ptr() as usize % 256, 0);
        assert_eq!(v.capacity(), len);
        assert_eq!(v.len(), 0);
    }
    */
}
