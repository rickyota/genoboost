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


pub fn or_bool_vec_mut(v1: &mut [bool], v2: &[bool]) {
    v1.iter_mut().zip(v2.iter()).for_each(|(b1, b2)| *b1=*b1 | *b2);
}


pub fn and_bool_vec_mut(v1: &mut [bool], v2: &[bool]) {
    v1.iter_mut().zip(v2.iter()).for_each(|(b1, b2)| *b1=*b1 & *b2);
}

pub fn and_bool_vec(v1: &[bool], v2: &[bool]) -> Vec<bool> {
    v1.iter().zip(v2.iter()).map(|(b1, b2)| b1 & b2).collect()
}

// TODO: cleaner
/// v2.len() = v1.count_true()
pub fn and_in_bool_vec(v1: &[bool], v2: &[bool]) -> Vec<bool> {
    assert_eq!(count_true(&v1), v2.len());
    //v1.iter().zip(v2.iter()).map(|(b1, b2)| b1 & b2).collect()

    let mut v = vec![false; v1.len()];

    let mut i2 = 0usize;
    for (i1, &x1) in v1.iter().enumerate() {
        if x1 {
            v[i1] = x1 & v2[i2];
            i2 += 1;
        }
    }

    v
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
