use maligned::{self, A256};
use sysinfo::{RefreshKind, SystemExt};

//[#[cfg(target_arch = "x86")]
//[use std::arch::x86::*;
//[#[cfg(target_arch = "x86_64")]
//[use std::arch::x86_64::*;

/// return aligned vector to 256 byte
/// capacity() is len but len() is 0
pub fn with_capacity_align<T>(capacity: usize) -> Vec<T> {
    // capacity cannot be 0 for align_first
    let cap = if capacity == 0 { 1 } else { capacity };
    // could be 0 for empty vec
    //if capacity == 0 {
    //    panic!("Indicate capacity>0.")
    //}
    let v: Vec<T> = maligned::align_first::<T, A256>(cap);
    assert_eq!(v.as_ptr() as usize % 256, 0);
    // for capacity=0, capacity is 1
    assert!(v.capacity() >= cap);
    //assert_eq!(v.capacity(), capacity);
    assert_eq!(v.len(), 0);
    v
}

/// return aligned vec of input vec
/// remaining elements are 0.0
// TODO: use num-traits and T: Num::Zero
pub fn vec_align_f64(v: Vec<f64>, capacity: usize) -> Vec<f64> {
    if v.len() > capacity {
        panic!("vec len should be the same or smaller than capacity")
    }
    let mut v_align = with_capacity_align::<f64>(capacity);
    v_align.resize(capacity, 0.0f64);

    v_align
        .iter_mut()
        .zip(v.into_iter())
        .for_each(|(x, y)| *x = y);

    v_align
}

pub fn get_available_memory() -> Option<usize> {
    if sysinfo::System::IS_SUPPORTED {
        // hang up here
        //let system = sysinfo::System::new_all();
        // worked!
        let system = sysinfo::System::new_with_specifics(RefreshKind::new().with_memory());

        Some(system.available_memory() as usize)
    } else {
        None
    }
}

/// memory as GB
pub fn mem_gb(mem: usize) -> f64 {
    (mem as f64) / ((1024usize * 1024 * 1024) as f64)
}


#[cfg(test)]
mod tests {
    use super::*;

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

    #[test]
    fn test_with_capacity_align_usize() {
        let len: usize = 5;
        let v = with_capacity_align::<usize>(len);
        assert_eq!(v.as_ptr() as usize % 256, 0);
        assert_eq!(v.capacity(), len);
        assert_eq!(v.len(), 0);
    }

    #[test]
    fn test_with_capacity_align_usize_0() {
        let len: usize = 0;
        let v: Vec<usize> = with_capacity_align::<usize>(len);
        assert_eq!(v.capacity(), 1);
    }

    #[test]
    fn test_with_capacity_align_f64() {
        let len: usize = 5;
        let v: Vec<f64> = with_capacity_align::<f64>(len);
        assert_eq!(v.as_ptr() as usize % 256, 0);
        assert_eq!(v.capacity(), len);
        assert_eq!(v.len(), 0);
    }

    #[test]
    fn test_mem_gb() {
        let mem: usize = 2 * 1024 * 1024 * 1024;
        assert_eq!(mem_gb(mem), 2.0f64);
    }
}
