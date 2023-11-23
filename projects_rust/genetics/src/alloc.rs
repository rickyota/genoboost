use maligned::{align_first, A256};
use sysinfo::{RefreshKind, SystemExt};

#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
//use std::convert::TryInto;

/// return aligned vector to 256 byte
/// capacity() is len but len() is 0
pub fn with_capacity_align_u8(capacity: usize) -> Vec<u8> {
    if capacity == 0 {
        panic!("Indicate capacity>0.")
    }
    let v: Vec<u8> = align_first::<u8, A256>(capacity);
    assert_eq!(v.as_ptr() as usize % 256, 0);
    assert_eq!(v.capacity(), capacity);
    assert_eq!(v.len(), 0);
    v
}

pub fn with_capacity_align_f64(capacity: usize) -> Vec<f64> {
    if capacity == 0 {
        panic!("Indicate capacity>0.")
    }
    let v: Vec<f64> = align_first::<f64, A256>(capacity);
    assert_eq!(v.as_ptr() as usize % 256, 0);
    assert_eq!(v.capacity(), capacity);
    assert_eq!(v.len(), 0);
    v
}

/*
// cannot call this func
pub fn with_capacity_align<T>(capacity: usize, dtype: T) -> Vec<T> {
    if capacity == 0 {
        panic!("Indicate capacity>0.")
    }
    let v: Vec<T> = align_first::<T, A256>(capacity);
    assert_eq!(v.as_ptr() as usize % 256, 0);
    assert_eq!(v.capacity(), capacity);
    assert_eq!(v.len(), 0);
    v
}
*/

#[target_feature(enable = "avx2")]
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
pub unsafe fn with_capacity_align_m256i(capacity: usize) -> Vec<__m256i> {
    if capacity == 0 {
        panic!("Indicate capacity>0.")
    }
    let v: Vec<__m256i> = align_first::<__m256i, A256>(capacity);
    assert_eq!(v.as_ptr() as usize % 256, 0);
    assert_eq!(v.capacity(), capacity);
    assert_eq!(v.len(), 0);

    v
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

/* // should create unused_memory()?
pub fn get_total_memory_old() -> usize {
    let mut system = sysinfo::System::new_all();
    system.refresh_all();

    1024 * system.total_memory() as usize
} */

#[cfg(test)]
mod tests {
    use super::*;

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
}
