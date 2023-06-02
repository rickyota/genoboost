use maligned::{align_first, A256};

// recover if necessary
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
