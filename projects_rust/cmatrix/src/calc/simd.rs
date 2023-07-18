#[cfg(target_arch = "x86_64")]
#[inline(always)]
pub fn popcnt(x: u32) -> usize {
	//pub fn popcnt(x: u32) -> usize {
    unsafe { core::arch::x86_64::_popcnt32(x as i32) as usize }
    //unsafe { core::arch::x86_64::_popcnt64(x as i64) as usize }
}

//#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
// -> already above
#[cfg(target_arch = "x86")]
#[inline(always)]
pub fn popcnt(x: u32) -> usize {
    unsafe { core::arch::x86::_popcnt32(x as i32) as usize }
    //unsafe { core::arch::x86::_popcnt32(x as i64) as usize }
}