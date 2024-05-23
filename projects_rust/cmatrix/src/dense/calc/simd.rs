#[cfg(target_arch = "x86_64")]
#[inline(always)]
pub fn popcnt(x: u32) -> usize {
    unsafe { core::arch::x86_64::_popcnt32(x as i32) as usize }
    //unsafe { core::arch::x86_64::_popcnt64(x as i64) as usize }
}

#[cfg(target_arch = "x86")]
#[inline(always)]
pub fn popcnt(x: u32) -> usize {
    unsafe { core::arch::x86::_popcnt32(x as i32) as usize }
    //unsafe { core::arch::x86::_popcnt32(x as i64) as usize }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn test_popcnt_nosimd() {
        let x: u32 = 0b1000_0000_0000_0000_0000_0000_0000_0001;
        assert_eq!(popcnt(x), 2);
    }

    #[cfg(target_arch = "x86")]
    #[test]
    fn test_popcnt_nosimd() {
        let x: u32 = 0b1000_0000_0000_0000_0000_0000_0000_0001;
        assert_eq!(popcnt(x), 2);
    }
}
