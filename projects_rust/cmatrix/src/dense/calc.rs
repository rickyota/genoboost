#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
mod simd;

#[allow(unreachable_code)]
#[inline(always)]
pub fn popcnt(x: u32) -> usize {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        return popcnt_simd(x);
    }
    popcnt_nosimd(x)
}

// split for test and this only is used somewhere
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline(always)]
pub fn popcnt_simd(x: u32) -> usize {
    simd::popcnt(x)
}

pub fn popcnt_nosimd(x: u32) -> usize {
    x.count_ones() as usize
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_popcnt_nosimd() {
        let x: u32 = 0b1000_0000_0000_0000_0000_0000_0000_0001;
        assert_eq!(popcnt_nosimd(x), 2);
    }
}
