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
