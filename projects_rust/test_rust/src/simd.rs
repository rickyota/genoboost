#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
unsafe fn test_simd_sse() -> f64 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;
    let _zerod: __m128d = _mm_setzero_pd();
    1.0
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
unsafe fn test_simd_avx() -> f64 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;
    let _zerod: __m256d = _mm256_setzero_pd();
    1.0
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
pub unsafe fn test() {
    unsafe {
        let t = test_simd_sse();
        println!("t: {}", t);

        let s = test_simd_avx();
        println!("s: {}", s);
    }
}
