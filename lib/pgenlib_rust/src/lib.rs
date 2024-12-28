#[cxx::bridge]
pub mod pgenlib_ffi {
    unsafe extern "C++" {

        // "test_cbind": the name of package defined in Cargo.toml
        include!("pgenlib_rust/lib/pgenlib/pgenlibr_wrapc.hpp");
        //include!("test_cbind_rust30/include/blobstore.h");
        //
        // cannot use '.' or '..' so path should in the crate
        // > #include relative to `.` or `..` is not supported in Cargo builds
        // > note: use a path starting with the crate name
        //include!("lib/pgenlib_rust/src/pgenlibr_wrapc.hpp");
        //include!("genetics/lib/pgenlib_rust/src/pgenlibr_wrapc.hpp");
        //include!("projects_rust/genetics/lib/pgenlib_rust/src/pgenlibr_wrapc.hpp");
        //include!("../../lib/pgenlib_rust/src/pgenlibr_wrapc.h");
        //include!("/nfs/data06/ricky/code/genetics/projects_rust/genetics/lib/pgenlib_rust/src/pgenlibr_wrapc.hpp")

        fn foo();
        fn gcd(a: i32, b: i32) -> i32;
        // https://stackoverflow.com/questions/74395595/cxx-translate-const-char-to-rust-equivalent
        unsafe fn ghi(a: *const c_char) -> i32;
        //unsafe fn ghi(a: *const i8) -> i32;
        //fn ghi(a: &str) -> i32;
        //fn foo_bar();

        unsafe fn pgenreader_load_snvs_extract(
            genot: *mut i8,
            filename: *const c_char,
            snv_start: i32,
            snv_end: i32,
            use_snvs: *const bool,
            sample_n_in: i32,
            sample_subset: *const i32,
            sample_subset_n: i32,
            nthr: i32,
        ) -> i32;

        // int pgenreader_load(double *genot, const char *filename, int sample_n, int *sample_subset, int sample_subset_n, int nthr);
        //unsafe fn pgenreader_load(genot: *mut f64,
        //unsafe fn pgenreader_load(
        //    genot: *mut f64,
        //    filename: *const c_char,
        //    sample_n: i32,
        //    sample_subset: *mut i32,
        //    sample_subset_n: i32,
        //    nthr: i32,
        //) -> i32;

        //unsafe fn pgenreader_load_arg(
        //    //filename: *const c_char,
        //    nthr: i32);

    }
}
