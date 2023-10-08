
#include "pgenlibr_wrapc.h"
#include "pgenlibr.h"
// how to avoid include .cpp
#include "pgenlibr.cpp"
// #include <stdio.h>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

void foo() { printf("Hello, pgenlib!!\n"); }
int gcd(int a, int b) {
    int t;
    while(b != 0) {
        t = a % b;
        a = b;
        b = t;
    }
    return a;
}

int ghi(const char *str) { return 89; }

/*
std::string ghi(){
    std::string str("ghi");
    return str;
}
*/

MyClass::MyClass() { abc = 67; }
int MyClass::get_abc() { return abc; }
MyClass::~MyClass() {}

MyClass *myclass_create() { return new MyClass(); }
int myclass_get_abc(MyClass *m) { return m->get_abc(); }
void myclass_delete(MyClass *m) {}
/*
// size_t pgenreader_load_whole(double *genot, const char *filename, int snv_n_in, int sample_n_in, int *sample_subset, int sample_subset_n, int nthr)
// {
int pgenreader_load_whole(double *genot, const char *filename, int snv_n_in, int sample_n_in, int *sample_subset, int sample_subset_n, int nthr) {
    PgenReader *pread = new PgenReader();
    std::string fname(filename);

    // might have better way
    std::vector<int> sample_subset_1based_v(sample_subset_n);
    std::copy(sample_subset, sample_subset + sample_subset_n, sample_subset_1based_v.begin());
    for(int &d : sample_subset_1based_v) {
        d += 1;
    }

    pread->Load(fname, sample_n_in, sample_subset_1based_v, nthr);

    size_t mi;
#ifdef _OPENMP
    omp_set_num_threads(nthr);
#endif

#ifdef _OPENMP
    size_t const thrn = omp_get_max_threads();
#else
    size_t const thrn = 1;
#endif
    double **buf_thrs = static_cast<double **>(malloc(sizeof(double *) * thrn));
    for(size_t thri = 0; thri < thrn; thri++) {
        double *buf_thr = static_cast<double *>(malloc(sizeof(double) * sample_subset_n));
        buf_thrs[thri] = buf_thr;
    }

#ifdef _OPENMP
#pragma omp parallel for private(mi)
#endif
    for(mi = 0; mi < snv_n_in; mi++) {
#ifdef _OPENMP
        size_t thri = omp_get_thread_num();
#else
        size_t thri = 0;
#endif

        double *buf_thr = buf_thrs[thri];

        pread->Read(buf_thr, sample_subset_n, thri, mi, 1);
        std::copy(buf_thr, buf_thr + sample_subset_n, genot + mi * sample_subset_n);

        // pread->Read(genot + mi * sample_subset_n, sample_subset_n, thread_num, mi, 1);
    }

    // run in parallel
    // [ref](https://github.com/rgcgithub/regenie/blob/master/src/Geno.cpp#L3748)
    // pread->Read(genot, sample_n, 0, 0, 1);
    //  pread->ReadHardcalls(genot, sample_n, 0, 1, 1);
    return nthr;
} */
/*
// TODO: use_snvs; this can speed up a lot by omitting Read()
// size_t pgenreader_load_snvs(double *genot, const char *filename, size_t snv_start, size_t snv_end, size_t sample_n_in, int *sample_subset,
//                            size_t sample_subset_n, size_t nthr) {
int pgenreader_load_snvs(double *genot, const char *filename, int snv_start, int snv_end, int sample_n_in, int const *const sample_subset,
                         int sample_subset_n, int nthr) {
    PgenReader *pread = new PgenReader();
    std::string fname(filename);

    // might have better way
    std::vector<int> sample_subset_1based_v(sample_subset_n);
    std::copy(sample_subset, sample_subset + sample_subset_n, sample_subset_1based_v.begin());
    for(int &d : sample_subset_1based_v) {
        d += 1;
    }
    // printf("sample_subset_1based %d\n", sample_subset_1based_v[0]);

    printf("cpp nthr %d\n", nthr);
    pread->Load(fname, sample_n_in, sample_subset_1based_v, nthr);


    size_t mi;
#ifdef _OPENMP
    //  disable dynamic
    //  https://stackoverflow.com/questions/11095309/openmp-set-num-threads-is-not-working
    // omp_set_dynamic(0);
    omp_set_num_threads(nthr);
#endif

#ifdef _OPENMP
    size_t const thrn = omp_get_max_threads();
#else
    size_t const thrn = 1;
#endif
    double **buf_thrs = static_cast<double **>(malloc(sizeof(double *) * thrn));
    for(size_t thri = 0; thri < thrn; thri++) {
        double *buf_thr = static_cast<double *>(malloc(sizeof(double) * sample_subset_n));
        buf_thrs[thri] = buf_thr;
    }

#ifdef _OPENMP
#pragma omp parallel for private(mi)
#endif
    for(mi = 0; mi < snv_end - snv_start; mi++) {
#ifdef _OPENMP
        size_t thri = omp_get_thread_num();
#else
        size_t thri = 0;
#endif

        double *buf_thr = buf_thrs[thri];
        // should use ReadHardcalls?
        pread->Read(buf_thr, sample_subset_n, thri, mi + snv_start, 1);
        std::copy(buf_thr, buf_thr + sample_subset_n, genot + mi * sample_subset_n);

        // pread->Read(genot + mi * sample_subset_n, sample_subset_n, thread_num, mi + snv_start, 1);
    }

    return nthr;
} */

// genot double* : #snvs(true in use_snvs) * sample_subset_n
int pgenreader_load_snvs_extract(int8_t *genot, const char *filename, int snv_start, int snv_end, bool const *const use_snvs, int sample_n_in,
                                 int const *const sample_subset, int sample_subset_n, int nthr) {

    PgenReader *pread = new PgenReader();
    std::string fname(filename);

    // might have better way
    std::vector<int> sample_subset_1based_v(sample_subset_n);
    std::copy(sample_subset, sample_subset + sample_subset_n, sample_subset_1based_v.begin());
    for(int &d : sample_subset_1based_v) {
        d += 1;
    }

    size_t const snv_len = snv_end - snv_start;

    // TODO: use map
    // std::unordered_map<size_t, size_t> snv_index;
    size_t *snv_index = static_cast<size_t *>(malloc(sizeof(size_t) * snv_len));
    size_t use_index = 0;
    for(size_t mi = 0; mi < snv_len; mi++) {
        if(use_snvs[mi]) {
            snv_index[mi] = use_index;
            use_index++;
        } else {
            // otherwise not initialized
            snv_index[mi] = SIZE_MAX;
        }
    }

    //printf("use snvs %zu\n", use_index);

    //printf("cpp nthr %d\n", nthr);
    pread->Load(fname, sample_n_in, sample_subset_1based_v, nthr);

    size_t mi;
#ifdef _OPENMP
    //  disable dynamic
    //  https://stackoverflow.com/questions/11095309/openmp-set-num-threads-is-not-working
    // omp_set_dynamic(0);
    omp_set_num_threads(nthr);
#endif

#ifdef _OPENMP
    size_t const thrn = omp_get_max_threads();
#else
    size_t const thrn = 1;
#endif
    double **buf_thrs = static_cast<double **>(malloc(sizeof(double *) * thrn));
    for(size_t thri = 0; thri < thrn; thri++) {
        buf_thrs[thri] = static_cast<double *>(malloc(sizeof(double) * sample_subset_n));
    }

#ifdef _OPENMP
#pragma omp parallel for private(mi)
#endif
    for(mi = 0; mi < snv_len; mi++) {
        // for(mi = 0; mi < snv_end - snv_start; mi++) {

        if(!use_snvs[mi]) {
            continue;
        }

#ifdef _OPENMP
        size_t const thri = omp_get_thread_num();
#else
        size_t const thri = 0;
#endif

        double *buf_thr = buf_thrs[thri];
        // should use ReadHardcalls?
        pread->Read(buf_thr, sample_subset_n, thri, mi + snv_start, 1);

        size_t const use_index = snv_index[mi];
        // printf("mi, use_index %zu, %zu\n", mi, use_index);

        for(size_t samplei = 0; samplei < sample_subset_n; samplei++) {
            genot[use_index * sample_subset_n + samplei] = (int8_t)buf_thr[samplei];
        }

        // pread->Read(genot + mi * sample_subset_n, sample_subset_n, thread_num, mi + snv_start, 1);
    }

    return nthr;
}

// Problem was due to int overflow
// mi * sample_n was overflowed
// int pgenreader_load_snvs_debug(double *genot, const char *filename, int snv_start, int snv_end, int sample_n_in, int *sample_subset, int
// sample_subset_n,
//                         int nthr) {
//    PgenReader *pread = new PgenReader();
//    std::string fname(filename);
//
//    // might have better way
//    std::vector<int> sample_subset_1based_v(sample_subset_n);
//    std::copy(sample_subset, sample_subset + sample_subset_n, sample_subset_1based_v.begin());
//    for(int &d : sample_subset_1based_v) {
//        d += 1;
//    }
//    // printf("sample_subset_1based %d\n", sample_subset_1based_v[0]);
//
//    printf("cpp nthr %d\n", nthr);
//    pread->Load(fname, sample_n_in, sample_subset_1based_v, nthr);
//
//    /*     // no error
//        int thread_num = 0;
//        int mi;
//        for(mi = 0; mi < snv_end - snv_start; mi++) {
//            pread->Read(genot + mi * sample_subset_n, sample_subset_n, thread_num, mi + snv_start, 1);
//        } */
//
//    // still error even for nthread=1
//    // int num_threads = 1;
//    size_t mi;
//    // TODO: create Vec<Vec<f64>> for buf
//    // omp_set_dynamic(0);
//    printf("bfr in openmp %d\n", nthr);
// #ifdef _OPENMP
//    //  disable dynamic
//    //  https://stackoverflow.com/questions/11095309/openmp-set-num-threads-is-not-working
//    omp_set_dynamic(0);
//    printf("in openmp %d\n", nthr);
//    printf("cpp nthr %d\n", nthr);
//    // value of OMP_NUM_THREADS=8
//    printf("bfr set omp max #threads %d\n", omp_get_max_threads());
//    // printf("bfr set omp #threads %d\n", omp_get_num_threads());
//    omp_set_num_threads(nthr);
//    printf("omp max #threads %d\n", omp_get_max_threads());
//    // always 1 since outside of parallel??
//    // printf("omp #threads %d\n", omp_get_num_threads());
// #endif
//
//    // THIS WORKED!!
//    //    int mtmp;
//    // #ifdef _OPENMP
//    // #pragma omp parallel for
//    // #endif
//    //    for(mtmp = 0; mtmp < 2; mtmp++) {
//    //        int thread_num = omp_get_thread_num();
//    //        printf("mtmp omp thread# %d / %d\n", thread_num, omp_get_num_threads());
//    //    }
//
//    // ok
//    //    int mitest;
//    // #ifdef _OPENMP
//    // #pragma omp parallel for
//    // #endif
//    //    for(mitest = 0; mitest < snv_end - snv_start; mitest++) {
//    //        int thread_num = omp_get_thread_num();
//    //        printf("mitest omp thread# %d / %d\n", thread_num, omp_get_num_threads());
//    //    }
//
//    // #ifdef _OPENMP
//    // size_t const thrn = omp_get_max_threads();
//    //    size_t const thrn = omp_get_num_threads();
//    // #else
//    //    size_t const thrn = 1;
//    // #endif
//    size_t thrn = nthr;
//    //printf("num thrn %d\n", omp_get_num_threads());
//    //printf("max thrn %d\n", omp_get_max_threads());
//    //printf("thrn %d\n", thrn);
//    // correct?
//    double **buf_thrs = static_cast<double **>(malloc(sizeof(double *) * thrn));
//    for(size_t thri = 0; thri < thrn; thri++) {
//        double *buf_thr = static_cast<double *>(malloc(sizeof(double) * sample_subset_n));
//        buf_thrs[thri] = buf_thr;
//        //printf("buf_thr buf_thr %d, %f\n", thri, buf_thr[0]);
//    }
//
//    // this is ok!
////    int mitest;
////#ifdef _OPENMP
////#pragma omp parallel for private(mitest)
////#endif
////    for(mitest = 0; mitest < 100; mitest++) {
////        int thread_num = omp_get_thread_num();
////        printf("mitest omp thread# %d / %d\n", thread_num, omp_get_num_threads());
////        // ok
////        // double *buf_thr = static_cast<double *>(malloc(sizeof(double) * sample_subset_n));
////        // ng
////        // size_t thri = omp_get_thread_num();
////        // double *buf_thr = buf_thrs[thri];
////        // printf("mitest buf_thr %f\n", buf_thr[0]);
////        // ok ; why ?
////        // std::copy(buf_thr, buf_thr + sample_subset_n, genot + mitest * sample_subset_n);
////    }
////
////    printf("mi max %d\n", snv_end - snv_start);
//
//// this section makes abortion when #pragma is on
// #ifdef _OPENMP
// #pragma omp parallel for private(mi)
// #endif
//     for(mi = 0; mi < snv_end - snv_start; mi++) {
//         int thread_num = 0;
// #ifdef _OPENMP
//         thread_num = omp_get_thread_num();
//         // if(mi % 10000 == 0) {
//         //if(mi < 20) {
//         //    printf("omp thread# %d / %d\n", thread_num, omp_get_num_threads());
//         //}
//         // printf("omp thread# %d / %d\n", thread_num, omp_get_num_threads());
// #endif
//         // printf("args %d, %d, %d, %d\n", sample_subset_n, thread_num, mi + snv_start, 1);
//
//         // ok; no omp error
//         // double *buf_thr = static_cast<double *>(malloc(sizeof(double) * sample_subset_n));
//         // ok; no omp error
//         size_t thri = omp_get_thread_num();
//         double *buf_thr = buf_thrs[thri];
//         pread->Read(buf_thr, sample_subset_n, thread_num, mi + snv_start, 1);
//
//         // pread->Read(genot + mi * sample_subset_n, sample_subset_n, thread_num, mi + snv_start, 1);
//
//         // ng
//         std::copy(buf_thr, buf_thr + sample_subset_n, genot + mi * sample_subset_n);
//         // if(mi < 20) {
//         //if(thri == 1) {
//         //    printf("thri buf_thr %d, %f\n", thri, buf_thr[0]);
//         //}
//         //}
//
//         // ok
//         // if(thri == 0) {
//         //    printf("thr0 copying mi  %d\n", mi);
//         //    std::copy(buf_thr, buf_thr + sample_subset_n, genot + mi * sample_subset_n);
//         //}
//
//         // ng -> due to int overflow
//         //  mi, sample_subset_n was int, so overflowed
//         //if(thri == 1) {
//         //    printf("thr1 copying mi  %d\n", mi);
//         //    // std::copy(buf_thr, buf_thr + sample_subset_n, genot + mi * sample_subset_n);
//         //    printf("index %zu", mi * sample_subset_n);
//         //    // printf("access %d", genot[mi * sample_subset_n]);
//         //}
//
//         // ng
//         // for(size_t samplei = 0; samplei < sample_subset_n; samplei++) {
//         //    genot[mi * sample_subset_n + samplei] = buf_thr[samplei];
//         //}
//
//         // if(thri == 0) {
//         //     printf("thr0 for copying mi  %d\n", mi);
//         // for(size_t samplei = 0; samplei < sample_subset_n; samplei++) {
//         //     genot[mi * sample_subset_n + samplei] = buf_thr[samplei];
//         // }
//         // }
//
//         // should use ReadHardcalls?
//     }
//
//     return nthr;
// }

/* int pgenreader_load_snv(double *genot, const char *filename, int snv_i, int sample_n_in, int *sample_subset, int sample_subset_n, int nthr) {
    PgenReader *pread = new PgenReader();
    std::string fname(filename);

    // might have better way
    std::vector<int> sample_subset_1based_v(sample_subset_n);
    std::copy(sample_subset, sample_subset + sample_subset_n, sample_subset_1based_v.begin());
    for(int &d : sample_subset_1based_v) {
        d += 1;
    }
    // printf("sample_subset_1based %d\n", sample_subset_1based_v[0]);

    pread->Load(fname, sample_n_in, sample_subset_1based_v, nthr);

    pread->Read(genot, sample_subset_n, 0, snv_i, 1);

    return nthr;
} */

int pgenreader_get_a() {
    PgenReader *pread = new PgenReader();
    return pread->get_a();
}

// tmp
// without this, raise error: undefined symbol PgenReader::PgenReader()
// PgenReader::PgenReader() {};

/* // necessary to wrap to avoid using std::string, std::vector
PgenReader *pgenreader_create() { return new PgenReader(); }
// PgenReader *pgenreaderc_create() { return new PgenReader{}; }
void pgenreader_delete(PgenReader *pread){}
// void pgenreaderc_delete(PgenReader *pread) { delete pread; }
int pgenreader_abc(PgenReader *pread, int a) { return 2 * a; }
int pgenreader_int_abc(PgenReader *pread, int a) { return get_int_abc(); }
int pgenreader_get_a(PgenReader *pread) {
    //const char* x="abc";
    //std::string str(x);

    //const char* x = "abc";
    //std::string str;
    //str.assign(x, x + 3);


    return pread->get_a();
} */
