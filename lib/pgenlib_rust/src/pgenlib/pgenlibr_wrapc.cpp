
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

int pgenreader_load_whole(double *genot, const char *filename, int snv_n_in, int sample_n_in, int *sample_subset, int sample_subset_n, int nthr) {
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

    // printf("outside of openmp\n");
    int num_threads = 1;
    int thread_num = 0;
    int mi;
#ifdef _OPENMP
    //omp_set_num_threads(nthr);
    //omp_get_num_threads();
    //printf("#threads %d\n", omp_get_num_threads());
#pragma omp parallel for private(mi)
#endif
    for(mi = 0; mi < snv_n_in; mi++) {
#ifdef _OPENMP
        thread_num = omp_get_thread_num();
        //printf("thread# %d\n", thread_num);
#endif
        pread->Read(genot + mi * sample_subset_n, sample_subset_n, thread_num, mi, 1);
    }

    // run in parallel
    // [ref](https://github.com/rgcgithub/regenie/blob/master/src/Geno.cpp#L3748)
    // pread->Read(genot, sample_n, 0, 0, 1);
    //  pread->ReadHardcalls(genot, sample_n, 0, 1, 1);
    return nthr;
}

int pgenreader_load_snvs(double *genot, const char *filename, int snv_start, int snv_end, int sample_n_in, int *sample_subset, int sample_subset_n,
                         int nthr) {
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


    /*     // no error
        int thread_num = 0;
        int mi;
        for(mi = 0; mi < snv_end - snv_start; mi++) {
            pread->Read(genot + mi * sample_subset_n, sample_subset_n, thread_num, mi + snv_start, 1);
        } */

    // still error even for nthread=1
    // int num_threads = 1;
    int mi;
    // TODO: create Vec<Vec<f64>> for buf
//#ifdef _OPENMP
//    omp_set_dynamic(0);
//    printf("cpp nthr %d\n", nthr);
//    omp_set_num_threads(nthr);
//    printf("omp #threads %d\n", omp_get_num_threads());
//#pragma omp parallel for private(mi)
//#endif
    for(mi = 0; mi < snv_end - snv_start; mi++) {
        int thread_num = 0;
//#ifdef _OPENMP
//        thread_num = omp_get_thread_num();
//        if(mi % 10000 == 0) {
//            printf("omp thread# %d\n", thread_num);
//        }
//#endif
        pread->Read(genot + mi * sample_subset_n, sample_subset_n, thread_num, mi + snv_start, 1);
        // should use ReadHardcalls?
        // load buf to genot
    }

    // this is wrong way
    /*     int mi;
        #ifdef _OPENMP
        printf("in openmp\n");
        #pragma omp parallel for private(mi)
        #endif
        for(mi = 0; mi < snv_end - snv_start; mi++) {
            pread->Read(genot + mi * sample_subset_n, sample_subset_n, 0, mi + snv_start, 1);
        } */

    return nthr;
}

int pgenreader_load_snv(double *genot, const char *filename, int snv_i, int sample_n_in, int *sample_subset, int sample_subset_n, int nthr) {
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
}

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
