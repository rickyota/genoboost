// C Wrapper for C++
// https://akitsu-sanae.hatenablog.com/entry/2016/12/21/010321
// https://nachtimwald.com/2017/08/18/wrapping-c-objects-in-()

//#pragma once
//#include "pgenlibr.h"
//#include <string>

void foo();
int gcd(int a, int b);
// std::string ghi();
int ghi(const char *);

int ghi_wrap(const char *str) { return 89; }

class MyClass {
  public:
    MyClass();

    int get_abc();
    ~MyClass();

  private:
    int abc;
};

MyClass *myclass_create();
void myclass_delete(MyClass *m);
int myclass_get_abc(MyClass *m);

int pgenreader_get_a();

int pgenreader_load_whole(double *genot, const char *filename,int snv_n, int sample_n, int *sample_subset, int sample_subset_n, int nthr);
int pgenreader_load_snvs(double *genot, const char *filename,int snv_start,int snv_end, int sample_n, int *sample_subset, int sample_subset_n, int nthr);
int pgenreader_load_snv(double *genot, const char *filename,int snv_i, int sample_n, int *sample_subset, int sample_subset_n, int nthr);

/* PgenReader *pgenreader_create();
void pgenreader_delete(PgenReader *pread);
int pgenreader_abc(PgenReader *pread, int a);
int pgenreader_int_abc(PgenReader *pread, int a);
int pgenreader_get_a(PgenReader *pread); */
/* void pgenreader_loadallsample(PgenReader *pread, const char *fname); */


/*
#ifdef __cplusplus
extern "C" {
#endif

struct PgenReaderC;
// what for? : work without using _t
typedef struct PgenReaderC PgenReaderC_t;
PgenReaderC_t *pgenreaderc_create();
int pgenreaderc_abc(PgenReaderC_t *m, int a);
void pgenreaderc_delete(PgenReaderC_t *m);  */

/* typedef struct {
    PgenReader impl;
} PgenReaderC;

// cannot run these two func
PgenReaderC *pgenreaderc_create();
void pgenreaderc_delete(PgenReaderC *pread);
int pgenreaderc_abc(PgenReaderC *pread, int a); */

/* #ifdef __cplusplus
}
#endif */
