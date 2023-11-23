/*
 *
 * File obtained from pgenlibr R library:
 * https://github.com/chrchang/plink-ng/tree/master/2.0/pgenlibr
 *
 * License info obtained from DESCRIPTION file:
 * https://github.com/chrchang/plink-ng/blob/master/2.0/pgenlibr/DESCRIPTION
 * -----------------------------------------------------
    Package: pgenlibr
    Type: Package
    Title: PLINK 2 Binary (.pgen) Reader
    Version: 0.2
    Date: 2019-07-10
    Author: Christopher Chang
    Maintainer: Christopher Chang <chrchang@alumni.caltech.edu>
    Description: A thin wrapper over PLINK 2's core libraries which provides an R
    interface for reading .pgen files.  A minimal .pvar loader is also included.
    License: LGPL (>= 3)
    Imports: Rcpp (>= 1.0.1)
    LinkingTo: Rcpp
 * -----------------------------------------------------

 *  Modified by Joelle Mbatchou - June 29 2020
 *  - removed functions that were for R
 *  - split file to header (added link to several standard C++ libraries)
 *  - modified remaining functions to be fully C/C++ compatible
 *  - multithreaded reading of pgen file (04/13/2021)
 *
 * This file remains under LGPL v3 license (license is in same directory as this file)
 */

#include "pgenlibr.h"

int get_int_abc() {
    std::string a;
    return 182;
}

int PgenReader::abc() { return 34; }
std::string PgenReader::def() {
    std::string str("pgen_def");
    return str;
}

PgenReader::PgenReader()
    : _info_ptr(nullptr),
      //_allele_idx_offsetsp(nullptr),
      _nonref_flagsp(nullptr)
//_state_ptr(nullptr)
{}

PgenReader::~PgenReader() { Close(); }

int PgenReader::get_a() {
    /*     const char *s = "abcdef";
        //std::string filename(s);
        std::string filename;
        //https://www.geeksforgeeks.org/convert-char-to-string-in-cpp/
        filename.assign(s, 6);
        // cout << "get_a(): " << filename << endl;
        printf("get_a %s", filename); */

    int a = 19;
    return a;
}

// looks like not using parallel in this func
// PgenReader::Read() can be run in parallel
void PgenReader::Load(std::string filename, uint32_t cur_sample_ct, std::vector<int> sample_subset_1based, int nthr) {
    if(_info_ptr) {
        Close();
    }
    _info_ptr = static_cast<plink2::PgenFileInfo *>(malloc(sizeof(plink2::PgenFileInfo)));
    if(!_info_ptr) {
        fprintf(stderr, "Out of memory");
        exit(-1);
    }

    plink2::PreinitPgfi(_info_ptr);
    uint32_t cur_variant_ct = UINT32_MAX;
    const char *fname = filename.c_str();
    plink2::PgenHeaderCtrl header_ctrl;
    uintptr_t pgfi_alloc_cacheline_ct;
    char errstr_buf[plink2::kPglErrstrBufBlen];
    if(PgfiInitPhase1(fname, cur_variant_ct, cur_sample_ct, 0, &header_ctrl, _info_ptr, &pgfi_alloc_cacheline_ct, errstr_buf) !=
       plink2::kPglRetSuccess) {
        fprintf(stderr, "%s\n", &(errstr_buf[7]));
        exit(-1);
    }
    const uint32_t raw_variant_ct = _info_ptr->raw_variant_ct;
    if(header_ctrl & 0x30) {
        fprintf(stderr, "Storing of allele count information is not supported (only bi-allelic variants should be present).");
        exit(-1);
        // no need to zero-initialize this
        //_allele_idx_offsetsp = plink2::CreateRefcountedWptr(raw_variant_ct + 1);
        //_info_ptr->allele_idx_offsets = _allele_idx_offsetsp->p;
        // _info_ptr->max_allele_ct updated by PgfiInitPhase2() in this case
    }
    _info_ptr->max_allele_ct = 2;
    if((header_ctrl & 0xc0) == 0xc0) {
        // todo: load this in pvar, to enable consistency check.  we use a
        // (manually implemented) shared_ptr in preparation for this.
        const uintptr_t raw_variant_ctl = plink2::DivUp(raw_variant_ct, plink2::kBitsPerWord);
        // no need to zero-initialize this
        _nonref_flagsp = plink2::CreateRefcountedWptr(raw_variant_ctl + 1);
        _info_ptr->nonref_flags = _nonref_flagsp->p;
    }
    const uint32_t file_sample_ct = _info_ptr->raw_sample_ct;
    unsigned char *pgfi_alloc = nullptr;
    if(plink2::cachealigned_malloc(pgfi_alloc_cacheline_ct * plink2::kCacheline, &pgfi_alloc)) {
        fprintf(stderr, "Out of memory");
        exit(-1);
    }
    uint32_t max_vrec_width;
    uintptr_t pgr_alloc_cacheline_ct;
    if(PgfiInitPhase2(header_ctrl, 1, 0, 0, 0, raw_variant_ct, &max_vrec_width, _info_ptr, pgfi_alloc, &pgr_alloc_cacheline_ct, errstr_buf)) {
        if(pgfi_alloc && (!_info_ptr->vrtypes)) {
            plink2::aligned_free(pgfi_alloc);
        }
        fprintf(stderr, "%s\n", &(errstr_buf[7]));
        exit(-1);
    }
    if((!_allele_idx_offsetsp) && (_info_ptr->gflags & 4)) {
        // Note that it's safe to be ignorant of multiallelic variants when
        // phase and dosage info aren't present; GetAlleleCt() then always returns
        // 2 when that isn't actually true, and all ALTs are treated as if they
        // were ALT1, but otherwise everything works properly.
        fprintf(stderr, "Multiallelic variants and phase/dosage info simultaneously present; pvar required in this case");
        exit(-1);
    }

    _state_ptr.resize(nthr);
    _subset_index.resize(nthr);
    _pgv.resize(nthr);
    _subset_include_interleaved_vec.resize(nthr);
    _subset_cumulative_popcounts.resize(nthr);
    _subset_size.resize(nthr);
    _subset_include_vec.resize(nthr);

    for(int i = 0; i < nthr; i++) {
        _state_ptr[i] = static_cast<plink2::PgenReader *>(malloc(sizeof(plink2::PgenReader)));
        if(!_state_ptr[i]) {
            fprintf(stderr, "Out of memory");
            exit(-1);
        }
        plink2::PreinitPgr(_state_ptr[i]);
        plink2::PgrSetFreadBuf(nullptr, _state_ptr[i]);
    }

    const uintptr_t pgr_alloc_main_byte_ct = pgr_alloc_cacheline_ct * plink2::kCacheline;
    const uintptr_t sample_subset_byte_ct = plink2::DivUp(file_sample_ct, plink2::kBitsPerVec) * plink2::kBytesPerVec;
    const uintptr_t cumulative_popcounts_byte_ct = plink2::DivUp(file_sample_ct, plink2::kBitsPerWord * plink2::kInt32PerVec) * plink2::kBytesPerVec;
    const uintptr_t genovec_byte_ct = plink2::DivUp(file_sample_ct, plink2::kNypsPerVec) * plink2::kBytesPerVec;
    // const uintptr_t ac_byte_ct = plink2::RoundUpPow2(file_sample_ct * sizeof(plink2::AlleleCode), plink2::kBytesPerVec);
    // const uintptr_t ac2_byte_ct = plink2::RoundUpPow2(file_sample_ct * 2 * sizeof(plink2::AlleleCode), plink2::kBytesPerVec);
    uintptr_t multiallelic_hc_byte_ct = 0;
    if(_info_ptr->max_allele_ct != 2) {
        fprintf(stderr, "no multiallelic vaariants allowed");
        exit(-1);
        // multiallelic_hc_byte_ct = 2 * sample_subset_byte_ct + ac_byte_ct + ac2_byte_ct;
    }
    const uintptr_t dosage_main_byte_ct = plink2::DivUp(file_sample_ct, (2 * plink2::kInt32PerVec)) * plink2::kBytesPerVec;

    for(int i = 0; i < nthr; i++) {
        unsigned char *pgr_alloc;
        if(plink2::cachealigned_malloc(pgr_alloc_main_byte_ct + (2 * plink2::kPglNypTransposeBatch + 5) * sample_subset_byte_ct +
                                           cumulative_popcounts_byte_ct + (1 + plink2::kPglNypTransposeBatch) * genovec_byte_ct +
                                           multiallelic_hc_byte_ct + dosage_main_byte_ct + plink2::kPglBitTransposeBufbytes +
                                           4 * (plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8),
                                       &pgr_alloc)) {
            fprintf(stderr, "Out of memory");
            exit(-1);
        }
        plink2::PglErr reterr = PgrInit(fname, max_vrec_width, _info_ptr, _state_ptr[i], pgr_alloc);
        if(reterr != plink2::kPglRetSuccess) {
            if(!plink2::PgrGetFreadBuf(_state_ptr[i])) {
                plink2::aligned_free(pgr_alloc);
            }
            // Reducing nthr might help.
            sprintf(errstr_buf, "PgrInit() error %d", static_cast<int>(reterr));
            fprintf(stderr, "%s\n", errstr_buf);
            exit(-1);
        }
        unsigned char *pgr_alloc_iter = &(pgr_alloc[pgr_alloc_main_byte_ct]);
        _subset_include_vec[i] = reinterpret_cast<uintptr_t *>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
        _subset_include_interleaved_vec[i] = reinterpret_cast<uintptr_t *>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);

#ifdef USE_AVX2
        _subset_include_interleaved_vec[i][-3] = 0;
        _subset_include_interleaved_vec[i][-2] = 0;
#endif
        _subset_include_interleaved_vec[i][-1] = 0;

        _subset_cumulative_popcounts[i] = reinterpret_cast<uint32_t *>(pgr_alloc_iter);

        pgr_alloc_iter = &(pgr_alloc_iter[cumulative_popcounts_byte_ct]);

        _pgv[i] = std::make_shared<plink2::PgenVariant>();
        _pgv[i]->genovec = reinterpret_cast<uintptr_t *>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[genovec_byte_ct]);
        /*
           if (multiallelic_hc_byte_ct) {
           _pgv.patch_01_set = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
           pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
           _pgv.patch_01_vals = reinterpret_cast<plink2::AlleleCode*>(pgr_alloc_iter);
           pgr_alloc_iter = &(pgr_alloc_iter[ac_byte_ct]);
           _pgv.patch_10_set = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
           pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
           _pgv.patch_10_vals = reinterpret_cast<plink2::AlleleCode*>(pgr_alloc_iter);
           pgr_alloc_iter = &(pgr_alloc_iter[ac2_byte_ct]);
           } else {
           */
        _pgv[i]->patch_01_set = nullptr;
        _pgv[i]->patch_01_vals = nullptr;
        _pgv[i]->patch_10_set = nullptr;
        _pgv[i]->patch_10_vals = nullptr;

        _pgv[i]->phasepresent = reinterpret_cast<uintptr_t *>(pgr_alloc_iter);

        // }
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
        _pgv[i]->phaseinfo = reinterpret_cast<uintptr_t *>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
        _pgv[i]->dosage_present = reinterpret_cast<uintptr_t *>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
        _pgv[i]->dosage_main = reinterpret_cast<uint16_t *>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[dosage_main_byte_ct]);

        if(sample_subset_1based.size() > 0) {
            SetSampleSubsetInternal(sample_subset_1based, i);
        } else {
            _subset_size[i] = file_sample_ct;
        }
    }

    /*
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglBitTransposeBufbytes]);
    _multivar_vmaj_geno_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * genovec_byte_ct]);
    _multivar_vmaj_phasepresent_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * sample_subset_byte_ct]);
    _multivar_vmaj_phaseinfo_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * sample_subset_byte_ct]);
    _multivar_smaj_geno_batch_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 4]);
    _multivar_smaj_phaseinfo_batch_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8]);
    _multivar_smaj_phasepresent_batch_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
     pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8]);
     */
}

uint32_t PgenReader::GetRawSampleCt() const {
    if(!_info_ptr) {
        fprintf(stderr, "pgen is closed");
        exit(-1);
    }
    return _info_ptr->raw_sample_ct;
}

uint32_t PgenReader::GetSubsetSize() const { return _subset_size[0]; }

uint32_t PgenReader::GetVariantCt() const {
    if(!_info_ptr) {
        fprintf(stderr, "pgen is closed");
        exit(-1);
    }
    return _info_ptr->raw_variant_ct;
}

uint32_t PgenReader::GetAlleleCt(uint32_t variant_idx) const {
    if(!_info_ptr) {
        fprintf(stderr, "pgen is closed");
        exit(-1);
    }
    if(variant_idx >= _info_ptr->raw_variant_ct) {
        char errstr_buf[256];
        sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
        fprintf(stderr, "%s\n", errstr_buf);
        exit(-1);
    }
    if(!_allele_idx_offsetsp) {
        return 2;
    }
    fprintf(stderr, "Error: only bi-allelic variants are supported");
    exit(-1);
    // const uintptr_t* allele_idx_offsets = _allele_idx_offsetsp->p;
    // return allele_idx_offsets[variant_idx + 1] - allele_idx_offsets[variant_idx];
}

uint32_t PgenReader::GetMaxAlleleCt() const {
    if(!_info_ptr) {
        fprintf(stderr, "pgen is closed");
        exit(-1);
    }
    return _info_ptr->max_allele_ct;
}

bool PgenReader::HardcallPhasePresent() const {
    if(!_info_ptr) {
        fprintf(stderr, "pgen is closed");
        exit(-1);
    }
    return ((_info_ptr->gflags & plink2::kfPgenGlobalHardcallPhasePresent) != 0);
}

// added by J.Mbatchou (09/22/20) to check if dosages are present in PGEN file
bool PgenReader::DosagePresent() const {
    if(!_info_ptr) {
        fprintf(stderr, "pgen is closed");
        exit(-1);
    }
    return ((_info_ptr->gflags & plink2::kfPgenGlobalDosagePresent) != 0);
}

// static const int32_t kGenoRInt32Quads[1024] ALIGNV16 = QUAD_TABLE256(0, 1, 2, -3);

static const double kGenoRDoublePairs[32] ALIGNV16 = PAIR_TABLE16(0.0, 1.0, 2.0, -3.0);

void PgenReader::ReadHardcalls(double *buf, size_t const &n, int const &thr, int variant_idx, int allele_idx) {
    if(!_info_ptr) {
        fprintf(stderr, "pgen is closed");
        exit(-1);
    }
    if(static_cast<uint32_t>(variant_idx) >= _info_ptr->raw_variant_ct) {
        char errstr_buf[256];
        sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
        fprintf(stderr, "%s\n", errstr_buf);
        exit(-1);
    }
    if(n != _subset_size[thr]) {
        char errstr_buf[256];
        sprintf(errstr_buf, "buf has wrong length (%" PRIdPTR "; %u expected)", n, _subset_size[thr]);
        fprintf(stderr, "%s\n", errstr_buf);
        exit(-1);
    }
    plink2::PglErr reterr =
        PgrGet1(_subset_include_vec[thr], _subset_index[thr], _subset_size[thr], variant_idx, allele_idx, _state_ptr[thr], _pgv[thr]->genovec);
    if(reterr != plink2::kPglRetSuccess) {
        char errstr_buf[256];
        sprintf(errstr_buf, "PgrGet1() error %d", static_cast<int>(reterr));
        fprintf(stderr, "%s\n", errstr_buf);
        exit(-1);
    }
    plink2::GenoarrLookup16x8bx2(_pgv[thr]->genovec, kGenoRDoublePairs, _subset_size[thr], buf);
}

// this can be run in parallel
void PgenReader::Read(double *buf, size_t const &n, int const &thr, int variant_idx, int allele_idx) {
    if(!_info_ptr) {
        fprintf(stderr, "pgen is closed");
        exit(-1);
    }
    if(static_cast<uint32_t>(variant_idx) >= _info_ptr->raw_variant_ct) {
        char errstr_buf[256];
        sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
        fprintf(stderr, "%s\n", errstr_buf);
        exit(-1);
    }
    if(n != _subset_size[thr]) {
        char errstr_buf[256];
        sprintf(errstr_buf, "buf has wrong length (%" PRIdPTR "; %u expected)", n, _subset_size[thr]);
        fprintf(stderr, "%s\n", errstr_buf);
        exit(-1);
    }

    uint32_t dosage_ct;
    plink2::PglErr reterr = PgrGet1D(_subset_include_vec[thr], _subset_index[thr], _subset_size[thr], variant_idx, allele_idx, _state_ptr[thr],
                                     _pgv[thr]->genovec, _pgv[thr]->dosage_present, _pgv[thr]->dosage_main, &dosage_ct);
    if(reterr != plink2::kPglRetSuccess) {
        char errstr_buf[256];
        sprintf(errstr_buf, "PgrGet1D() error %d", static_cast<int>(reterr));
        fprintf(stderr, "%s\n", errstr_buf);
        exit(-1);
    }

    // this part raise error on omp; -> due to shared buf
    plink2::Dosage16ToDoubles(kGenoRDoublePairs, _pgv[thr]->genovec, _pgv[thr]->dosage_present, _pgv[thr]->dosage_main, _subset_size[thr], dosage_ct,
                              buf);
}

void PgenReader::Close() {
    // don't bother propagating file close errors for now
    if(_info_ptr) {
        // CondReleaseRefcountedWptr(&_allele_idx_offsetsp);
        CondReleaseRefcountedWptr(&_nonref_flagsp);
        if(_info_ptr->vrtypes) {
            plink2::aligned_free(_info_ptr->vrtypes);
        }
        plink2::PglErr reterr = plink2::kPglRetSuccess;
        plink2::CleanupPgfi(_info_ptr, &reterr);
        free(_info_ptr);
        _info_ptr = nullptr;
    }

    for(size_t i = 0; i < _state_ptr.size(); i++) {
        if(_state_ptr[i]) {
            plink2::PglErr reterr = plink2::kPglRetSuccess;
            plink2::CleanupPgr(_state_ptr[i], &reterr);
            if(PgrGetFreadBuf(_state_ptr[i])) {
                plink2::aligned_free(PgrGetFreadBuf(_state_ptr[i]));
            }
            free(_state_ptr[i]);
            _state_ptr[i] = nullptr;
        }
        _subset_size[i] = 0;
    }
}

void PgenReader::SetSampleSubsetInternal(std::vector<int> &sample_subset_1based, int const &thr) {
    const uint32_t raw_sample_ct = _info_ptr->raw_sample_ct;
    const uint32_t raw_sample_ctv = plink2::DivUp(raw_sample_ct, plink2::kBitsPerVec);
    const uint32_t raw_sample_ctaw = raw_sample_ctv * plink2::kWordsPerVec;
    uintptr_t *sample_include = _subset_include_vec[thr];
    plink2::ZeroWArr(raw_sample_ctaw, sample_include);
    const uint32_t subset_size = sample_subset_1based.size();
    if(subset_size == 0) {
        fprintf(stderr, "Empty sample_subset is not currently permitted");
        exit(-1);
    }
    uint32_t sample_uidx = sample_subset_1based[0] - 1;
    uint32_t idx = 0;
    uint32_t next_uidx;
    while(1) {
        if(sample_uidx >= raw_sample_ct) {
            char errstr_buf[256];
            sprintf(errstr_buf, "sample number out of range (%d; must be 1..%u)", static_cast<int>(sample_uidx + 1), raw_sample_ct);
            fprintf(stderr, "%s\n", errstr_buf);
            exit(-1);
        }
        plink2::SetBit(sample_uidx, sample_include);
        if(++idx == subset_size) {
            break;
        }
        next_uidx = sample_subset_1based[idx] - 1;

        // prohibit this since it implies that the caller expects genotypes to be
        // returned in a different order
        if(next_uidx <= sample_uidx) {
            fprintf(stderr, "sample_subset is not in strictly increasing order");
            exit(-1);
        }
        sample_uidx = next_uidx;
    }

    plink2::FillInterleavedMaskVec(sample_include, raw_sample_ctv, _subset_include_interleaved_vec[thr]);
    const uint32_t raw_sample_ctl = plink2::DivUp(raw_sample_ct, plink2::kBitsPerWord);
    plink2::FillCumulativePopcounts(sample_include, raw_sample_ctl, _subset_cumulative_popcounts[thr]);
    plink2::PgrSetSampleSubsetIndex(_subset_cumulative_popcounts[thr], _state_ptr[thr], &_subset_index[thr]);
    _subset_size[thr] = subset_size;
}
