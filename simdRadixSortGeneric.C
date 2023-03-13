// ===========================================================================
//
// simdRadixSortGeneric.C --
// test of SIMD implementation of generic bitwise MSB radix sort
//
// This source code file is part of the following software:
//
//    - the low-level C++ template SIMD library
//    - the SIMD implementation of the MinWarping and the 2D-Warping methods
//      for local visual homing.
//
// The software is provided based on the accompanying license agreement
// in the file LICENSE or LICENSE.doc. The software is provided "as is"
// without any warranty by the licensor and without any liability of the
// licensor, and the software may not be distributed by the licensee; see
// the license agreement for details.
//
// (C) Ralf Möller
//     Computer Engineering
//     Faculty of Technology
//     Bielefeld University
//     www.ti.uni-bielefeld.de
//
// ===========================================================================

#include "SIMDAlloc.H"
#include "SIMDRadixSortGeneric.H"
#include "SIMDRadixSortGenericThreads.H"
#include "TimeMeasurement.H"

#include <algorithm> // std::sort
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <type_traits>
#include <vector> // std::vector

// determine if parallel version of std::sort is available
// g++ 9 is supposed to support it, but not clang++ at the moment
// check whether this has changed meanwhile:
// https://en.cppreference.com/w/cpp/compiler_support
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) &&  \
  (__GNUC__ >= 9)
// 30. Sep 22 (rm): TODO: commented out, weird linker errors on g++ 10, maybe
// tbb-lib is missing
// #define HAS_PARALLEL_STD_SORT
#endif

#ifdef HAS_PARALLEL_STD_SORT
#include <execution>
#endif

// thread-based version produces and prints statistics on thread usage
// #define THREAD_STATS

// for TimeMeasurement.H
using namespace simd;
using namespace radix;

// =========================================================================
// random number generation
// =========================================================================

// -------------------------------------------------------------------------
// single element random numbers, wide uniform distribution
// -------------------------------------------------------------------------

// random bytes
template <int NUMBYTES>
void randBytes(uint8_t b[NUMBYTES])
{
  for (int i = 0; i < NUMBYTES; i++) b[i] = rand() & 0xff;
}

// distinguish between integer and floating point
template <typename T, bool IsFloatingPoint>
struct _RandWideUniform;

// integer: no restrictions
template <typename T>
struct _RandWideUniform<T, false>
{
  T generate()
  {
    T res;
    uint8_t b[sizeof(T)];
    randBytes<sizeof(T)>(b);
    memcpy((void *) &res, (void *) b, sizeof(T));
    return res;
  }
};

// floating point: we only allow finite numbers (excludes infinite, NaN)
template <typename T>
struct _RandWideUniform<T, true>
{
  T generate()
  {
    T res;
    uint8_t b[sizeof(T)];
    do {
      randBytes<sizeof(T)>(b);
      memcpy((void *) &res, (void *) b, sizeof(T));
    } while (!std::isfinite(res));
    return res;
  }
};

// hub: add constructor
template <typename T>
struct RandWideUniform : _RandWideUniform<T, std::is_floating_point<T>::value>
{
  RandWideUniform(unsigned int seed) { srand(seed); }
};

// -------------------------------------------------------------------------
// normal distribution generators
// -------------------------------------------------------------------------

const double normalMean   = 1000.0;
const double normalStdDev = 100.0;

// alternatively, we could use enable_if_t to specialize for groups of
// types

template <typename T, bool IsFloatingPoint>
struct _RandNormal;

// floating point types
template <typename T>
struct _RandNormal<T, true>
{
  std::mt19937 gen;
  std::normal_distribution<T> normal {T(normalMean), T(normalStdDev)};
  T generate() { return normal(gen); }
};

// integer types
template <typename T>
struct _RandNormal<T, false>
{
  std::mt19937 gen;
  std::normal_distribution<double> normal {normalMean, normalStdDev};
  T generate()
  {
    double minv = double(std::numeric_limits<T>::min());
    double maxv = double(std::numeric_limits<T>::max());
    double v    = std::round(normal(gen));
    return T(std::max(minv, std::min(maxv, v)));
  }
};

// hub
template <typename T>
struct RandNormal : _RandNormal<T, std::is_floating_point<T>::value>
{
  RandNormal(unsigned int seed) { this->gen.seed(seed); }
};

// -------------------------------------------------------------------------
// generate data
// -------------------------------------------------------------------------

template <typename KEYTYPE, bool WithPayload>
struct _PayloadSortIndex;

template <typename KEYTYPE>
struct _PayloadSortIndex<KEYTYPE, false>
{
  template <typename ELEMENTTYPE>
  static INLINE void set(ELEMENTTYPE &, SortIndex)
  {}
};

template <typename KEYTYPE>
struct _PayloadSortIndex<KEYTYPE, true>
{
  template <typename ELEMENTTYPE>
  static INLINE void set(ELEMENTTYPE &e, SortIndex i)
  {
    typedef typename KeyPayloadInfo<KEYTYPE, true>::UIntPayloadType PayloadType;
    PayloadType p = PayloadType(i);
    setPayload<KEYTYPE>(e, p);
  }
};

// hub
template <typename KEYTYPE, bool WithPayload>
struct PayloadSortIndex : _PayloadSortIndex<KEYTYPE, WithPayload>
{};

// noDuplicates: avoids identical keys; is slow, use only for small num
template <bool WithPayload, typename KEYTYPE,
          template <typename> class GENERATOR>
typename KeyPayloadInfo<KEYTYPE, WithPayload>::UIntElementType *generateData(
  int repeats, SortIndex num, bool noDuplicates, GENERATOR<KEYTYPE> &generator)
{
  typedef
    typename KeyPayloadInfo<KEYTYPE, WithPayload>::UIntElementType ElemType;
  // allocate contiguous data for multiple repeats
  ElemType *d =
    (ElemType *) simd_aligned_malloc(64, repeats * num * sizeof(ElemType));
  if (d == NULL) {
    fprintf(stderr, "failed to allocate memory (%s)\n", strerror(errno));
    exit(-1);
  }
  SortIndex i = 0, j;
  bool dup;
  while (i < num) {
    setKey(generator.generate(), d[i]);
    dup = false;
    if (noDuplicates) {
      for (j = 0; j < i; j++)
        if (getKey<KEYTYPE>(d[i]) == getKey<KEYTYPE>(d[j])) {
          dup = true;
          break;
        }
    }
    if (!dup) {
      PayloadSortIndex<KEYTYPE, WithPayload>::set(d[i], i);
      i++;
    }
  }
  // we have initialized first of repeats, now duplicate data to rest
  // of repeats
  for (int r = 1; r < repeats; r++)
    memcpy((void *) (d + r * num), (void *) d, num * sizeof(ElemType));
  return d;
}

template <bool WithPayload, typename KEYTYPE>
typename KeyPayloadInfo<KEYTYPE, WithPayload>::UIntElementType *generateData(
  int rndMode, unsigned int seed, int repeats, SortIndex num, bool noDuplicates)
{
  RandWideUniform<KEYTYPE> randWideUniform(seed);
  RandNormal<KEYTYPE> randNormal(seed);
  switch (rndMode) {
  case 0:
    return generateData<WithPayload>(repeats, num, noDuplicates,
                                     randWideUniform);
  case 1:
    return generateData<WithPayload>(repeats, num, noDuplicates, randNormal);
  default: fprintf(stderr, "invalid rndMode %d\n", rndMode); exit(-1);
  }
}

// =========================================================================
// check if keys are sorted
// =========================================================================

template <typename KEYTYPE, int UP, typename ELEMENTTYPE>
bool keysAreSorted(ELEMENTTYPE *d, SortIndex num)
{
  for (SortIndex i = 1; i < num; i++)
    // compareKeys has to work for std::sort, where the comparison
    // function has to return false for equal values, see
    // https://stackoverflow.com/questions/45929474, so instead of
    // !compareKeys<...,UP> we have to use compareKeys<...,1-UP> here
    // if we want to reuse compareKeys
    if (compareKeys<KEYTYPE, 1 - UP>(d[i - 1], d[i])) { return false; }
  return true;
}

// =========================================================================
// check if all payloads are present (overwrites keys!)
// =========================================================================

template <typename KEYTYPE, bool WithPayload>
struct CheckPayloads;

// without payloads
template <typename KEYTYPE>
struct CheckPayloads<KEYTYPE, false>
{
  static bool payloadsAreOk(
    typename KeyPayloadInfo<KEYTYPE, false>::UIntElementType *, SortIndex)
  {
    return true;
  }
};

// with payloads
template <typename KEYTYPE>
struct CheckPayloads<KEYTYPE, true>
{
  // NOTE: this destroys the keys!!!
  static bool payloadsAreOk(
    typename KeyPayloadInfo<KEYTYPE, true>::UIntElementType *d, SortIndex num)
  {
    // this assumes that UIntKeyType and UIntPayloadType are the same
    static_assert(
      sizeof(typename KeyPayloadInfo<KEYTYPE, true>::UIntPayloadType) ==
        sizeof(typename KeyPayloadInfo<KEYTYPE, true>::UIntKeyType),
      "CheckPayloads: size of key and payload have to be the same");
    typedef
      typename KeyPayloadInfo<KEYTYPE, true>::UIntPayloadType KeyAndPayloadType;
    KeyAndPayloadType uIntKey, uIntPayload;
    KeyAndPayloadType invalid = std::numeric_limits<KeyAndPayloadType>::max();
    if (KeyAndPayloadType(num) >= invalid) {
      fprintf(stderr, "num too large for correct payload check");
      exit(-1);
    }
    // overwrite all keys with an "invalid" value
    for (SortIndex i = 0; i < num; i++)
      memcpy((void *) (d + i), (void *) &invalid, sizeof(KeyAndPayloadType));
    // transfer d[i].payload to d.key[d[i].payload]
    for (SortIndex i = 0; i < num; i++) {
      getPayload<KEYTYPE>(d[i], uIntPayload);
      // payload could be invalid, check
      if (uIntPayload > KeyAndPayloadType(num)) return false;
      memcpy((void *) (d + uIntPayload), (void *) &uIntPayload,
             sizeof(KeyAndPayloadType));
    }
    // check whether all payloads are there
    for (SortIndex i = 0; i < num; i++) {
      memcpy((void *) &uIntKey, (void *) (d + i), sizeof(KeyAndPayloadType));
      if (SortIndex(uIntKey) != i) return false;
    }
    return true;
  }
};

// =========================================================================
// configurations
// =========================================================================

// select one from table below
// or use e.g. in tcsh
// % setenv EXTRA_USER_DEFINES "-DRADIX_CONFIG=2"
// and use EXTRA_USER_DEFINES for compiler otions in Makefile
#ifndef RADIX_CONFIG
#define RADIX_CONFIG 0
#endif

template <typename KEYTYPE, bool WITHPAYLOAD>
struct _Config
{
  typedef KEYTYPE KeyType;
  static constexpr bool WithPayload = WITHPAYLOAD;
};

template <int>
struct Config;

// 10. Feb 21 (rm): changed order

// ----- float -----
template <>
struct Config<0> : _Config<float, false>
{};
template <>
struct Config<1> : _Config<float, true>
{};
// ----- double -----
template <>
struct Config<2> : _Config<double, false>
{};
template <>
struct Config<3> : _Config<double, true>
{};
// ----- uint32_t -----
template <>
struct Config<4> : _Config<uint32_t, false>
{};
template <>
struct Config<5> : _Config<uint32_t, true>
{};
// ----- uint64_t -----
template <>
struct Config<6> : _Config<uint64_t, false>
{};
template <>
struct Config<7> : _Config<uint64_t, true>
{};
// ----- int32_t -----
template <>
struct Config<8> : _Config<int32_t, false>
{};
template <>
struct Config<9> : _Config<int32_t, true>
{};
// ----- int64_t -----
template <>
struct Config<10> : _Config<int64_t, false>
{};
template <>
struct Config<11> : _Config<int64_t, true>
{};

// =========================================================================
// aux
// =========================================================================

void printRadixThreadStats(RadixThreadStats *threadStats)
{
  printf("maxListSize %zu\n", threadStats->maxListSize);
  for (size_t i = 0; i < threadStats->elements.size(); i++)
    printf("%zu\t%ld\t%ld\n", i, threadStats->chunks[i],
           threadStats->elements[i]);
}

// =========================================================================
// main
// =========================================================================

int main(int argc, char *argv[])
{
  // argument processing
  if (argc != 10) {
    fprintf(stderr, "simdRadixSortGeneric "
                    "<rndMode> <seed> <rep> <num> <nodup> <meth> <up> <thresh> "
                    "<nthreads>\n");
    exit(-1);
  }
  int rndMode       = atoi(argv[1]);
  unsigned int seed = (unsigned int) atol(argv[2]);
  if (seed == 0) {
    seed = time(NULL);
    printf("random seed %u\n", seed);
  }
  static_assert(sizeof(long) == 8, "long should have 8 bytes");
  int rep          = atol(argv[3]);
  SortIndex num    = atol(argv[4]);
  int nodup        = atoi(argv[5]);
  int meth         = atoi(argv[6]);
  int up           = atoi(argv[7]);
  SortIndex thresh = atol(argv[8]);
  // thread parameters
  unsigned nthreads = atoi(argv[9]);
  if (nthreads < 1) {
    nthreads = std::thread::hardware_concurrency();
    printf("automatic nthreads = %u\n", nthreads);
  }
  // shorthands
  typedef Config<RADIX_CONFIG>::KeyType KeyType;
  constexpr bool WithPayload = Config<RADIX_CONFIG>::WithPayload;
  typedef typename KeyPayloadInfo<KeyType, WithPayload>::UIntElementType Data;
  // print config
  printf("RADIX_CONFIG: %d, WithPayload %d, sizeof: Data %zu KeyType %zu\n",
         RADIX_CONFIG, WithPayload, sizeof(Data), sizeof(KeyType));
  // sort
  const char *dir = up ? "upwards" : "downwards";
  // time measurements (Prep: preparation phase, Sort: summation phase)
  // none of the methods have a preparation phase, so we set it to zero
  double dtPrep = 0.0;
  // generate data for multiple repeats
  Data *dAll =
    generateData<WithPayload, KeyType>(rndMode, seed, rep, num, nodup);
  // save first 100 elements
  std::ofstream rndSampleFile;
  rndSampleFile.open(std::string("rndSample") + "_config" +
                     std::to_string(RADIX_CONFIG) + "_rndMode" +
                     std::to_string(rndMode) + ".dat");
  for (SortIndex i = 0; i < std::min(SortIndex(100), num); i++)
    rndSampleFile << getKey<KeyType>(dAll[i]) << "\n";
  rndSampleFile.close();
  // stats for thread version
#ifdef THREAD_STATS
  RadixThreadStats *threadStats = new RadixThreadStats(nthreads);
#else
  RadixThreadStats *threadStats = 0;
#endif
  printf("sorting, %d repetitions\n", rep);
  fflush(stdout);
  // multiple repeats
  Data *d                         = dAll;
  struct timespec t0Sort          = getTimeSpec();
  struct timespec t0SortMonotonic = getTimeSpecMonotonic();

  for (int r = 0; r < rep; r++, d += num) {
    // method numbers are irregular to be compatible with older code version

    // ======================================================================
    // non-threaded
    // ======================================================================

    if (meth == 0) {
      // ----- sequential radix sort -----
      if (up)
        seqRadixSort<KeyType, 1>(d, 0, num - 1, thresh);
      else
        seqRadixSort<KeyType, 0>(d, 0, num - 1, thresh);

    }

    else if (meth == 1) {
      // ----- sequential radix sort, experimental -----
      if (up)
        seqRadixSort2<KeyType, 1>(d, 0, num - 1, thresh);
      else
        seqRadixSort2<KeyType, 0>(d, 0, num - 1, thresh);

    }

    else if (meth == 20) {
      // ----- std::sort -----
      if (up)
        std::sort(d, d + num, compareKeys<KeyType, 1, Data>);
      else
        std::sort(d, d + num, compareKeys<KeyType, 0, Data>);

    }
#ifdef SIMD_RADIX_HAS_AVX512

    else if (meth == 42) {

      // ----- SIMD radix sort with compress instructions
      if (up)
        simdRadixSortCompress<KeyType, 1>(d, 0, num - 1, thresh);
      else
        simdRadixSortCompress<KeyType, 0>(d, 0, num - 1, thresh);

    }
#endif // SIMD_RADIX_HAS_AVX512

    else if (meth == 50) {

      // ----- baseline radix sort (no bit sorting at all)
      if (up)
        baselineRadixSort<KeyType, 1>(d, 0, num - 1, thresh);
      else
        baselineRadixSort<KeyType, 0>(d, 0, num - 1, thresh);

    }

    // ======================================================================
    // threaded
    // ======================================================================

    else if (meth == 100) {
      // ----- sequential radix sort with threads, no slaves -----
      if (up)
        seqRadixSortThreads<KeyType, 1>(
          RadixThreadConfig(nthreads, RadixThreadConfig::RADIX_FIFO_QUEUE, 0,
                            1.0),
          threadStats, d, 0, num - 1, thresh);
      else
        seqRadixSortThreads<KeyType, 0>(
          RadixThreadConfig(nthreads, RadixThreadConfig::RADIX_FIFO_QUEUE, 0,
                            1.0),
          threadStats, d, 0, num - 1, thresh);
    }

    else if (meth == 101) {
      // ----- sequential radix sort with threads, with slaves -----
      if (up)
        seqRadixSortThreads<KeyType, 1>(
          RadixThreadConfig(nthreads, RadixThreadConfig::RADIX_FIFO_QUEUE, 1,
                            1.0),
          threadStats, d, 0, num - 1, thresh);
      else
        seqRadixSortThreads<KeyType, 0>(
          RadixThreadConfig(nthreads, RadixThreadConfig::RADIX_FIFO_QUEUE, 1,
                            1.0),
          threadStats, d, 0, num - 1, thresh);
    }

    else if (meth == 102) {
      // ----- sequential radix sort with threads, with slaves, factor -----
      if (up)
        seqRadixSortThreads<KeyType, 1>(
          RadixThreadConfig(nthreads, RadixThreadConfig::RADIX_FIFO_QUEUE, 1,
                            2.0),
          threadStats, d, 0, num - 1, thresh);
      else
        seqRadixSortThreads<KeyType, 0>(
          RadixThreadConfig(nthreads, RadixThreadConfig::RADIX_FIFO_QUEUE, 1,
                            2.0),
          threadStats, d, 0, num - 1, thresh);
    }
#ifdef SIMD_RADIX_HAS_AVX512

    else if (meth == 142) {

      // ----- SIMD radix sort with compress instructions, no slaves ----
      if (up)
        simdRadixSortCompressThreads<KeyType, 1>(
          RadixThreadConfig(nthreads, RadixThreadConfig::RADIX_FIFO_QUEUE, 0,
                            1.0),
          threadStats, d, 0, num - 1, thresh);
      else
        simdRadixSortCompressThreads<KeyType, 0>(
          RadixThreadConfig(nthreads, RadixThreadConfig::RADIX_FIFO_QUEUE, 0,
                            1.0),
          threadStats, d, 0, num - 1, thresh);
    }

    else if (meth == 143) {
      // ----- SIMD radix sort with compress instructions, with slaves ----
      if (up)
        simdRadixSortCompressThreads<KeyType, 1>(
          RadixThreadConfig(nthreads, RadixThreadConfig::RADIX_FIFO_QUEUE, 1,
                            1.0),
          threadStats, d, 0, num - 1, thresh);
      else
        simdRadixSortCompressThreads<KeyType, 0>(
          RadixThreadConfig(nthreads, RadixThreadConfig::RADIX_FIFO_QUEUE, 1,
                            1.0),
          threadStats, d, 0, num - 1, thresh);
    }

    else if (meth == 144) {
      // ----- SIMD radix sort with compress instructions, with slaves ----
      if (up)
        simdRadixSortCompressThreads<KeyType, 1>(
          RadixThreadConfig(nthreads, RadixThreadConfig::RADIX_FIFO_QUEUE, 1,
                            2.0),
          threadStats, d, 0, num - 1, thresh);
      else
        simdRadixSortCompressThreads<KeyType, 0>(
          RadixThreadConfig(nthreads, RadixThreadConfig::RADIX_FIFO_QUEUE, 1,
                            2.0),
          threadStats, d, 0, num - 1, thresh);
    }

    else if (meth == 145) {
      // ----- SIMD radix sort with compress instructions, with slaves ----
      if (up)
        simdRadixSortCompressThreads<KeyType, 1>(
          RadixThreadConfig(nthreads, RadixThreadConfig::RADIX_FIFO_QUEUE, 1,
                            4.0),
          threadStats, d, 0, num - 1, thresh);
      else
        simdRadixSortCompressThreads<KeyType, 0>(
          RadixThreadConfig(nthreads, RadixThreadConfig::RADIX_FIFO_QUEUE, 1,
                            4.0),
          threadStats, d, 0, num - 1, thresh);
    }

    else if (meth == 146) {
      // ----- SIMD radix sort with compress instructions, with slaves ----
      if (up)
        simdRadixSortCompressThreads<KeyType, 1>(
          RadixThreadConfig(nthreads, RadixThreadConfig::RADIX_FIFO_QUEUE, 1,
                            8.0),
          threadStats, d, 0, num - 1, thresh);
      else
        simdRadixSortCompressThreads<KeyType, 0>(
          RadixThreadConfig(nthreads, RadixThreadConfig::RADIX_FIFO_QUEUE, 1,
                            8.0),
          threadStats, d, 0, num - 1, thresh);
    }
#endif // SIMD_RADIX_HAS_AVX512

#ifdef HAS_PARALLEL_STD_SORT
    else if (meth == 120) {

      // ----- std::sort -----
      if (up)
        std::sort(std::execution::seq, d, d + num,
                  compareKeys<KeyType, 1, Data>);
      else
        std::sort(std::execution::seq, d, d + num,
                  compareKeys<KeyType, 0, Data>);

    }

    else if (meth == 121) {
      // ----- std::sort -----
      if (up)
        std::sort(std::execution::par, d, d + num,
                  compareKeys<KeyType, 1, Data>);
      else
        std::sort(std::execution::par, d, d + num,
                  compareKeys<KeyType, 0, Data>);

    }

    else if (meth == 122) {
      // ----- std::sort -----
      if (up)
        std::sort(std::execution::par_unseq, d, d + num,
                  compareKeys<KeyType, 1, Data>);
      else
        std::sort(std::execution::par_unseq, d, d + num,
                  compareKeys<KeyType, 0, Data>);

    }

    else if (meth == 123) {
      // ----- std::sort -----
      if (up)
        std::sort(std::execution::unseq, d, d + num,
                  compareKeys<KeyType, 1, Data>);
      else
        std::sort(std::execution::unseq, d, d + num,
                  compareKeys<KeyType, 0, Data>);

    }
#endif

    else {

      fprintf(stderr, "invalid meth parameter %d\n", meth);
#ifndef SIMD_RADIX_HAS_AVX512
      fprintf(stderr, "possible reason: not compiled for AVX-512\n");
#endif // SIMD_RADIX_HAS_AVX512
#ifndef HAS_PARALLEL_STD_SORT
      fprintf(stderr, "possible reason: parallel std::sort not avaiable\n");
#endif
      exit(-1);
    }
  }
  // average time
  double dtSort = timeSpecDiffUsec(getTimeSpec(), t0Sort) / rep;
  double dtSortMonotonic =
    timeSpecDiffUsec(getTimeSpecMonotonic(), t0SortMonotonic) / rep;
  // check if sorted (only for the first repeat)
  bool sortOk = up ? keysAreSorted<KeyType, 1>(dAll, num) :
                     keysAreSorted<KeyType, 0>(dAll, num);
  // check payloads
  bool payloadOk =
    CheckPayloads<KeyType, WithPayload>::payloadsAreOk(dAll, num);
  if (!sortOk) printf("ERROR: is not sorted %s !!!\n", dir);
  if (!payloadOk) printf("ERROR: payloads error !!!\n");
  printf("RESULT: rndMode %d seed %u rep %d num %ld nodup %d "
         "meth %d up %d thresh %ld "
         "nthreads %u "
         "prep %f sort %f %f mono %f %f "
         "ok %d %d\n",
         rndMode, seed, rep, num, nodup, meth, up, thresh, nthreads, dtPrep,
         dtSort, dtSort / num, dtSortMonotonic, dtSortMonotonic / num, sortOk,
         payloadOk);
#ifdef THREAD_STATS
  printRadixThreadStats(threadStats);
#endif
  fflush(stdout);
  return 0;
}
