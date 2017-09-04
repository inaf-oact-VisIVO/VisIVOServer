#ifndef PLANCK_SSE_UTILS_H
#define PLANCK_SSE_UTILS_H

#if (defined(__SSE__))

#include <xmmintrin.h>

#define PLANCK_HAVE_SSE

#ifdef __cplusplus
extern "C" {
#endif

typedef __m128 v4sf; /* vector of 4 floats (SSE1) */

typedef union {
  float f[4];
  v4sf v;
} V4SF;

static inline v4sf build_v4sf (float a, float b, float c, float d)
  { return _mm_set_ps(d,c,b,a); }
static inline void read_v4sf (v4sf v, float *a, float *b, float *c, float *d)
  {
  V4SF tmp;
  tmp.v = v;
  if (a) *a=tmp.f[0];
  if (b) *b=tmp.f[1];
  if (c) *c=tmp.f[2];
  if (d) *d=tmp.f[3];
  }

#ifdef __cplusplus
}
#endif

#endif

#if (defined(__SSE2__))

#include <emmintrin.h>

#define PLANCK_HAVE_SSE2

#ifdef __cplusplus
extern "C" {
#endif

typedef __m128d v2df; /* vector of 2 doubles (SSE2) */

typedef union {
  double d[2];
  unsigned long long ull[2];
  v2df v;
} V2DF;

typedef struct {
  v2df a,b;
} v2df2;
typedef struct {
  V2DF a,b;
} V2DF2;

static inline v2df build_v2df (double a, double b)
  { return _mm_set_pd(b,a); }
static inline v2df build_v2ull (unsigned long long a, unsigned long long b)
  {
  V2DF tmp;
  tmp.ull[0]=a; tmp.ull[1]=b;
  return tmp.v;
  }
static inline void read_v2df (v2df v, double *a, double *b)
  { _mm_store_sd(a,v); _mm_storeh_pd(b,v); }

static inline double minabs_v2df(v2df val, v2df inv_signmask)
  {
  V2DF v2;
  v2.v=val;
  v2.v=_mm_and_pd(v2.v,inv_signmask); /* absolute value */
  v2.v=_mm_min_pd(v2.v,_mm_shuffle_pd (v2.v,v2.v,_MM_SHUFFLE2(0,1)));
  return v2.d[0];
  }

static inline double maxabs_v2df(v2df val, v2df inv_signmask)
  {
  V2DF v2;
  v2.v=val;
  v2.v=_mm_and_pd(v2.v,inv_signmask); /* absolute value */
  v2.v=_mm_max_pd(v2.v,_mm_shuffle_pd (v2.v,v2.v,_MM_SHUFFLE2(0,1)));
  return v2.d[0];
  }

#ifdef __cplusplus
}
#endif

#endif

#endif
