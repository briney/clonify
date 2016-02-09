#include <map>
#include <unordered_map>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>

#include <fstream>
#include <iostream>
#include <sys/time.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <netdb.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/sysinfo.h>
#include <cstdlib>
#include <cfloat>
#include <cstddef>
#include <stdexcept>
#include <cstdio>
#include <map>
#include <sstream>
#include <cctype>
#include <iomanip>
#include <sstream>
#include <vector>
#include <limits>



#define fc_isnan(X) ((X)!=(X))


#ifndef NO_INCLUDE_FENV
#include <fenv.h>
#endif

#ifndef DBL_MANT_DIG
#error The constant DBL_MANT_DIG could not be defined.
#endif
#define T_FLOAT_MANT_DIG DBL_MANT_DIG

#ifndef LONG_MAX
#include <climits>
#endif
#ifndef LONG_MAX
#error The constant LONG_MAX could not be defined.
#endif
#ifndef INT_MAX
#error The constant INT_MAX could not be defined.
#endif

#ifndef INT32_MAX
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#endif

#ifndef HAVE_DIAGNOSTIC
#if __GNUC__ > 4 || (__GNUC__ == 4 && (__GNUC_MINOR__ >= 6))
#define HAVE_DIAGNOSTIC 1
#endif
#endif

#ifndef HAVE_VISIBILITY
#if __GNUC__ >= 4
#define HAVE_VISIBILITY 1
#endif
#endif

/* Since the public interface is given by the Python respectively R interface,
 * we do not want other symbols than the interface initalization routines to be
 * visible in the shared object file. The "visibility" switch is a GCC concept.
 * Hiding symbols keeps the relocation table small and decreases startup time.
 * See http://gcc.gnu.org/wiki/Visibility
 */
#if HAVE_VISIBILITY
#pragma GCC visibility push(hidden)
#endif

typedef int_fast32_t t_index;
#ifndef INT32_MAX
#define MAX_INDEX 0x7fffffffL
#else
#define MAX_INDEX INT32_MAX
#endif
#if (LONG_MAX < MAX_INDEX)
#error The integer format "t_index" must not have a greater range than "long int".
#endif
#if (INT_MAX > MAX_INDEX)
#error The integer format "int" must not have a greater range than "t_index".
#endif
typedef double t_float;

// self-destructing array pointer
template <typename type>
class auto_array_ptr{
private:
  type * ptr;
  auto_array_ptr(auto_array_ptr const &); // non construction-copyable
  auto_array_ptr& operator=(auto_array_ptr const &); // non copyable
public:
  auto_array_ptr()
    : ptr(NULL)
  { }
  template <typename index>
  auto_array_ptr(index const size)
    : ptr(new type[size])
  { }
  template <typename index, typename value>
  auto_array_ptr(index const size, value const val)
    : ptr(new type[size])
  {
    std::fill_n(ptr, size, val);
  }
  ~auto_array_ptr() {
    delete [] ptr; }
  void free() {
    delete [] ptr;
    ptr = NULL;
  }
  template <typename index>
  void init(index const size) {
    ptr = new type [size];
  }
  template <typename index, typename value>
  void init(index const size, value const val) {
    init(size);
    std::fill_n(ptr, size, val);
  }
  inline operator type *() const { return ptr; }
};

struct node {
  t_index node1, node2;
  t_float dist;

  /*
  inline bool operator< (const node a) const {
    return this->dist < a.dist;
  }
  */

  inline friend bool operator< (const node a, const node b) {
    return (a.dist < b.dist);
  }
};

class cluster_result {
private:
  auto_array_ptr<node> Z;
  t_index pos;

public:
  cluster_result(const t_index size)
    : Z(size)
    , pos(0)
  {}

  void append(const t_index node1, const t_index node2, const t_float dist) {
    Z[pos].node1 = node1;
    Z[pos].node2 = node2;
    Z[pos].dist  = dist;
    ++pos;
  }

  node * operator[] (const t_index idx) const { return Z + idx; }

  /* Define several methods to postprocess the distances. All these functions
     are monotone, so they do not change the sorted order of distances. */

  void sqrt() const {
    for (node * ZZ=Z; ZZ!=Z+pos; ++ZZ) {
      ZZ->dist = ::sqrt(ZZ->dist);
    }
  }

  void sqrt(const t_float) const { // ignore the argument
    sqrt();
  }

  void sqrtdouble(const t_float) const { // ignore the argument
    for (node * ZZ=Z; ZZ!=Z+pos; ++ZZ) {
      ZZ->dist = ::sqrt(2*ZZ->dist);
    }
  }

  #ifdef R_pow
  #define my_pow R_pow
  #else
  #define my_pow pow
  #endif

  void power(const t_float p) const {
    t_float const q = 1/p;
    for (node * ZZ=Z; ZZ!=Z+pos; ++ZZ) {
      ZZ->dist = my_pow(ZZ->dist,q);
    }
  }

  void plusone(const t_float) const { // ignore the argument
    for (node * ZZ=Z; ZZ!=Z+pos; ++ZZ) {
      ZZ->dist += 1;
    }
  }

  void divide(const t_float denom) const {
    for (node * ZZ=Z; ZZ!=Z+pos; ++ZZ) {
      ZZ->dist /= denom;
    }
  }
};

class doubly_linked_list {
  /*
    Class for a doubly linked list. Initially, the list is the integer range
    [0, size]. We provide a forward iterator and a method to delete an index
    from the list.

    Typical use: for (i=L.start; L<size; i=L.succ[I])
    or
    for (i=somevalue; L<size; i=L.succ[I])
  */
public:
  t_index start;
  auto_array_ptr<t_index> succ;

private:
  auto_array_ptr<t_index> pred;
  // Not necessarily private, we just do not need it in this instance.

public:
  doubly_linked_list(const t_index size)
    // Initialize to the given size.
    : start(0)
    , succ(size+1)
    , pred(size+1)
  {
    for (t_index i=0; i<size; ++i) {
      pred[i+1] = i;
      succ[i] = i+1;
    }
    // pred[0] is never accessed!
    //succ[size] is never accessed!
  }

  ~doubly_linked_list() {}

  void remove(const t_index idx) {
    // Remove an index from the list.
    if (idx==start) {
      start = succ[idx];
    }
    else {
      succ[pred[idx]] = succ[idx];
      pred[succ[idx]] = pred[idx];
    }
    succ[idx] = 0; // Mark as inactive
  }

  bool is_inactive(t_index idx) const {
    return (succ[idx]==0);
  }
};

// Indexing functions
// D is the upper triangular part of a symmetric (NxN)-matrix
// We require r_ < c_ !
#define D_(r_,c_) ( D[(static_cast<std::ptrdiff_t>(2*N-3-(r_))*(r_)>>1)+(c_)-1] )
// Z is an ((N-1)x4)-array
#define Z_(_r, _c) (Z[(_r)*4 + (_c)])

/*
  Lookup function for a union-find data structure.

  The function finds the root of idx by going iteratively through all
  parent elements until a root is found. An element i is a root if
  nodes[i] is zero. To make subsequent searches faster, the entry for
  idx and all its parents is updated with the root element.
 */
class union_find {
private:
  auto_array_ptr<t_index> parent;
  t_index nextparent;

public:
  union_find(const t_index size)
    : parent(size>0 ? 2*size-1 : 0, 0)
    , nextparent(size)
  { }

  t_index Find (t_index idx) const {
    if (parent[idx] != 0 ) { // a → b
      t_index p = idx;
      idx = parent[idx];
      if (parent[idx] != 0 ) { // a → b → c
        do {
          idx = parent[idx];
        } while (parent[idx] != 0);
        do {
          t_index tmp = parent[p];
          parent[p] = idx;
          p = tmp;
        } while (parent[p] != idx);
      }
    }
    return idx;
  }

  void Union (const t_index node1, const t_index node2) {
    parent[node1] = parent[node2] = nextparent++;
  }
};

class nan_error{};
#ifdef FE_INVALID
class fenv_error{};
#endif

/* Functions for the update of the dissimilarity array */

inline static void f_average( t_float * const b, const t_float a, const t_float s, const t_float t) {
  *b = s*a + t*(*b);
  #ifndef FE_INVALID
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif
  if (fc_isnan(*b)) {
    throw(nan_error());
  }
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif
  #endif
}

template <typename t_members>
static void NN_chain_core(const t_index N, t_float * const D, t_members * const members, cluster_result & Z2) {
/*
    N: integer
    D: condensed distance matrix N*(N-1)/2
    Z2: output data structure

    This is the NN-chain algorithm, described on page 86 in the following book:

    Fionn Murtagh, Multidimensional Clustering Algorithms,
    Vienna, Würzburg: Physica-Verlag, 1985.
*/
    t_index i;

    auto_array_ptr<t_index> NN_chain(N);
    t_index NN_chain_tip = 0;

    t_index idx1, idx2;

    t_float size1, size2;
    doubly_linked_list active_nodes(N);

    t_float min;


#ifdef FE_INVALID
    if (feclearexcept(FE_INVALID)) throw fenv_error();
#endif

    for (t_index j=0; j<N-1; ++j) {
	if (NN_chain_tip <= 3) {
	    NN_chain[0] = idx1 = active_nodes.start;
	    NN_chain_tip = 1;

	    idx2 = active_nodes.succ[idx1];
	    min = D_(idx1,idx2);
	    for (i=active_nodes.succ[idx2]; i<N; i=active_nodes.succ[i]) {
		if (D_(idx1,i) < min) {
		    min = D_(idx1,i);
		    idx2 = i;
		}
	    }
	}  // a: idx1   b: idx2
	else {
	    NN_chain_tip -= 3;
	    idx1 = NN_chain[NN_chain_tip-1];
	    idx2 = NN_chain[NN_chain_tip];
	    min = idx1<idx2 ? D_(idx1,idx2) : D_(idx2,idx1);
	}  // a: idx1   b: idx2

	do {
	    NN_chain[NN_chain_tip] = idx2;

	    for (i=active_nodes.start; i<idx2; i=active_nodes.succ[i]) {
		if (D_(i,idx2) < min) {
		    min = D_(i,idx2);
		    idx1 = i;
		}
	    }
	    for (i=active_nodes.succ[idx2]; i<N; i=active_nodes.succ[i]) {
		if (D_(idx2,i) < min) {
		    min = D_(idx2,i);
		    idx1 = i;
		}
	    }

	    idx2 = idx1;
	    idx1 = NN_chain[NN_chain_tip++];

	} while (idx2 != NN_chain[NN_chain_tip-2]);

	Z2.append(idx1, idx2, min);

	if (idx1>idx2) {
	    t_index tmp = idx1;
	    idx1 = idx2;
	    idx2 = tmp;
	}

	size1 = static_cast<t_float>(members[idx1]);
	size2 = static_cast<t_float>(members[idx2]);
	members[idx2] += members[idx1];

	// Remove the smaller index from the valid indices (active_nodes).
	active_nodes.remove(idx1);

	t_float s = size1/(size1+size2);
	t_float t = size2/(size1+size2);
	for (i=active_nodes.start; i<idx1; i=active_nodes.succ[i])
	    f_average(&D_(i, idx2), D_(i, idx1), s, t );
	// Update the distance matrix in the range (idx1, idx2).
	for (; i<idx2; i=active_nodes.succ[i])
	    f_average(&D_(i, idx2), D_(idx1, i), s, t );
	// Update the distance matrix in the range (idx2, N).
	for (i=active_nodes.succ[idx2]; i<N; i=active_nodes.succ[i])
	    f_average(&D_(idx2, i), D_(idx1, i), s, t );

    }
#ifdef FE_INVALID
    if (fetestexcept(FE_INVALID)) throw fenv_error();
#endif
}

#if HAVE_VISIBILITY
#pragma GCC visibility pop
#endif


#define size_(r_) ( ((r_<N) ? 1 : Z_(r_-N,3)) )

class linkage_output {
private:
  t_float * Z;

public:
  linkage_output(t_float * const Z_)
    : Z(Z_)
  {}

  void append(const t_index node1, const t_index node2, const t_float dist,
              const t_float size) {
    if (node1<node2) {
      *(Z++) = static_cast<t_float>(node1);
      *(Z++) = static_cast<t_float>(node2);
    }
    else {
      *(Z++) = static_cast<t_float>(node2);
      *(Z++) = static_cast<t_float>(node1);
    }
    *(Z++) = dist;
    *(Z++) = size;
  }
};

template <const bool sorted>
static void generate_SciPy_dendrogram(t_float * const Z, cluster_result & Z2, const t_index N) {
  // The array "nodes" is a union-find data structure for the cluster
  // identities (only needed for unsorted cluster_result input).
  union_find nodes(sorted ? 0 : N);
  if (!sorted) {
    std::stable_sort(Z2[0], Z2[N-1]);
  }

  linkage_output output(Z);
  t_index node1, node2;

  for (node const * NN=Z2[0]; NN!=Z2[N-1]; ++NN) {
    // Get two data points whose clusters are merged in step i.
    if (sorted) {
      node1 = NN->node1;
      node2 = NN->node2;
    }
    else {
      // Find the cluster identifiers for these points.
      node1 = nodes.Find(NN->node1);
      node2 = nodes.Find(NN->node2);
      // Merge the nodes in the union-find data structure by making them
      // children of a new node.
      nodes.Union(node1, node2);
    }
    output.append(node1, node2, NN->dist, size_(node1)+size_(node2));
  }
}

void linkage(const size_t N, double* matrix, t_index* members, double* Z) {
  cluster_result Z2(N - 1);
  NN_chain_core<t_index>(N, matrix, members, Z2);
  generate_SciPy_dendrogram<false>(Z, Z2, N);
}

#define CPY_MAX(_x, _y) ((_x > _y) ? (_x) : (_y))
#define CPY_MIN(_x, _y) ((_x < _y) ? (_x) : (_y))

#define NCHOOSE2(_n) ((_n)*(_n-1)/2)

#define CPY_BITS_PER_CHAR (sizeof(unsigned char) * 8)
#define CPY_FLAG_ARRAY_SIZE_BYTES(num_bits) (CPY_CEIL_DIV((num_bits), \
                                                          CPY_BITS_PER_CHAR))
#define CPY_GET_BIT(_xx, i) (((_xx)[(i) / CPY_BITS_PER_CHAR] >> \
                             ((CPY_BITS_PER_CHAR-1) - \
                              ((i) % CPY_BITS_PER_CHAR))) & 0x1)
#define CPY_SET_BIT(_xx, i) ((_xx)[(i) / CPY_BITS_PER_CHAR] |= \
                              ((0x1) << ((CPY_BITS_PER_CHAR-1) \
                                         -((i) % CPY_BITS_PER_CHAR))))
#define CPY_CLEAR_BIT(_xx, i) ((_xx)[(i) / CPY_BITS_PER_CHAR] &= \
                              ~((0x1) << ((CPY_BITS_PER_CHAR-1) \
                                         -((i) % CPY_BITS_PER_CHAR))))

#ifndef CPY_CEIL_DIV
#define CPY_CEIL_DIV(x, y) ((((double)x)/(double)y) == \
                            ((double)((x)/(y))) ? ((x)/(y)) : ((x)/(y) + 1))
#endif

#ifdef CPY_DEBUG
#define CPY_DEBUG_MSG(...) fprintf(stderr, __VA_ARGS__)
#else
#define CPY_DEBUG_MSG(...)
#endif


#define ISCLUSTER(_nd) ((_nd)->id >= n)
#define GETCLUSTER(_id) ((lists + _id - n))

/** The number of link stats (for the inconsistency computation) for each
    cluster. */

#define CPY_NIS 4

/** The column offsets for the different link stats for the inconsistency
    computation. */
#define CPY_INS_MEAN 0
#define CPY_INS_STD 1
#define CPY_INS_N 2
#define CPY_INS_INS 3

/** The number of linkage stats for each cluster. */
#define CPY_LIS 4

/** The column offsets for the different link stats for the linkage matrix. */
#define CPY_LIN_LEFT 0
#define CPY_LIN_RIGHT 1
#define CPY_LIN_DIST 2
#define CPY_LIN_CNT 3

void get_max_dist_for_each_cluster(const double *Z, double *max_dists, int n) {
  int *curNode;
  int ndid, lid, rid, k;
  unsigned char *lvisited, *rvisited;
  const double *Zrow;
  double max_dist;
  const int bff = CPY_FLAG_ARRAY_SIZE_BYTES(n);

  k = 0;
  curNode = (int*)malloc(n * sizeof(int));
  lvisited = (unsigned char*)malloc(bff);
  rvisited = (unsigned char*)malloc(bff);
  curNode[k] = (n * 2) - 2;
  memset(lvisited, 0, bff);
  memset(rvisited, 0, bff);
  while (k >= 0) {
    ndid = curNode[k];
    Zrow = Z + ((ndid-n) * CPY_LIS);
    lid = (int)Zrow[CPY_LIN_LEFT];
    rid = (int)Zrow[CPY_LIN_RIGHT];
    if (lid >= n && !CPY_GET_BIT(lvisited, ndid-n)) {
      CPY_SET_BIT(lvisited, ndid-n);
      curNode[k+1] = lid;
      k++;
      continue;
    }
    if (rid >= n && !CPY_GET_BIT(rvisited, ndid-n)) {
      CPY_SET_BIT(rvisited, ndid-n);
      curNode[k+1] = rid;
      k++;
      continue;
    }
    max_dist = Zrow[CPY_LIN_DIST];
    if (lid >= n) {
      max_dist = CPY_MAX(max_dist, max_dists[lid-n]);
    }
    if (rid >= n) {
      max_dist = CPY_MAX(max_dist, max_dists[rid-n]);
    }
    max_dists[ndid-n] = max_dist;
    CPY_DEBUG_MSG("i=%d maxdist[i]=%5.5f verif=%5.5f\n",
		  ndid-n, max_dist, max_dists[ndid-n]);
    k--;
  }
  free(curNode);
  free(lvisited);
  free(rvisited);
}

void form_flat_clusters_from_monotonic_criterion(const double *Z,
						 const double *mono_crit,
						 int *T, double cutoff, int n) {
  int *curNode;
  int ndid, lid, rid, k, ms, nc;
  unsigned char *lvisited, *rvisited;
  double max_crit;
  const double *Zrow;
  const int bff = CPY_FLAG_ARRAY_SIZE_BYTES(n);

  curNode = (int*)malloc(n * sizeof(int));
  lvisited = (unsigned char*)malloc(bff);
  rvisited = (unsigned char*)malloc(bff);

  /** number of clusters formed so far. */
  nc = 0;
  /** are we in part of a tree below the cutoff? .*/
  ms = -1;
  k = 0;
  curNode[k] = (n * 2) - 2;
  memset(lvisited, 0, bff);
  memset(rvisited, 0, bff);
  ms = -1;
  while (k >= 0) {
    ndid = curNode[k];
    Zrow = Z + ((ndid-n) * CPY_LIS);
    lid = (int)Zrow[CPY_LIN_LEFT];
    rid = (int)Zrow[CPY_LIN_RIGHT];
    max_crit = mono_crit[ndid-n];
    CPY_DEBUG_MSG("cutoff: %5.5f maxc: %5.5f nc: %d\n", cutoff, max_crit, nc);
    if (ms == -1 && max_crit <= cutoff) {
      CPY_DEBUG_MSG("leader: i=%d\n", ndid);
      ms = k;
      nc++;
    }
    if (lid >= n && !CPY_GET_BIT(lvisited, ndid-n)) {
      CPY_SET_BIT(lvisited, ndid-n);
      curNode[k+1] = lid;
      k++;
      continue;
    }
    if (rid >= n && !CPY_GET_BIT(rvisited, ndid-n)) {
      CPY_SET_BIT(rvisited, ndid-n);
      curNode[k+1] = rid;
      k++;
      continue;
    }
    if (ndid >= n) {
      if (lid < n) {
	if (ms == -1) {
	  nc++;
	  T[lid] = nc;
	}
	else {
	  T[lid] = nc;
	}
      }
      if (rid < n) {
	if (ms == -1) {
	  nc++;
	  T[rid] = nc;
	}
	else {
	  T[rid] = nc;
	}
      }
      if (ms == k) {
	ms = -1;
      }
    }
    k--;
  }

  free(curNode);
  free(lvisited);
  free(rvisited);  
}

void form_flat_clusters_from_dist(const double *Z, int *T,
				  double cutoff, int n) {
  double *max_dists = (double*)malloc(sizeof(double) * n);
  get_max_dist_for_each_cluster(Z, max_dists, n);
  //CPY_DEBUG_MSG("cupid: n=%d cutoff=%5.5f MD[0]=%5.5f MD[n-1]=%5.5f\n", n, cutoff, max_dists[0], max_dists[n-2]);
  form_flat_clusters_from_monotonic_criterion(Z, max_dists, T, cutoff, n);
  free(max_dists);
}

void linkage(const size_t N, double* matrix, double* Z);
void form_flat_clusters_from_dist(const double *Z, int *T,
				  double cutoff, int n) ;
void form_flat_clusters_maxclust_monocrit(const double *Z,
					  const double *mono_crit,
					  int *T, int n, int mc);
void get_max_dist_for_each_cluster(const double *Z, double *max_dists, int n);

using namespace std;


/*
 * Compilation options.
 */

#define SIMD // use optimized SSE implementation (SSE4 required)
#define NO_TIMING // remove ALL calls to timing functions

//#define PROFILING // enable profiling
//#define SLOW_MODE // disable megaclustering (much slower but could be more accurate)


/*
 * Problem definition.
 */
namespace cluster_param {
    const double cutoff = 0.35;

    const double mut_value = 0.35;
    const double epsilon = 0.001;
    const int len_penalty = 2;
}

/*
 * Input data limits. They must NEVER be exceeded.
 */
const int MAX_SEQUENCES = 3500000;
const int MAX_AVERAGE_MUTATIONS_PER_SEQUENCE = 32;
const int MAX_PARTITIONS = 7 + 1;
const int MAX_SEQUENCES_PER_PARTITION = MAX_SEQUENCES/2;
const int MAX_ESSENCES_PER_PARTITION = MAX_SEQUENCES_PER_PARTITION/5;
const int MAX_CLUSTERS_PER_PARTITION = 256 + MAX_SEQUENCES_PER_PARTITION/50;
const int MAX_AA_LENGTH = 64;
const int MAX_MUTATION_LOC = 4096;

/*
 * Tweakable parameters.
 */
#ifndef SLOW_MODE
const int MIN_CENTER_SIZE_10K = 3; // minimum number of sequences for a center (10K dataset)
const int MIN_CENTER_SIZE_100K = 9; // minimum number of sequences for a center (100K dataset)
#else
const int MIN_CENTER_SIZE_10K = INT_MAX;
const int MIN_CENTER_SIZE_100K = INT_MAX;
#endif
const double MIN_MEGACLUSTER_DISSIMILARITY = 0.40; // should probably be a bit above cluster_param::cutoff

const int CANONICAL_SAMPLES = 17; // number of samples to compute canonical_mutlist

const int MUTLIST_BLOCKS = 3; // keep track of 3*128 mutation locations: this is enough for the dataset (384 > 316)

/*
 * Constants.
 */
const int MUTTYPE_BITS = 4; // there are 16 = 1<<4 distinct types of mutations
const int SIMD_BITS_PER_BLOCK = 128; // size of SIMD register
const int SIMD_BYTES = SIMD_BITS_PER_BLOCK / 8;


#ifdef SIMD
#pragma GCC target ("sse4")
#endif

#ifdef SIMD
#define SIMD_LEV // SIMD-optimized Levenshtein distance
#define SIMD_MUT // SIMD-optimized mutation intersection

union ByteSimd {
    __m128i simd;
    uint8_t byte[SIMD_BYTES];
};
#endif



#define NOINLINE __attribute__((noinline))


#define FOR(i,a,b) for(int i=(a); i<(b); i++)
#define REP(i,n) FOR(i,0,n)

#define BLACK_BOX(x) asm("" :: "r" (x))

static inline void assert_invariant(bool condition) {
    if (!condition)
	throw std::logic_error("assert_invariant failed");
}
static inline void _range_check(bool condition, const char *message) {
    if (!condition)
	throw std::out_of_range((string)"range_check("+message+") failed");
}
#define range_check(condition) _range_check(condition, #condition)

map<string, double> globalStats;

class Stopwatch {
    struct timeval tv_startup;
public:
    Stopwatch() {
        reset();
    }
    void reset() {
#ifndef NO_TIMING
        gettimeofday(&tv_startup, NULL);
#endif
    }
    double time() const {
#ifndef NO_TIMING
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return (tv.tv_sec - tv_startup.tv_sec) + (tv.tv_usec - tv_startup.tv_usec) * 1e-6;
#else
        return 0;
#endif
    }
#ifndef NO_TIMING
    void print(string prefix) const {
        cerr << prefix << ": " << (time()*1000) << " ms" << endl;
    }
#else
    void print(const char *) const {}
#endif
#ifdef PROFILING
    void store(string prefix) const {
        double t = time();
        globalStats["time for " + prefix + " (ms)"] += t * 1000;
    }
#else
    void store(const char *) const {}
#endif
} globalStopwatch;


template<class T>
struct Slice {
    typedef T value_type;

    const T* _head;
    int _size;
    Slice<T>(const T* _head, int _size): _head(_head), _size(_size) {}
    const T& operator[](int i) const {
        return _head[i];
    }
    int size() const { return _size; }
    const T* begin() const {
        return _head;
    }
    const T* end() const {
        return _head + size();
    }
};

struct StringSlice : Slice<char> {
    StringSlice(const Slice<char>& slice): Slice<char>(slice._head, slice._size) {}
    StringSlice(const char* _head, int _size): Slice<char>(_head, _size) {}
    bool operator<(const StringSlice& o) const {
        return lexicographical_compare(begin(), end(), o.begin(), o.end());
    }
    bool operator==(const StringSlice& o) const {
        return size() == o.size() && !memcmp(begin(), o.begin(), size());
    }
    bool operator==(const string& o) const {
        return (size_t)size() == o.size() && !memcmp(begin(), o.c_str(), size());
    }
    string to_string() const {
        return string(_head, _size);
    }
};

template<class T, class S>
S& operator<<(S& stream, Slice<T> slice) {
    stream << "[";
    REP(i, slice.size())
        stream << (i ? ", " : "") << slice[i];
    return stream << "]";
}
template<class T, class S>
S& operator<<(S& stream, StringSlice slice) {
    return stream << slice.to_string();
}
template<class T, class S>
S& operator<<(S& stream, const vector<T>& vec) {
    stream << "[";
    REP(i, vec.size())
        stream << (i ? ", " : "") << vec[i];
    return stream << "]";
}
template<class T, class U, class S>
S& operator<<(S& stream, const pair<T,U>& xy) {
    return stream << "(" << xy.first << ", " << xy.second << ")";
}


StringSlice make_persistent(StringSlice slice) {
    static vector<string> archive;
    archive.push_back(slice.to_string());
    const string& s = archive.back();
    return StringSlice(s.c_str(), s.size());
}


template<class K, class V>
class Interning {
    unordered_map<K,V> mapping;
    vector<K> values;
    V next;
public:
    Interning(): next(V()) {}
    V intern(const K& key) {
        auto it = mapping.find(key);
        if (it != mapping.end())
            return it->second;
        V value = next++;
        K persistent_key = make_persistent(key);
        mapping[persistent_key] = value;
        values.push_back(persistent_key);
        return value;
    }
    const K& lookup(V id) const {
        return values[id];
    }
};

typedef uint16_t Mutation;
#ifdef SIMD_MUT
int popcnt128(__m128i x) {
    int lo = _popcnt64(_mm_cvtsi128_si64(x));
    int hi = _popcnt64(_mm_cvtsi128_si64(_mm_unpackhi_epi64(x, x)));
    return lo + hi;
}

vector<int16_t> mutloc_map(MAX_MUTATION_LOC, -1);
int16_t mutloc_count = 0;

string str(__m128i x) {
    union {
        __m128i simd;
        uint8_t byte[16];
    } v;
    v.simd = x;
    ostringstream oss;
    oss << "0x";
    REP(i, 16) {
        char buf[3];
        sprintf(buf, "%02x", v.byte[15 - i]);
        oss << buf;
    }
    return oss.str();
}

struct FastMutList {
    ByteSimd bitmap[MUTLIST_BLOCKS][MUTTYPE_BITS + 1];

    FastMutList() {
        REP(i, MUTLIST_BLOCKS)
            REP(k, MUTTYPE_BITS + 1)
                bitmap[i][k].simd = _mm_set1_epi8(0);
    }
    FastMutList(const FastMutList& o) {
        *this = o;
    }
    FastMutList& operator=(const FastMutList& o) {
        REP(i, MUTLIST_BLOCKS)
            REP(k, MUTTYPE_BITS + 1)
                bitmap[i][k].simd = o.bitmap[i][k].simd;
        return *this;
    }

    void set_mutations(Slice<Mutation> mutlist) {
        for (Mutation m : mutlist) {
            unsigned int bit = mutloc_map[m >> MUTTYPE_BITS];
            unsigned int i = bit / SIMD_BITS_PER_BLOCK;
            unsigned int j = bit % SIMD_BITS_PER_BLOCK;
            unsigned int data = m | 1 << MUTTYPE_BITS;
            REP(k, MUTTYPE_BITS + 1) {
                //bitmap[i][k].byte[j/8] &= ~(1 << (j%8));
                bitmap[i][k].byte[j/8] |= ((data>>k) & 1) << (j%8);
            }
        }
    }

    int count_intersection(const FastMutList& o) const {
        range_check(MUTTYPE_BITS == 4);
        int shared = 0;
        REP(i, MUTLIST_BLOCKS) {
            __m128i d0 = _mm_xor_si128(bitmap[i][0].simd, o.bitmap[i][0].simd);
            __m128i d1 = _mm_xor_si128(bitmap[i][1].simd, o.bitmap[i][1].simd);
            __m128i d01 = _mm_or_si128(d0, d1);
            __m128i d2 = _mm_xor_si128(bitmap[i][2].simd, o.bitmap[i][2].simd);
            __m128i d3 = _mm_xor_si128(bitmap[i][3].simd, o.bitmap[i][3].simd);
            __m128i d23 = _mm_or_si128(d2, d3);
            __m128i d = _mm_or_si128(d01, d23);
            __m128i enable = _mm_and_si128(bitmap[i][4].simd, o.bitmap[i][4].simd);
            __m128i common = _mm_andnot_si128(d, enable);
            shared += popcnt128(common);
        }
        return shared;
    }
};
#endif

int min_center_size = -1;

FastMutList* fast_pool;
int fast_pool_count = 0;

struct MutList {
#ifdef SIMD_MUT
    FastMutList* fast;
#endif
    Slice<Mutation> data;
    int _weight;

    MutList(): data(Slice<Mutation>(nullptr, 0)), _weight(0) {}
    MutList(Slice<Mutation> data, bool should_finalize): data(data), _weight(0) {
        if (should_finalize)
            finalize();
    }

    MutList clone() const {
        return *this;
    }

    int weight() const {
        return _weight;
    }
    bool empty() const {
        return data.size() == 0;
    }
    void finalize() {
#ifdef SIMD_MUT
        fast = &fast_pool[fast_pool_count++];
        fast->set_mutations(data);
#endif
    }
};
int NumSharedMuts(const MutList& m1, const MutList& m2) {
#ifdef SIMD_MUT
    return m1.fast->count_intersection(*m2.fast);
#else
    int n = 0;
    int p1 = 0, p2 = 0;
    while (p1 < m1.data.size() && p2 < m2.data.size()) {
        if (m1.data[p1] < m2.data[p2]) {
            p1++;
        } else if (m2.data[p2] < m1.data[p1]) {
            p2++;
        } else {
            p1++;
            n++;
        }
    }
    return n;
#endif
}

struct MutBag {
    vector<pair<Mutation,int> > count;
    vector<Mutation> quantized;

    MutBag& operator+=(const MutList& m2) {
        int initial_size = count.size();
        int p1 = 0, p2 = 0;
        while (p1 < initial_size && p2 < m2.data.size()) {
            if (count[p1].first < m2.data[p2]) {
                p1++;
            } else if (m2.data[p2] < count[p1].first) {
                count.push_back(pair<Mutation,int> { m2.data[p2], 1 });
                p2++;
            } else {
                count[p1].second++;
                p1++;
                p2++;
            }
        }
        while (p2 < m2.data.size()) {
            count.push_back(pair<Mutation,int> { m2.data[p2], 1 });
            p2++;
        }
        std::inplace_merge(count.begin(), count.begin() + initial_size, count.end());
        return *this;
    }
    MutList quantize(int threshold) {
        quantized.clear();
        for (auto kv : count) {
            if (kv.second >= threshold)
                quantized.push_back(kv.first);
        }
        return MutList(Slice<Mutation>(&quantized[0], (int)quantized.size()), true);
    }
};

struct EssenceKey {
    StringSlice junc;
    uint8_t v_gene, j_gene;
    uint16_t v_gene_all, j_gene_all;
    bool operator==(const EssenceKey& o) const {
        return v_gene_all == o.v_gene_all &&
            j_gene_all == o.j_gene_all &&
            junc == o.junc;
    }
    bool operator<(const EssenceKey& o) const {
        return v_gene_all < o.v_gene_all || (
                v_gene_all == o.v_gene_all && (
                    j_gene_all < o.j_gene_all || (
                        j_gene_all == o.j_gene_all &&
                            junc < o.junc)));
    }
    string to_string() const;
};

namespace std {
    template<>
    struct hash<StringSlice> {
        size_t operator()(const StringSlice& str) const {
            size_t x = 0;
            REP(i, str.size())
                x = x * 31415926535 + str[i];
            return x;
        }
    };
    template<>
    struct hash<EssenceKey> {
        size_t operator()(const EssenceKey& key) const {
            return hash<StringSlice>()(key.junc) + key.v_gene_all * 123546789 + key.j_gene_all * 987;
        }
    };
}

typedef uint64_t MutListHash;
MutListHash hash_mutations(Slice<Mutation> coll) {
    MutListHash h = 1;
    for (Mutation m : coll)
        h = h * 547129405631827 + m;
    return h;
}

struct Essence {
    EssenceKey key;
    MutList canonical_mutlist;
    unordered_map<MutListHash, MutList> mutlists;
    int _weight;
    string cluster_string;
    MutBag mutsum;

    Essence(EssenceKey key):
        key(key),
        canonical_mutlist(MutList(Slice<Mutation>(nullptr, 0), false)),
        _weight(0)
        {}
    Essence(EssenceKey key, MutList&& canonical_mutlist):
        key(key),
        canonical_mutlist(std::move(canonical_mutlist)),
        _weight(0)
        {}

    void push_mutlist(Slice<Mutation> mutlist) {
        MutList& m = mutlists[hash_mutations(mutlist)];
        if (m._weight)
            m._weight++;
        else {
            m.data = mutlist;
            m._weight = 1;
            m.finalize();
        }
        if (weight() > 1 && weight() <= CANONICAL_SAMPLES) {
            if (weight() == 2)
                for (auto& kv : mutlists)
                    mutsum += kv.second;
            else
                mutsum += m;
        }
        _weight++;
    }
    bool finalize_parsing() {
        int n = weight();
        if (n == 1)
            canonical_mutlist = mutlists.begin()->second.clone();
        else {
            int p = min(n, CANONICAL_SAMPLES);
            canonical_mutlist = mutsum.quantize(p/2);
        }
        return n >= min_center_size;
    }
    int weight() const {
        return _weight;
    }
};


int LevenshteinDistance(StringSlice s1, StringSlice s2) {
  int dp[s1.size() + 1][s2.size() + 1];
  for (int i = 0; i <= s1.size(); i++) {
    dp[i][0] = i;
  }
  for (int i = 0; i <= s2.size(); i++) {
    dp[0][i] = i;
  }
  for (int i = 1; i <= s1.size(); i++) {
    for (int j = 1; j <= s2.size(); j++) {
      int temp = min(dp[i - 1][j] + 1,
                     dp[i][j - 1] + 1);
      dp[i][j] = min(temp,
                     dp[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1));
    }
  }
  return dp[s1.size()][s2.size()];
}


#ifdef SIMD_LEV
const int MULTI_CAPACITY = SIMD_BYTES;
struct MultiString {
    ByteSimd data[MAX_AA_LENGTH];
    ByteSimd exact_size;
    int width = 0;
    int max_size = 0;
    int size() const { return max_size; }

    bool full() const { return width >= MULTI_CAPACITY; }

    void load(StringSlice s) {
        max_size = max(max_size, s.size());
        exact_size.byte[width] = s.size();
        REP(i, size())
            data[i].byte[width] = s[i];
        width++;
    }
};

void multi_levenshtein(StringSlice s1, const MultiString& s2, uint8_t *result) {
    const int bias = 1;
    __m128i dp[MAX_AA_LENGTH + 1];
    for (int i = 0; i <= s2.size(); i++) {
        dp[i] = _mm_set1_epi8(i + bias);
    }
    for (int i = 1; i <= s1.size(); i++) {
        const __m128i s1_char = _mm_set1_epi8(s1[i - 1]);
        __m128i acc = _mm_set1_epi8(i + bias);
        __m128i prev_acc = _mm_set1_epi8(i - 1 + bias);
        dp[0] = acc;
        for (int j = 1; j <= s2.size(); j++) {
            __m128i replace = _mm_add_epi8(prev_acc, _mm_cmpeq_epi8(s1_char, s2.data[j - 1].simd));
            prev_acc = dp[j];
            __m128i ins_del = _mm_min_epu8(prev_acc, acc);
            acc = _mm_add_epi8(_mm_min_epu8(ins_del, replace), _mm_set1_epi8(1));
            dp[j] = acc;
        }
    }
    REP(k, s2.width)
        result[k] = ((unsigned char*)&dp[s2.exact_size.byte[k]])[k] - bias;
}
#endif

int hamming_distance(const char* s1, const char* s2, int n) {
    int d = 0;
    REP(i, n)
        d += s1[i] != s2[i];
    return d;
}

int hamming_distance(const Essence& s1, const Essence& s2) {
    return hamming_distance(&s1.key.junc[0], &s2.key.junc[0], s1.key.junc.size());
}

int GetLD(const Essence& s1, const Essence& s2) {
    if (s1.key.junc.size() == s2.key.junc.size()) {
        return hamming_distance(&s1.key.junc[0], &s2.key.junc[0], s1.key.junc.size());
    } else {
        return LevenshteinDistance(s1.key.junc, s2.key.junc);
    }
}

int vCompare(const Essence& s1, const Essence& s2) {
    return 7 * (s1.key.v_gene != s2.key.v_gene) +
           1 * (s1.key.v_gene_all != s2.key.v_gene_all);
}

int jCompare(const Essence& s1, const Essence& s2) {
    return 7 * (s1.key.j_gene != s2.key.j_gene) +
           1 * (s1.key.j_gene_all != s2.key.j_gene_all);
}

double MutBonus(const Essence& s1, const Essence& s2, double ceiling) {
    int n1 = s1.weight();
    int n2 = s2.weight();
    int weight = n1 * n2;
    double s = 0;
    for (auto& kv1 : s1.mutlists) {
        double w1 = (double)kv1.second.weight();
        for (auto& kv2 : s2.mutlists) {
            int w2 = kv2.second.weight();
            s += min(ceiling, (double)NumSharedMuts(kv1.second, kv2.second)) * w1 * w2;
        }
    }
    return cluster_param::mut_value * (s / weight);
}

double merging_dissimilarity(const Essence& s1, const Essence& s2, double cutoff) {
    // Note: Requires that s1 != s2.
    double lenPenalty = fabs(0. + s1.key.junc.size() - s2.key.junc.size()) * cluster_param::len_penalty;
    double editLength = (double)min(s1.key.junc.size(), s2.key.junc.size());
    if (lenPenalty >= editLength*cutoff)
        return cutoff;
    int LD = GetLD(s1, s2);
    int vPenalty = vCompare(s1, s2);
    int jPenalty = jCompare(s1, s2);
    int basic = LD + vPenalty + jPenalty;
    double mutBonus = cluster_param::mut_value * NumSharedMuts(s1.canonical_mutlist, s2.canonical_mutlist);
    double withBonus = max(cluster_param::epsilon, basic - mutBonus);
    return (withBonus + lenPenalty) / editLength;
}
double fast_dissimilarity(const Essence& s1, const Essence& s2, int LD) {
    // Note: Requires that s1 != s2.
    int vPenalty = vCompare(s1, s2);
    int jPenalty = jCompare(s1, s2);
    int basic = LD + vPenalty + jPenalty;
    double mutBonus = cluster_param::mut_value * NumSharedMuts(s1.canonical_mutlist, s2.canonical_mutlist);
    double withBonus = max(cluster_param::epsilon, basic - mutBonus);
    double lenPenalty = fabs(0. + s1.key.junc.size() - s2.key.junc.size()) * cluster_param::len_penalty;
    double editLength = (double)min(s1.key.junc.size(), s2.key.junc.size());
    return (withBonus + lenPenalty) / editLength;
}

double full_dissimilarity(const Essence& s1, const Essence& s2, int LD) {
    // Note: Requires that s1 != s2.
    int vPenalty = vCompare(s1, s2);
    int jPenalty = jCompare(s1, s2);
    int basic = LD + vPenalty + jPenalty;
    double mutBonus = MutBonus(s1, s2, basic * (1/cluster_param::mut_value) - cluster_param::epsilon/cluster_param::mut_value);
    double withBonus = basic - mutBonus;
    double lenPenalty = fabs(0. + s1.key.junc.size() - s2.key.junc.size()) * cluster_param::len_penalty;
    double editLength = (double)min(s1.key.junc.size(), s2.key.junc.size());
    return (withBonus + lenPenalty) / editLength;
}


template<class T>
T* pool_new(vector<T>& pool, T&& value) {
    T* ptr = &*pool.end();
    pool.push_back(std::move(value));
    return ptr;
}


StringSlice quoted_slice(const char *ptr) {
    const char *end = strchr(ptr, '"');
    //range_check(end != 0);
    return Slice<char> { ptr, (int)(end-ptr) };
}

string itoa(int n) {
    ostringstream oss;
    oss << n;
    return oss.str();
}

int next_cluster = 1;

struct Megacluster {
    vector<Essence*> essences;

    NOINLINE
    void fill_base_matrix(vector<uint8_t>& base_matrix) {
        Stopwatch watch;
        int n = essences.size();
        base_matrix.reserve(n * (n - 1) / 2);
        REP(i, n) {
            int ni = essences[i]->key.junc.size();
            FOR(j, i+1, n) {
                int nj = essences[j]->key.junc.size();
                if (ni==nj) {
                    base_matrix.push_back(hamming_distance(&essences[i]->key.junc[0], &essences[j]->key.junc[0], ni));
                } else {
                    MultiString multi;
                    while (!multi.full() && j < n) {
                        multi.load(essences[j]->key.junc);
                        j++;
                    }
                    int size = base_matrix.size();
                    base_matrix.resize(size + multi.width);
                    multi_levenshtein(essences[i]->key.junc, multi, &base_matrix[size]);
                    j--;
                }
            }
        }
        watch.store("3a. base dissimilarity");
    }

    void cluster() {
        int n = essences.size();
        if (!n)
            return;

        vector<uint8_t> base_matrix;
        fill_base_matrix(base_matrix);
        vector<double> dist_matrix;
        {
            Stopwatch watch;
            dist_matrix.resize(n * (n - 1) / 2);
            const int STEP = 16;
            for (int i1 = 0; i1 < n; i1 += STEP) {
                int i2 = min(i1 + STEP, n);
                for (int j1 = i1+1; j1 < n; j1 += STEP) {
                    int j2 = min(j1 + STEP, n);
                    FOR(i, i1, i2)
                        FOR(j, j1, j2)
                            if (j > i) {
                                int pos = n*(n-1)/2 - (n-i)*(n-i-1)/2 + j-i-1;
                                dist_matrix[pos] = full_dissimilarity(*essences[i], *essences[j], base_matrix[pos]);
                            }
                }
            }
            watch.store("3b. final dissimilarity");
        }

        vector<t_index> weight(n);
        REP(i, n)
            weight[i] = essences[i]->weight();

        vector<int> flat_cluster(n);
        if (n > 1) {
            Stopwatch watch;
            double *linkage_matrix = new double[4 * (n - 1)];
            linkage(n, &dist_matrix[0], &weight[0], linkage_matrix);
            form_flat_clusters_from_dist(linkage_matrix, &flat_cluster[0],
                    cluster_param::cutoff, n);
            delete []linkage_matrix;
            watch.store("4. linkage");
        } else {
            flat_cluster[0] = 1;
        }

        // Note: flat_cluster has values from 1 to n_clusters inclusive.

        {
            Stopwatch watch;
            int n_clusters = *std::max_element(flat_cluster.begin(), flat_cluster.end());
            string name[n_clusters];
            REP(i, n_clusters)
                name[i] = itoa(next_cluster + i);
            REP(i, n)
                essences[i]->cluster_string = name[flat_cluster[i] - 1];
            next_cluster += n_clusters;
            watch.store("5. pre-emit");
        }
    }
};

struct MegaclusterCandidate {
    Essence* ess_ptr;
    double best_dist;
    Megacluster* best_mega;

    bool finished(int ds, int s0) {
        return cluster_param::len_penalty * ds >= best_dist * s0;
    }
    void finish() {
        best_mega->essences.push_back(ess_ptr);
    }
    void try_match(const pair<Essence, Megacluster*>& center) {
        try_match(center, GetLD(*ess_ptr, center.first));
    }
    void try_match(const pair<Essence, Megacluster*>& center, uint8_t LD) {
        double dist = fast_dissimilarity(*ess_ptr, center.first, LD);
        if (best_dist > dist) {
            best_dist = dist;
            best_mega = center.second;
        }
    }
};

struct Partition {
    vector<Mutation> mutation_pool;
    vector<Essence> essence_pool;
    vector<vector<Essence*> > essence_by_length;
    unordered_map<EssenceKey, Essence*> essence_map;

    vector<Essence*> center_candidates;
    vector<Megacluster*> candidate_mega;

    vector<vector<pair<Essence, Megacluster*> > > centers;
    int n_centers = 0;
    vector<Megacluster> megaclusters;

    Partition() {
        mutation_pool.reserve(MAX_SEQUENCES_PER_PARTITION * MAX_AVERAGE_MUTATIONS_PER_SEQUENCE);
        essence_pool.reserve(MAX_ESSENCES_PER_PARTITION);
        center_candidates.reserve(MAX_ESSENCES_PER_PARTITION);
        candidate_mega.reserve(MAX_CLUSTERS_PER_PARTITION);
        megaclusters.reserve(MAX_CLUSTERS_PER_PARTITION);

        essence_by_length.resize(MAX_AA_LENGTH + 1);
        centers.resize(MAX_AA_LENGTH + 1);
    }

    Essence& essence_lookup(const EssenceKey& key) {
        auto it = essence_map.find(key);
        if (it == essence_map.end()) {
            Essence* ess = pool_new(essence_pool, Essence(key));
            essence_by_length[key.junc.size()].push_back(ess);
            return *essence_map.insert(pair<EssenceKey, Essence*>(key, ess)).first->second;
        } else {
            return *it->second;
        }
    }

    void finalize_parsing() {
        Stopwatch timer;
        for (Essence& ess : essence_pool)
            if (ess.finalize_parsing())
                center_candidates.push_back(&ess);

        timer.store("1. finalize_parsing");
    }

    Megacluster* merge(Megacluster* a, Megacluster* b) {
        if (a == b)
            return a;

        // Merge b into a
        a->essences.insert(a->essences.end(), b->essences.begin(), b->essences.end());
        b->essences.clear();

        // Replace all references to b by a
        for (auto& row : centers)
            for (auto& pair : row)
                if (pair.second == b)
                    pair.second = a;
        for (auto& ptr : candidate_mega)
            if (ptr == b)
                ptr = a;

        return a;
    }

    void create_centers() {
        Stopwatch timer;
        std::sort(center_candidates.begin(), center_candidates.end(), [](Essence* i, Essence* j) { return i->weight() > j->weight(); });

        candidate_mega.reserve(center_candidates.size());
        REP(i, (int)center_candidates.size()) {
            bool good = true;
            Megacluster* mega = nullptr;
            REP(j, i) {
                double dist = merging_dissimilarity(*center_candidates[i], *center_candidates[j], MIN_MEGACLUSTER_DISSIMILARITY);
#if 0
                if (dist < MIN_CENTER_DISSIMILARITY) {
                    mega = candidate_mega[j];
                    good = false;
                    break;
                }
#endif
                if (dist < MIN_MEGACLUSTER_DISSIMILARITY) {
                    if (mega != 0 && mega != candidate_mega[j])
                        mega = merge(mega, candidate_mega[j]);
                    else
                        mega = candidate_mega[j];
                }
            }
            if (!mega) {
                megaclusters.push_back(Megacluster());
                mega = &megaclusters.back();
            }
            if (good) {
                int bucket = center_candidates[i]->key.junc.size();
                centers[bucket].push_back(pair<Essence, Megacluster*>(Essence(center_candidates[i]->key, center_candidates[i]->canonical_mutlist.clone()), mega));
                n_centers++;
            }
            candidate_mega.push_back(mega);
        }
        timer.store("1b. create centers");
    }

    void megacluster() {
        if (n_centers <= 1) {
            megaclusters.clear();
            megaclusters.push_back(Megacluster());
            Megacluster& mega = megaclusters.back();
            mega.essences.reserve(essence_pool.size());
            REP(len, MAX_AA_LENGTH + 1)
                for (Essence* ess_ptr : essence_by_length[len])
                    mega.essences.push_back(ess_ptr);
        } else {
            Stopwatch timer;
            REP(len, MAX_AA_LENGTH + 1) {
#ifndef SIMD_LEV
                for (Essence* ess_ptr : essence_by_length[len]) {
                    Essence& ess = *ess_ptr;
                    int s0 = len;
                    double best_dist = numeric_limits<double>::infinity();
                    Megacluster* best_mega = nullptr;
                    REP(ds, MAX_AA_LENGTH+1) {
                        if (cluster_param::len_penalty * ds >= best_dist * s0)
                            break;
                        if (s0+ds <= MAX_AA_LENGTH)
                            for (const pair<Essence, Megacluster*>& center : centers[s0+ds]) {
                                double dist = fast_dissimilarity(ess, center.first, GetLD(ess, center.first));
                                if (best_dist > dist) {
                                    best_dist = dist;
                                    best_mega = center.second;
                                }
                            }
                        if (ds > 0 && s0-ds >= 0 && cluster_param::len_penalty * ds < best_dist * (s0-ds))
                            for (const pair<Essence, Megacluster*>& center : centers[s0-ds]) {
                                double dist = fast_dissimilarity(ess, center.first, GetLD(ess, center.first));
                                if (best_dist > dist) {
                                    best_dist = dist;
                                    best_mega = center.second;
                                }
                            }
                    }
                    best_mega->essences.push_back(&ess);
                }
#else
                vector<MegaclusterCandidate> candidate;
                for (Essence* ess_ptr : essence_by_length[len]) {
                    MegaclusterCandidate cand = { ess_ptr, numeric_limits<double>::infinity(), nullptr };
                    for (const pair<Essence, Megacluster*>& center : centers[len])
                        cand.try_match(center, hamming_distance(*cand.ess_ptr, center.first));

                    if (cand.finished(1, len))
                        cand.finish();
                    else
                        candidate.push_back(cand);
                }
                vector<MegaclusterCandidate> next;
                FOR(ds, 1, MAX_AA_LENGTH+1) {
                    for (int i = len-ds; i <= len+ds; i += 2*ds) {
                        if (i < 0 || i > MAX_AA_LENGTH)
                            continue;
                        for (const auto& center : centers[i]) {
                            int j = 0;
                            while (j < (int)candidate.size()) {
                                MultiString multi;
                                while (!multi.full() && j+multi.width < (int)candidate.size())
                                    multi.load(candidate[j+multi.width].ess_ptr->key.junc);

                                uint8_t LD[MULTI_CAPACITY];
                                multi_levenshtein(center.first.key.junc, multi, LD);

                                REP(k, multi.width)
                                    candidate[j+k].try_match(center, LD[k]);
                                j += multi.width;
                            }
                        }
                    }
                    next.clear();
                    for (MegaclusterCandidate& cand : candidate)
                        if (cand.finished(ds+1, len))
                            cand.finish();
                        else
                            next.push_back(cand);
                    if (next.empty())
                        break;
                    candidate = next;
                }
                for (MegaclusterCandidate& cand : candidate)
                    cand.finish();
#endif
            }
            timer.store("2. megacluster");
        }
    }

    void cluster() {
        for (Megacluster& mega : megaclusters)
            mega.cluster();
    }

    void process() {
        finalize_parsing();
        create_centers();
        megacluster();
        cluster();
        int links = 0;
        for (Megacluster& mega : megaclusters)
            links += mega.essences.size() * (mega.essences.size() - 1) / 2;
    }
};

struct Dataset {
    vector<Partition> partition;
    vector<Essence*> seq_essence;
    Interning<StringSlice, uint8_t> v_gene_map, j_gene_map;

    void initialize() {
        partition.resize(MAX_PARTITIONS);
        seq_essence.reserve(MAX_SEQUENCES);
        next_cluster = 1;
        fast_pool = new FastMutList[MAX_SEQUENCES + MAX_SEQUENCES/5];
        fast_pool_count = 0;
    }
    void process() {
        if (min_center_size < 0) {
            min_center_size =
                seq_essence.size() >= 30000 ? MIN_CENTER_SIZE_100K :
                                              MIN_CENTER_SIZE_10K;
        }

        for (Partition& p : partition)
            p.process();
    }
} dataset;

string EssenceKey::to_string() const {
    return junc.to_string() + "/" + dataset.v_gene_map.lookup(v_gene).to_string() + "/" + dataset.j_gene_map.lookup(j_gene).to_string() + "/" + itoa(v_gene_all/256) + "/" + itoa(j_gene_all/256);
}

int digit_to_int(char c) {
    return c - '0';
}
int quoted_fast_atoi(const char *str) {
    // Assumes that the terminating character is a quote.
    int x = digit_to_int(*str++);
    while (digit_to_int(*str) >= 0)
        x = x * 10 + digit_to_int(*str++);
    return x;
}

int parse_mut_code(const char *mut) {
    uint8_t m1 = (uint8_t)mut[0] / 2 % 4;
    uint8_t m2 = (uint8_t)mut[2] / 2 % 4;
    m1 -= m2;
    if (mut[0]=='-')
        m1 = 0;
    int mut_code = m1%4 * 4 + m2;
    return mut_code;
}

size_t read_lines(vector<string>& lines, istream& input, int count) {
    string line;
    while (lines.size() < count && std::getline(input, line).good()) {
        if (line.size() >= 1 && line[line.size()-1] == '\r')
            line.resize(line.size()-1);
        lines.push_back(line);
    }
    return lines.size();
}

void throw_parse_error(const char* message) {
    cerr << "parse error: " << message << endl;
    std::exit(1);
}

bool parse_single(istream& input) {
    /* JSON ordering assumptions, used by FIELD() macro */
    const int row_v_fam = 5;
    const int row_v_all = 2;
    const int row_v_gene = 3;
    const int row_j_all = 9;
    const int row_j_gene = 11;
    const int row_junc_aa = 13;
    const int row_first_mut = 17;

    /* for each mutation */
    const int rows_per_mut = 4;
    const int mut_row_loc = 0;
    const int mut_row_mut = 1;
    const int mut_column_loc = 18;
    const int mut_column_mut = 18;
    /* end for */

    const int rows_final = 2;
    /* end JSON ordering assumptions */

    // Parse fixed-size header.
    vector<string> rows;
    if (read_lines(rows, input, row_first_mut) < row_first_mut) {
        if (rows.size() != 1)
            throw_parse_error("unexpected end of file");
        return false;
    }
#define LEN(s) (sizeof(#s) - 1)
#define FIELD(name) (rows[(row_##name)].c_str() + (LEN(name)+9))
    int v_fam = quoted_fast_atoi(FIELD(v_fam));
    Partition& partition = dataset.partition[v_fam];
    int v_all = quoted_fast_atoi(FIELD(v_all));
    uint8_t v_gene = dataset.v_gene_map.intern(quoted_slice(FIELD(v_gene)));
    int j_all = quoted_fast_atoi(FIELD(j_all));
    uint8_t j_gene = dataset.j_gene_map.intern(quoted_slice(FIELD(j_gene)));
    StringSlice junc = make_persistent(quoted_slice(FIELD(junc_aa)));
    if (junc.size() > MAX_AA_LENGTH)
        junc._size = MAX_AA_LENGTH;
    uint16_t v_gene_all = v_gene + 256 * v_all;
    uint16_t j_gene_all = j_gene + 256 * j_all;
    EssenceKey key = { junc, v_gene, j_gene, v_gene_all, j_gene_all };
    Essence& essence = partition.essence_lookup(key);
    if (dataset.seq_essence.size() >= MAX_SEQUENCES)
        throw_parse_error("MAX_SEQUENCES must be increased to process this dataset");
    dataset.seq_essence.push_back(&essence);

    // Parse mutation list.
    rows.clear();
    if (!read_lines(rows, input, 1))
        throw_parse_error("unexpected end of file");
    Slice<Mutation> mutlist(&*partition.mutation_pool.end(), 0);
    while (rows[0].size() > mut_column_loc) {
        if (read_lines(rows, input, rows_per_mut) < rows_per_mut)
            throw_parse_error("unexpected end of file");
        const int step = 3 * rows_per_mut;
        __builtin_prefetch(rows[step].c_str());

        uint16_t loc = quoted_fast_atoi(rows[mut_row_loc].c_str() + mut_column_loc);
        loc %= MAX_MUTATION_LOC;
        int mut_code = parse_mut_code(rows[mut_row_mut].c_str() + mut_column_mut);
#ifdef SIMD_LEV
        if (mutloc_map[loc] < 0) {
            mutloc_map[loc] = mutloc_count++;
            if (mutloc_count == MUTLIST_BLOCKS * SIMD_BITS_PER_BLOCK) {
                // There aren't any bits left, we must reuse some of them.
                // This shouldn't happen if MUTLIST_BLOCKS is high enough.
                mutloc_count = MUTLIST_BLOCKS * SIMD_BITS_PER_BLOCK / 2;
            }
        }
#endif
        Mutation mutation = loc << MUTTYPE_BITS | mut_code;
        partition.mutation_pool.push_back(mutation);
        mutlist._size++;
        rows.clear();
        if (!read_lines(rows, input, 1))
            throw_parse_error("unexpected end of file");
    }
    essence.push_mutlist(mutlist);

    int n_skip = mutlist.size() != 0 ? rows_final + 1 : rows_final;
    if (read_lines(rows, input, n_skip) < n_skip)
        throw_parse_error("unexpected end of file");

    return true;
}

void parse(istream& input) {
    Stopwatch timer;
    vector<string> dummy_lines;
    if (!read_lines(dummy_lines, input, 1))
        throw_parse_error("unexpected end of file");
    while (parse_single(input));
    timer.store("0. parsing");
}

struct ABCSpeedup {
    ABCSpeedup() {
        dataset.initialize();
    }
    vector<string> cluster(istream& input) {
        parse(input);

        dataset.process();

        Stopwatch timer;
        vector<string> result;
        result.reserve(dataset.seq_essence.size());
        REP(i, (int)dataset.seq_essence.size())
            result.push_back(dataset.seq_essence[i]->cluster_string);
        timer.store("5b. emit");
        return result;
    }
};

int main(int argc, char **argv) {
    std::ios::sync_with_stdio(false);
    if (argc > 1 && !strncmp(argv[1], "--min-center-size=", strlen("--min-center-size="))) {
        if (!strcmp(argv[1], "--min-center-size=disabled"))
            min_center_size = INT_MAX;
        else
            min_center_size = atoi(argv[1] + strlen("--min-center-size="));
        if (min_center_size < 1) {
            cerr << "invalid parameter for min-center-size" << endl;
            std::exit(1);
        }
        argc--, argv++;
    }
    if (argc != 3) {
        cerr << "usage: " << argv[0] << " [options...] input.json output.txt" << endl;
        cerr << endl;
        cerr << "Options:" << endl;
        cerr << "  --min-center-size=INTEGER:  set the minimum center size to INTEGER." << endl;
        cerr << "                              Controls the speed-vs-accuracy tradeoff." << endl;
        cerr << "                              Reasonable values range from 10 to N/10000." << endl;
        cerr << "                              (e.g. 10 to 100 when processing 1M antibodies)" << endl;
        cerr << "                              Lower values are generally more accurate." << endl;
        cerr << "                              Higher values are generally faster, but values" << endl;
        cerr << "                              that are too high run extremely slowly." << endl;
        cerr << "  --min-center-size=disabled: disable megaclustering." << endl;
        cerr << "                              Warning: very slow and memory hungry!" << endl;
        cerr << "                              Don't use with large datasets." << endl;
        std::exit(1);
    }
    string input_path = (string)argv[1];
    string output_path = (string)argv[2];

    ABCSpeedup abc;
    ifstream input(input_path);
    if (!input.good()) {
        cerr << "couldn't open JSON file '" << input_path << "'" << endl;
        std::exit(1);
    }

    struct timeval tv_startup, tv;
    gettimeofday(&tv_startup, NULL);
    vector<string> ret = abc.cluster(input);
    gettimeofday(&tv, NULL);
    double t = ((tv.tv_sec - tv_startup.tv_sec) + (tv.tv_usec - tv_startup.tv_usec) * 1e-6) * 1000;
#ifndef QUIET
    cout << "parsing and clustering time: " << (int)t << " ms" << endl;
#endif

    {
        ofstream of(output_path);
        for (string line : ret)
            of << line << endl;
    }
}
