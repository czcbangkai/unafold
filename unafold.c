#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>



// util.h
#define TURN 3

#ifdef NO_GU_BASEPAIRS
const int BPI[6][6] = {{6, 6, 6, 0, 6, 6},
		       {6, 6, 1, 6, 6, 6},
		       {6, 2, 6, 6, 6, 6},
		       {3, 6, 6, 6, 6, 6},
		       {6, 6, 6, 6, 6, 6},
		       {6, 6, 6, 6, 6, 6}};
#else
const int BPI[6][6] = {{6, 6, 6, 0, 6, 6},
		       {6, 6, 1, 6, 6, 6},
		       {6, 2, 6, 4, 6, 6},
		       {3, 6, 5, 6, 6, 6},
		       {6, 6, 6, 6, 6, 6},
		       {6, 6, 6, 6, 6, 6}};
#endif
#define basePairIndex(a, b) BPI[a][b]

int min3(int a, int b, int c)
{
  if (a <= b && a <= c)
    return a;
  if (b <= c)
    return b;
  return c;
}



// energy.h
#ifndef ENERGY
# define ENERGY double
#endif

#ifndef PRECISION
# define PRECISION 1
#endif

#ifdef INTEGER
# define isFinite(x) (x < INFINITY / 2)
#else
# define isFinite(x) (x < DBL_MAX)
#endif

#ifdef INFINITY
# undef INFINITY
#endif

#ifdef INTEGER
//#define scale(d) (isinf(d) ? INFINITY : floor((d) * PRECISION + 0.5))
const ENERGY INFINITY = 999999;
#else
//#define scale(d) ((d) * PRECISION)
const ENERGY INFINITY = 1.0 / 0.0;
#endif

extern const ENERGY INFINITY;
extern const double R;
extern const char BASES[5];
extern const char BASE_PAIRS[6][4];

struct triloop { char loop[5]; ENERGY energy; };
struct hexaloop { char loop[8]; ENERGY energy; };
struct tloop { char loop[6]; ENERGY energy; };

int triloopcmp(const void*, const void*);
int hexaloopcmp(const void*, const void*);
int tloopcmp(const void*, const void*);



// energy.c
int hexaloopcmp(const void* loop1, const void* loop2)
{
  int i;
  const unsigned char* h1 = loop1;
  const struct hexaloop *h2 = loop2;

  for (i = 0; i < 8; ++i)
    if (h1[i] < h2->loop[i])
      return -1;
    else if (h1[i] > h2->loop[i])
      return 1;

  return 0;
}
int tloopcmp(const void* loop1, const void* loop2)
{
  int i;
  const unsigned char* h1 = loop1;
  const struct tloop *h2 = loop2;

  for (i = 0; i < 6; ++i)
    if (h1[i] < h2->loop[i])
      return -1;
    else if (h1[i] > h2->loop[i])
      return 1;

  return 0;
}
int triloopcmp(const void* loop1, const void* loop2)
{
  int i;
  const unsigned char* h1 = loop1;
  const struct triloop *h2 = loop2;

  for (i = 0; i < 5; ++i)
    if (h1[i] < h2->loop[i])
      return -1;
    else if (h1[i] > h2->loop[i])
      return 1;

  return 0;
}



// hybrid-ss-min.c
struct stackNode
{
    int i;
    int j;
    int matrix; /* [0, 1, 2, 3, 4] ~ [Q', QM, Q, Q5, Q3] */
    struct stackNode* next;
};

struct constraintListNode
{
    int i, j, k, l;
    struct constraintListNode* next;
} *prohibitList, *forceList;
#if ENABLE_FORCE
char* g_ssok;
#define ssOK(i, j) g_ssok[(i) * (g_len + 2) + j]
#else
#define ssOK(i, j) 1
#endif

struct pairListNode
{
    int i, j, length;
    ENERGY E;
    struct pairListNode* next;
} *pairList;

#define auPenalty(i, j) g_aup[g_seq[i]][g_seq[j]]

// global variables
int g_len;
ENERGY *q, *qprime, *qm, *q5, *q3;
ENERGY *qprime_ip;
double RT;

int g_debug, g_nodangle, g_allPairs, g_maxLoop, g_simple, g_noisolate, g_prefilter1, g_prefilter2, g_mfoldMax, g_mfoldP, g_mfoldW, g_quiet, g_maxBP, g_circular;
int g_numSeqs;
char *g_name, *g_string;
unsigned char* g_seq; /* [0-4] for [A,C,G,TU,N] */
char *g_file, *g_prefix, *g_bpFile;
int g_oneTemp, g_firstSeq;

ENERGY g_dangle3[5][5][6];
ENERGY g_dangle5[5][5][6];
ENERGY g_stack[5][5][5][5];
ENERGY g_hairpinLoop[30];
ENERGY g_interiorLoop[30];
ENERGY g_bulgeLoop[30];
ENERGY g_sint2[7][7][5][5];
ENERGY g_asint1x2[7][7][5][5][5];
ENERGY g_sint4[7][7][5][5][5][5];
ENERGY g_tstackh[5][5][5][5];
ENERGY g_tstacki[5][5][5][5];
ENERGY g_tstackm[5][5][6][6];
ENERGY g_tstacke[5][5][6][6];
struct triloop* g_triloop; int numTriloops;
struct tloop* g_tloop; int numTloops;
struct hexaloop* g_hexaloop; int numHexaloops;
ENERGY g_multi[3];
ENERGY g_misc[13];
ENERGY g_aup[5][5];

ENERGY Eh(int, int);
ENERGY Es(int, int);
ENERGY a(int); 
ENERGY b(int); 
ENERGY c(int); 
ENERGY INFINITY_VAL(int); 
int Eval_isFinite(ENERGY);
int Eval_ssOK(int,int);
ENERGY Ebi_sizePenalty(int); 
ENERGY Ebi_stacking(int,int); 
ENERGY Ebi_asymmetry(int,int); 
ENERGY Ebi_Bulge1(int,int,int,int); 
ENERGY Ebi_Bulge(int,int,int,int,int); 
ENERGY Ebi_iloop1x1(int,int,int,int); 
ENERGY Ebi_iloop1x2(int,int,int,int); 
ENERGY Ebi_iloop2x1(int,int,int,int); 
ENERGY Ebi_iloop2x2(int,int,int,int);

ENERGY Eh(int i, int j)
{
    ENERGY energy;
    int loopSize = j - i - 1;
    int k;
    
    if (loopSize < TURN)
        return INFINITY;
    
    if (i <= g_len && g_len < j)
        return INFINITY;
    else if (i > g_len)
    {
        i -= g_len;
        j -= g_len;
    }
    
#if ENABLE_FORCE
    if (!ssOK(i + 1, j - 1))
        return INFINITY;
#endif
    
    if (loopSize <= 30)
        energy = g_hairpinLoop[loopSize - 1];
    else
        energy = g_hairpinLoop[29] + g_misc[12] * log((double) loopSize / 30);
    
    if (loopSize > 3)
        energy += g_tstackh[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]];
    else
        energy += auPenalty(i, j);
    
    if (loopSize == 3)
    {
        struct triloop* loop;
        if (numTriloops)
            if ((loop = bsearch(g_seq + i, g_triloop, numTriloops, sizeof(struct triloop), triloopcmp)))
                energy += loop->energy;
    }
    else if (loopSize == 4)
    {
        struct tloop* loop;
        if (numTloops)
            if ((loop = bsearch(g_seq + i, g_tloop, numTloops, sizeof(struct tloop), tloopcmp)))
                energy += loop->energy;
    }
    else if (loopSize == 6)
    {
        struct hexaloop* loop;
        if (numHexaloops)
            if ((loop = bsearch(g_seq + i, g_hexaloop, numHexaloops, sizeof(struct hexaloop), hexaloopcmp)))
                energy += loop->energy;
    }
    
    /* GGG */
    if (i >= 3 && g_seq[i - 2] == 2 && g_seq[i - 1] == 2 && g_seq[i] == 2 && g_seq[j] == 3)
        energy += g_misc[8];
    
    /* poly-C */
    if (loopSize == 3 && g_seq[i + 1] == 1 && g_seq[i + 2] == 1 && g_seq[i + 3] == 1)
        energy += g_misc[11];
    else
    {
        for (k = 1; k <= loopSize; ++k)
            if (g_seq[i + k] != 1)
                return energy;
        energy += g_misc[9] * loopSize + g_misc[10];
    }
    
    return energy;
}

ENERGY Es(int i, int j)
{
    if (i >= j)
        return INFINITY;
    /* fputs("Error in Es(): i isn't less than j\n", stderr); */
    
    if (i == g_len || j == g_len + 1)
        return INFINITY;
    
    if (i > g_len)
        i -= g_len;
    if (j > g_len)
        j -= g_len;
    
    return g_stack[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]];
}

// Added by Tanveer
ENERGY End(int i, int j)
{
    return auPenalty(i, j);
}
ENERGY a(int i)
{
    return g_multi[0];
}

ENERGY b(int i)
{
    return g_multi[1];
}

ENERGY c(int i)
{
    return g_multi[2];
}

ENERGY INFINITY_VAL(int i)
{
    return INFINITY;
}

int Eval_isFinite(ENERGY E)
{
    return isFinite(E);
}

int Eval_ssOK(int i, int j)
{
    return (ssOK(i,j));
}

ENERGY Ebi_sizePenalty(int iloopSize)
{
    ENERGY loopEnergy;
    
    if(iloopSize <= 30)
    {
        loopEnergy = g_interiorLoop[iloopSize-1];
    }
    else
    {
        loopEnergy = g_interiorLoop[29] + g_misc[12] * log((double)(iloopSize)/30);
    }
    return loopEnergy;
}

ENERGY Ebi_stacking(int i, int j)
{
    ENERGY loopEnergy;
    
    loopEnergy = g_tstacki[g_seq[i]][g_seq[j]][g_seq[i+1]][g_seq[j-1]];
    return loopEnergy;
}


ENERGY Ebi_asymmetry(int loopSize1, int loopSize2)
{
    ENERGY asPenalty;
    //loopSize1 = [ip-i-1] : ker[j,jp,i+ip]
    //loopSize2 = [j-jp-1] : ker[i,ip,j+jp]
    
    //abs(ip-i-1 - (j-jp-1)) = abs(ip+jp-i-j) : ker[i-j,i+ip,j+jp]
    
    //min(ip-i-1, j-jp-1) = ker[i+ip,j+jp]
    
    asPenalty = abs(loopSize1 - loopSize2) * g_misc[min3(4, loopSize1, loopSize2) - 1];
    if (asPenalty > g_misc[4])
    {
        asPenalty = g_misc[4];
    }
    return asPenalty;
}

ENERGY Ebi_Bulge1(int i, int j, int ii, int jj)
{
    return (g_bulgeLoop[0] + g_stack[g_seq[i]][g_seq[j]][g_seq[ii]][g_seq[jj]]);
}

ENERGY Ebi_Bulge(int i, int j, int ii, int jj, int loopSize)
{
    if (loopSize <= 30)
    {
        return (g_bulgeLoop[loopSize-1] + auPenalty(i, j) + auPenalty(ii, jj));
    }
    else
    {
        return (g_bulgeLoop[29] + g_misc[12] * log((double) loopSize / 30) + auPenalty(i, j) + auPenalty(ii, jj));
    }
}

ENERGY Ebi_iloop1x1(int i, int j, int ii, int jj)
{
    return g_sint2[basePairIndex(g_seq[i], g_seq[j])][basePairIndex(g_seq[ii], g_seq[jj])][g_seq[i+1]][g_seq[j-1]];
}

ENERGY Ebi_iloop1x2(int i, int j, int ii, int jj)
{
    return g_asint1x2[basePairIndex(g_seq[i], g_seq[j])][basePairIndex(g_seq[ii], g_seq[jj])][g_seq[i+1]][g_seq[j-1]][g_seq[j-2]];
}

ENERGY Ebi_iloop2x1(int i, int j, int ii, int jj)
{
    return g_asint1x2[basePairIndex(g_seq[jj], g_seq[ii])][basePairIndex(g_seq[j], g_seq[i])][g_seq[jj + 1]][g_seq[ii - 1]][g_seq[ii - 2]];
}

ENERGY Ebi_iloop2x2(int i, int j, int ii, int jj)
{
    return g_sint4[basePairIndex(g_seq[i], g_seq[j])][basePairIndex(g_seq[ii], g_seq[jj])][g_seq[i + 1]][g_seq[j - 1]][g_seq[i + 2]][g_seq[j - 2]];
} 



// fillMatrices1_unafold.c
// Common Macros
#define MAX(x,y)  ((x)>(y) ? (x) : (y))
#define max(x,y)  ((x)>(y) ? (x) : (y))
#define MIN(x,y)  ((x)>(y) ? (y) : (x))
#define min(x,y)  ((x)>(y) ? (y) : (x))
#define CEILD(n,d) (int)ceil(((double)(n))/((double)(d)))
#define ceild(n,d) (int)ceil(((double)(n))/((double)(d)))
#define FLOORD(n,d) (int)floor(((double)(n))/((double)(d)))
#define floord(n,d) (int)floor(((double)(n))/((double)(d)))
#define CDIV(x,y)  CEILD((x),(y))
#define div(x,y)  CDIV((x),(y))
#define FDIV(x,y)  FLOORD((x),(y))
#define LB_SHIFT(b,s) ((int)ceild(b,s) * s)
#define MOD(i,j)  ( (i)>=0 ? (i)%(j) : (j-1)-((-(i))%(j)) )
// Reduction Operators
#define RADD(x,y)  ((x)+=(y))
#define RMUL(x,y)  ((x)*=(y))
#define RMAX(x,y)  ((x)=MAX((x),(y)))
#define RMIN(x,y)  ((x)=MIN((x),(y)))
 
//Local Function Declarations
int fillMatrices1_unafold_QBI_SR2_reduce_1(long,long,int,int,int,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*);
int fillMatrices1_unafold_Qprime_body1_reduce_1(long,long,int,int,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*);
int fillMatrices1_unafold_Qprime_body2_reduce_1(long,long,int,int,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*);
int fillMatrices1_unafold_Qprime_body3_reduce_1(long,long,int,int,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*);
int fillMatrices1_unafold_QM_body_reduce_1(long,long,int,int,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*);
int fillMatrices1_unafold_QBI_SR1_init_reduce_1(long,long,int,int,int,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*);
int fillMatrices1_unafold_QBI_SR1_add_reduce_1(long,long,int,int,int,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*);
//Memory Macros

#define Q(i, j) Q[(N - 1) * (i - 1) + j - 1]
#define Qprime(i, j) Qprime[(N - 1) * (i - 1) + j - 1]
#define QM(i, j) QM[(N - 1) * (i - 1) + j - 1]
#define Qprime_ip(i, j) Qprime_ip[(N - 1) * (i - 1) + j - 1]

//#define Qprime_ip(i,j) Qprime_ip[(i)*(N+1)+j]
//#define Q(i,j) Q[(i)*(N+1)+j]
#define QBI_SR1(i,j,ip) QBI_SR1[(i)*((N+1)*(N-6))+(j)*(N-6)+ip]
#define QBI_SR1_init(i,j,ip) QBI_SR1_init[(i)*((N+1)*(N-6))+(j)*(N-6)+ip]
#define QBI_SR1_add(i,j,ip) QBI_SR1_add[(i)*((N+1)*(N-8))+(j)*(N-8)+ip]
//#define Qprime(i,j) Qprime[(i)*(N+1)+j]
#define Qprime_body2() Qprime_body2[0]
#define Qprime_body3() Qprime_body3[0]
//#define QM(i,j) QM[(i)*(N+1)+j]
#define QBI_SR2(i) QBI_SR2[i]
//Function bodies
void fillMatrices1_unafold(long N,long MAXLOOP,int* Qprime_ip,int* Q,int* Qprime,int* QM){
    // Parameter checking
    if (!((N-8>= 0 && MAXLOOP-8>= 0))) {
        printf("The value of parameters are not vaild.\n");
        exit(-1);
    }
    int* QBI_SR1 = malloc(sizeof(int)*((N-9)*(N+1)*(N-6)));
    if (QBI_SR1 == NULL) {
        printf("Failed to allocate memory for QBI_SR1 : size=%ld\n", ((N-9)*(N+1)*(N-6)) * sizeof(int));
        exit(-1);
    }
    int* QBI_SR2 = malloc(sizeof(int)*(N-6));
    if (QBI_SR2 == NULL) {
        printf("Failed to allocate memory for QBI_SR2 : size=%ld\n", (N-6) * sizeof(int));
        exit(-1);
    }
    int* Qprime_body2 = malloc(sizeof(int)*(1));
    if (Qprime_body2 == NULL) {
        printf("Failed to allocate memory for Qprime_body2 : size=%ld\n", (1) * sizeof(int));
        exit(-1);
    }
    int* Qprime_body3 = malloc(sizeof(int)*(1));
    if (Qprime_body3 == NULL) {
        printf("Failed to allocate memory for Qprime_body3 : size=%ld\n", (1) * sizeof(int));
        exit(-1);
    }
    int* QBI_SR1_init = malloc(sizeof(int)*((N-9)*(N+1)*(N-6)));
    if (QBI_SR1_init == NULL) {
        printf("Failed to allocate memory for QBI_SR1_init : size=%ld\n", ((N-9)*(N+1)*(N-6)) * sizeof(int));
        exit(-1);
    }
    int* QBI_SR1_add = malloc(sizeof(int)*((N-11)*(N+1)*(N-8)));
    if (QBI_SR1_add == NULL) {
        printf("Failed to allocate memory for QBI_SR1_add : size=%ld\n", ((N-11)*(N+1)*(N-8)) * sizeof(int));
        exit(-1);
    }
    #define S0(i,j,ip) QBI_SR2(ip) = fillMatrices1_unafold_QBI_SR2_reduce_1(N,MAXLOOP,-i+N+1,j,ip,Qprime_ip,Q,QBI_SR1,QBI_SR1_init,QBI_SR1_add,Qprime,Qprime_body2,Qprime_body3,QM,QBI_SR2);
    #define S1(i,j,i2) Qprime(-i+N+1,j) = fillMatrices1_unafold_Qprime_body1_reduce_1(N,MAXLOOP,-i+N+1,j,Qprime_ip,Q,QBI_SR1,QBI_SR1_init,QBI_SR1_add,Qprime,Qprime_body2,Qprime_body3,QM,QBI_SR2);
    #define S2(i,j,i2) Qprime_body2() = fillMatrices1_unafold_Qprime_body2_reduce_1(N,MAXLOOP,-i+N+1,j,Qprime_ip,Q,QBI_SR1,QBI_SR1_init,QBI_SR1_add,Qprime,Qprime_body2,Qprime_body3,QM,QBI_SR2);
    #define S3(i,j,i2) Qprime_body3() = fillMatrices1_unafold_Qprime_body3_reduce_1(N,MAXLOOP,-i+N+1,j,Qprime_ip,Q,QBI_SR1,QBI_SR1_init,QBI_SR1_add,Qprime,Qprime_body2,Qprime_body3,QM,QBI_SR2);
    #define S4(i,j,i2) QM(-i+N+1,j) = fillMatrices1_unafold_QM_body_reduce_1(N,MAXLOOP,-i+N+1,j,Qprime_ip,Q,QBI_SR1,QBI_SR1_init,QBI_SR1_add,Qprime,Qprime_body2,Qprime_body3,QM,QBI_SR2);
    #define S5(i,j,ip) QBI_SR1_init(-i+N+1,j,ip) = fillMatrices1_unafold_QBI_SR1_init_reduce_1(N,MAXLOOP,-i+N+1,j,ip,Qprime_ip,Q,QBI_SR1,QBI_SR1_init,QBI_SR1_add,Qprime,Qprime_body2,Qprime_body3,QM,QBI_SR2);
    #define S6(i,j,ip) QBI_SR1_add(-i+N+1,j,ip) = fillMatrices1_unafold_QBI_SR1_add_reduce_1(N,MAXLOOP,-i+N+1,j,ip,Qprime_ip,Q,QBI_SR1,QBI_SR1_init,QBI_SR1_add,Qprime,Qprime_body2,Qprime_body3,QM,QBI_SR2);
    #define S7(i,j,i2) QM(-i+N+1,j) = INFINITY_VAL(0);
    #define S8(i,j,i2) QM(-i+N+1,j) = QM(-i+N+1,j);
    #define S9(i,j,ip) QBI_SR1(-i+N+1,j,ip) = QBI_SR1_init(-i+N+1,j,ip);
    #define S10(i,j,ip) QBI_SR1(-i+N+1,j,ip) = min(QBI_SR1_add(-i+N+1,j,ip),QBI_SR1(-i+N+2,j-1,ip));
    #define S11(i,j,i2) Qprime(-i+N+1,j) = INFINITY_VAL(0);
    #define S12(i,j,i2) Qprime(-i+N+1,j) = (((Eval_isFinite(Qprime_ip(-i+N+1,j)))>(0))?(min(Eh(-i+N+1,j),min((Es(-i+N+1,j))+(Qprime(-i+N+2,j-1)),min((((N-i-j+7>= 0 && -N+i+j-5>= 0))?INFINITY_VAL(0):(((N-i-j+10>= 0 && -N+i+j-8>= 0))?Qprime(-i+N+1,j):(min(Qprime_body2(),(Qprime_body3())+(Ebi_stacking(-i+N+1,j)))))),(((a(0))+(c(0)))+(End(-i+N+1,j)))+(QM(-i+N+2,j-1)))))):(INFINITY_VAL(0)));
    #define S13(i,j,i2) Q(-i+N+1,j) = INFINITY_VAL(0);
    #define S14(i,j,i2) Q(-i+N+1,j) = min((b(0))+(Q(-i+N+2,j)),min((b(0))+(Q(-i+N+1,j-1)),min(((c(0))+(End(-i+N+1,j)))+(Qprime(-i+N+1,j)),QM(-i+N+1,j))));
    {
        //Domain
    //{i,j,ip|N-8>= 0 && MAXLOOP-8>= 0 && i-5>= 0 && -N+i+j-ip-7>= 0 && j-4>= 0 && N-j+ip+3>= 0 && N+ip-2>= 0 && N-j>= 0 && N-ip-1>= 0 && j-ip-3>= 0 && i-ip-5>= 0 && N+MAXLOOP-i-j+ip+3>= 0 && ip-4>= 0 && N-i>= 0}
    //{i,j,i2|-N+i2== 0 && N-i>= 0 && MAXLOOP-8>= 0 && N-i-j+10>= 0 && -N+i+j-8>= 0 && N-j>= 0}
    //{i,j,i2|-N+i2== 0 && N-j>= 0 && MAXLOOP-8>= 0 && N-i>= 0 && -N+i+j-11>= 0}
    //{i,j,i2|-N+i2== 0 && -N+i+j-11>= 0 && MAXLOOP-8>= 0 && N-i>= 0 && N-j>= 0 && N-8>= 0}
    //{i,j,i2|-N+i2== 0 && N-j>= 0 && MAXLOOP-8>= 0 && N-i>= 0 && -N+i+j-10>= 0 && N-8>= 0}
    //{i,j,ip|N-j>= 0 && MAXLOOP-8>= 0 && N-i-j+ip+8>= 0 && ip-4>= 0 && -N+i+j-ip-7>= 0 && N-i>= 0}
    //{i,j,ip|N-j>= 0 && MAXLOOP-8>= 0 && N+MAXLOOP-i-j+ip+3>= 0 && ip-4>= 0 && -N+i+j-ip-9>= 0 && N-i>= 0}
    //{i,j,i2|-N+i2== 0 && N-8>= 0 && MAXLOOP-8>= 0 && N-i-j+9>= 0 && N-i>= 0 && N-j>= 0 && i-1>= 0 && j-2>= 0}
    //{i,j,i2|-N+i2== 0 && N-8>= 0 && MAXLOOP-8>= 0 && -N+i+j-10>= 0 && N-i>= 0 && N-j>= 0 && j-2>= 0 && i-1>= 0}
    //{i,j,ip|N-8>= 0 && MAXLOOP-8>= 0 && N-i-j+ip+8>= 0 && N-j>= 0 && ip-4>= 0 && N-i>= 0 && -N+i+j-ip-7>= 0 && N+MAXLOOP-i-j+ip+3>= 0}
    //{i,j,ip|N-8>= 0 && MAXLOOP-8>= 0 && -N+i+j-ip-9>= 0 && N-j>= 0 && N+MAXLOOP-i-j+ip+3>= 0 && N-i>= 0 && ip-4>= 0}
    //{i,j,i2|-N+i2== 0 && N-8>= 0 && MAXLOOP-8>= 0 && N-i-j+4>= 0 && N-i>= 0 && N-j>= 0 && i-1>= 0 && j-2>= 0}
    //{i,j,i2|-N+i2== 0 && N-8>= 0 && MAXLOOP-8>= 0 && -N+i+j-5>= 0 && N-j>= 0 && N-i>= 0}
    //{i,j,i2|-N+i2== 0 && N-8>= 0 && MAXLOOP-8>= 0 && N-i-j+4>= 0 && N-j>= 0 && N-i>= 0 && i-1>= 0 && j-2>= 0}
    //{i,j,i2|-N+i2== 0 && N-8>= 0 && MAXLOOP-8>= 0 && -N+i+j-5>= 0 && N-i>= 0 && N-j>= 0 && i-2>= 0 && j-3>= 0}
        int c1,c2,c3;
        for(c1=1;c1 <= 4;c1+=1){
            for(c2=2;c2 <= N;c2+=1){
                S7((c1),(c2),(N));
                S11((c1),(c2),(N));
                S13((c1),(c2),(N));
            }
        }
        for(c1=5;c1 <= 7;c1+=1){
            for(c2=2;c2 <= N+4-c1;c2+=1){
                S7((c1),(c2),(N));
                S11((c1),(c2),(N));
                S13((c1),(c2),(N));
            }
            for(c2=N+5-c1;c2 <= N;c2+=1){
                S7((c1),(c2),(N));
                S12((c1),(c2),(N));
                S14((c1),(c2),(N));
            }
        }
        for(c1=8;c1 <= MIN(9,N);c1+=1){
            for(c2=2;c2 <= N+4-c1;c2+=1){
                S7((c1),(c2),(N));
                S11((c1),(c2),(N));
                S13((c1),(c2),(N));
            }
            for(c2=N+5-c1;c2 <= N+7-c1;c2+=1){
                S7((c1),(c2),(N));
                S12((c1),(c2),(N));
                S14((c1),(c2),(N));
            }
            for(c2=N+8-c1;c2 <= N;c2+=1){
                S1((c1),(c2),(N));
                S7((c1),(c2),(N));
                S12((c1),(c2),(N));
                S14((c1),(c2),(N));
            }
        }
        if(N>=10){
            for(c2=2;c2 <= N-6;c2+=1){
                S7((10),(c2),(N));
                S11((10),(c2),(N));
                S13((10),(c2),(N));
            }
            for(c2=N-5;c2 <= N-3;c2+=1){
                S7((10),(c2),(N));
                S12((10),(c2),(N));
                S14((10),(c2),(N));
            }
            for(c2=N-2;c2 <= N-1;c2+=1){
                S1((10),(c2),(N));
                S7((10),(c2),(N));
                S12((10),(c2),(N));
                S14((10),(c2),(N));
            }
            S1((10),(N),(N));
            S4((10),(N),(N));
            S8((10),(N),(N));
            S12((10),(N),(N));
            S14((10),(N),(N));
        }
        for(c1=11;c1 <= MIN(12,N);c1+=1){
            for(c2=2;c2 <= N+4-c1;c2+=1){
                S7((c1),(c2),(N));
                S11((c1),(c2),(N));
                S13((c1),(c2),(N));
            }
            for(c2=N+5-c1;c2 <= N+7-c1;c2+=1){
                S7((c1),(c2),(N));
                S12((c1),(c2),(N));
                S14((c1),(c2),(N));
            }
            for(c2=N+8-c1;c2 <= N+9-c1;c2+=1){
                S1((c1),(c2),(N));
                S7((c1),(c2),(N));
                S12((c1),(c2),(N));
                S14((c1),(c2),(N));
            }
            S1((c1),(N+10-c1),(N));
            S4((c1),(N+10-c1),(N));
            S8((c1),(N+10-c1),(N));
            S12((c1),(N+10-c1),(N));
            S14((c1),(N+10-c1),(N));
            for(c2=N+11-c1;c2 <= N;c2+=1){
                for(c3=4;c3 <= c1+c2-N-7;c3+=1){
                    S0((c1),(c2),(c3));
                    S5((c1),(c2),(c3));
                    S9((c1),(c2),(c3));
                }
                S2((c1),(c2),(N));
                S3((c1),(c2),(N));
                S4((c1),(c2),(N));
                S8((c1),(c2),(N));
                S12((c1),(c2),(N));
                S14((c1),(c2),(N));
            }
        }
        for(c1=13;c1 <= N;c1+=1){
            for(c2=2;c2 <= N+4-c1;c2+=1){
                S7((c1),(c2),(N));
                S11((c1),(c2),(N));
                S13((c1),(c2),(N));
            }
            for(c2=N+5-c1;c2 <= N+7-c1;c2+=1){
                S7((c1),(c2),(N));
                S12((c1),(c2),(N));
                S14((c1),(c2),(N));
            }
            for(c2=N+8-c1;c2 <= N+9-c1;c2+=1){
                S1((c1),(c2),(N));
                S7((c1),(c2),(N));
                S12((c1),(c2),(N));
                S14((c1),(c2),(N));
            }
            S1((c1),(N+10-c1),(N));
            S4((c1),(N+10-c1),(N));
            S8((c1),(N+10-c1),(N));
            S12((c1),(N+10-c1),(N));
            S14((c1),(N+10-c1),(N));
            for(c2=N+11-c1;c2 <= N+12-c1;c2+=1){
                for(c3=4;c3 <= c1+c2-N-7;c3+=1){
                    S0((c1),(c2),(c3));
                    S5((c1),(c2),(c3));
                    S9((c1),(c2),(c3));
                }
                S2((c1),(c2),(N));
                S3((c1),(c2),(N));
                S4((c1),(c2),(N));
                S8((c1),(c2),(N));
                S12((c1),(c2),(N));
                S14((c1),(c2),(N));
            }
            for(c2=N+13-c1;c2 <= N;c2+=1){
                for(c3=MAX(4,c1+c2-N-MAXLOOP-3);c3 <= c1+c2-N-9;c3+=1){
                    S0((c1),(c2),(c3));
                    S6((c1),(c2),(c3));
                    S10((c1),(c2),(c3));
                }
                for(c3=c1+c2-N-8;c3 <= c1+c2-N-7;c3+=1){
                    S0((c1),(c2),(c3));
                    S5((c1),(c2),(c3));
                    S9((c1),(c2),(c3));
                }
                S2((c1),(c2),(N));
                S3((c1),(c2),(N));
                S4((c1),(c2),(N));
                S8((c1),(c2),(N));
                S12((c1),(c2),(N));
                S14((c1),(c2),(N));
            }
        }
    }
    #undef S0
    #undef S1
    #undef S2
    #undef S3
    #undef S4
    #undef S5
    #undef S6
    #undef S7
    #undef S8
    #undef S9
    #undef S10
    #undef S11
    #undef S12
    #undef S13
    #undef S14
    //Memory Free
    free(QBI_SR1);
    free(QBI_SR2);
    free(Qprime_body2);
    free(Qprime_body3);
    free(QBI_SR1_init);
    free(QBI_SR1_add);
}
int fillMatrices1_unafold_QBI_SR2_reduce_1(long N,long MAXLOOP,int i_p,int j_p,int ip_p,int* Qprime_ip,int* Q,int* QBI_SR1,int* QBI_SR1_init,int* QBI_SR1_add,int* Qprime,int* Qprime_body2,int* Qprime_body3,int* QM,int* QBI_SR2){
    int reduceVar = INT_MAX;
#define S0(i,j,ip,jp) {int __temp__ = ((Ebi_stacking(jp,ip))+(Ebi_asymmetry(-i+ip-1,j-jp-1)))+(Qprime(ip,jp)); reduceVar = min(reduceVar,__temp__); }
    {
        //Domain
    //{i,j,ip,jp|-ip_p-ip+jp== 0 && -j_p+j== 0 && -i_p+i== 0 && N-8>= 0 && MAXLOOP-8>= 0 && N-i_p-4>= 0 && -i_p+j_p-ip_p-6>= 0 && j_p-4>= 0 && N-j_p+ip_p+3>= 0 && N+ip_p-2>= 0 && N-j_p>= 0 && N-ip_p-1>= 0 && j_p-ip_p-3>= 0 && N-i_p-ip_p-4>= 0 && i_p-1>= 0 && ip_p-4>= 0 && MAXLOOP+i_p-j_p+ip_p+2>= 0 && ip-1>= 0 && N-ip_p-ip>= 0 && -i_p+ip-4>= 0 && -j_p+ip_p+ip+3>= 0 && j_p-ip_p-ip-2>= 0 && N-ip>= 0 && ip_p+ip-2>= 0}
        int c3;
        for(c3=MAX(i_p+4,j_p-ip_p-3);c3 <= j_p-ip_p-2;c3+=1){
            S0((i_p),(j_p),(c3),(c3+ip_p));
        }
    }
    #undef S0
    return reduceVar;
}
int fillMatrices1_unafold_Qprime_body1_reduce_1(long N,long MAXLOOP,int i_p,int j_p,int* Qprime_ip,int* Q,int* QBI_SR1,int* QBI_SR1_init,int* QBI_SR1_add,int* Qprime,int* Qprime_body2,int* Qprime_body3,int* QM,int* QBI_SR2){
    int reduceVar = INT_MAX;
#define S0(i,j,ip,jp) {int __temp__ = (((-j+jp+2== 0 && -i+ip-1== 0 && i-1>= 0 && -i+j-7>= 0 && i-j+9>= 0 && N-j>= 0) || (-j+jp+1== 0 && -i+ip-2== 0 && i-1>= 0 && -i+j-7>= 0 && i-j+9>= 0 && N-j>= 0))?(Ebi_Bulge1(i,j,ip,jp))+(Qprime(ip,jp)):(((-i+ip-1== 0 && -i+j-7>= 0 && j-jp-3>= 0 && jp-2>= 0 && i-j+9>= 0 && N-j>= 0 && i-1>= 0))?(Ebi_Bulge(i,j,ip,jp,j-jp-1))+(Qprime(ip,jp)):(((-j+jp+1== 0 && N-j>= 0 && i-j+9>= 0 && N-ip>= 0 && -i+ip-3>= 0 && i-1>= 0 && -i+j-7>= 0))?(Ebi_Bulge(i,j,ip,jp,-i+ip-1))+(Qprime(ip,jp)):(((-j+jp+2== 0 && -i+ip-2== 0 && i-1>= 0 && -i+j-7>= 0 && i-j+9>= 0 && N-j>= 0))?(Ebi_iloop1x1(i,j,ip,jp))+(Qprime(ip,jp)):(((-i+jp-6== 0 && -i+ip-2== 0 && -i+j-9== 0 && i-1>= 0 && N-i-9>= 0))?(Ebi_iloop1x2(i,j,ip,jp))+(Qprime(ip,jp)):((Ebi_iloop2x1(i,j,ip,jp))+(Qprime(ip,jp)))))))); reduceVar = min(reduceVar,__temp__); }
    {
        //Domain
//{i,j,ip,jp|-j_p+jp+1== 0 && -j_p+j== 0 && -i_p+i== 0 && i_p-1>= 0 && MAXLOOP-8>= 0 && -i_p+j_p-7>= 0 && i_p-j_p+9>= 0 && N-j_p>= 0 && -i_p+ip-2>= 0 && N-ip>= 0} || {i,j,ip,jp|-j_p+j== 0 && -i_p+i== 0 && i_p-1>= 0 && MAXLOOP-8>= 0 && -i_p-j_p+ip+jp+1>= 0 && i_p-j_p+9>= 0 && N-j_p>= 0 && j_p-jp-2>= 0 && -i_p-j_p+2jp-3>= 0 && i_p-j_p-2ip+2jp+1>= 0} || {i,j,ip,jp|-i_p+ip-1== 0 && -j_p+j== 0 && -i_p+i== 0 && i_p-1>= 0 && MAXLOOP-8>= 0 && -i_p+j_p-7>= 0 && i_p-j_p+9>= 0 && N-j_p>= 0 && j_p-jp-3>= 0 && jp-2>= 0}
        int c3,c4;
        for(c3=i_p+1;c3 <= N;c3+=1){
            if(c3==i_p+1){
                for(c4=2;c4 <= j_p-3;c4+=1){
                    S0((i_p),(j_p),(i_p+1),(c4));
                }
            }
            for(c4=MAX(CDIV(2*c3+j_p-i_p-1,2),i_p+j_p-c3-1);c4 <= j_p-2;c4+=1){
                S0((i_p),(j_p),(c3),(c4));
            }
            if(c3>=i_p+2){
                S0((i_p),(j_p),(c3),(j_p-1));
            }
        }
    }
    #undef S0
    return reduceVar;
}
int fillMatrices1_unafold_Qprime_body2_reduce_1(long N,long MAXLOOP,int i_p,int j_p,int* Qprime_ip,int* Q,int* QBI_SR1,int* QBI_SR1_init,int* QBI_SR1_add,int* Qprime,int* Qprime_body2,int* Qprime_body3,int* QM,int* QBI_SR2){
    int reduceVar = INT_MAX;
#define S0(i,j,ip,jp) {int __temp__ = (((-j+jp+2== 0 && -i+ip-1== 0 && -i+j-10>= 0 && i-1>= 0 && N-j>= 0) || (-j+jp+1== 0 && -i+ip-2== 0 && -i+j-10>= 0 && i-1>= 0 && N-j>= 0))?(Ebi_Bulge1(i,j,ip,jp))+(Qprime(ip,jp)):(((-i+ip-1== 0 && i-1>= 0 && -i+j-10>= 0 && jp-2>= 0 && N-j>= 0 && j-jp-3>= 0))?(Ebi_Bulge(i,j,ip,jp,j-jp-1))+(Qprime(ip,jp)):(((-j+jp+1== 0 && -i+j-10>= 0 && N-j>= 0 && N-ip>= 0 && -i+ip-3>= 0 && i-1>= 0))?(Ebi_Bulge(i,j,ip,jp,-i+ip-1))+(Qprime(ip,jp)):(((-j+jp+2== 0 && -i+ip-2== 0 && -i+j-10>= 0 && i-1>= 0 && N-j>= 0))?(Ebi_iloop1x1(i,j,ip,jp))+(Qprime(ip,jp)):(((-j+jp+3== 0 && -i+ip-2== 0 && -i+j-10>= 0 && i-1>= 0 && N-j>= 0))?(Ebi_iloop1x2(i,j,ip,jp))+(Qprime(ip,jp)):(((-j+jp+2== 0 && -i+ip-3== 0 && -i+j-10>= 0 && i-1>= 0 && N-j>= 0))?(Ebi_iloop2x1(i,j,ip,jp))+(Qprime(ip,jp)):((Ebi_iloop2x2(i,j,ip,jp))+(Qprime(ip,jp))))))))); reduceVar = min(reduceVar,__temp__); }
    {
        //Domain
//{i,j,ip,jp|-j_p+jp+2== 0 && -i_p+ip-2== 0 && -j_p+j== 0 && -i_p+i== 0 && N-j_p>= 0 && MAXLOOP-8>= 0 && i_p-1>= 0 && -i_p+j_p-10>= 0} || {i,j,ip,jp|-j_p+j== 0 && -i_p+i== 0 && N-j_p>= 0 && MAXLOOP-8>= 0 && i_p-1>= 0 && -i_p+j_p-10>= 0 && -j_p+jp+3>= 0 && -i_p+j_p+ip-jp-5>= 0 && i_p-ip+3>= 0} || {i,j,ip,jp|-i_p+ip-1== 0 && -j_p+j== 0 && -i_p+i== 0 && N-j_p>= 0 && MAXLOOP-8>= 0 && i_p-1>= 0 && -i_p+j_p-10>= 0 && jp-2>= 0 && j_p-jp-2>= 0} || {i,j,ip,jp|-j_p+jp+1== 0 && -j_p+j== 0 && -i_p+i== 0 && N-j_p>= 0 && MAXLOOP-8>= 0 && i_p-1>= 0 && -i_p+j_p-10>= 0 && -i_p+ip-2>= 0 && N-ip>= 0}
        int c3,c4;
        for(c3=i_p+1;c3 <= N;c3+=1){
            if(c3==i_p+1){
                for(c4=2;c4 <= j_p-2;c4+=1){
                    S0((i_p),(j_p),(i_p+1),(c4));
                }
            }
            if(c3<=i_p+3){
                for(c4=j_p-3;c4 <= c3+j_p-i_p-5;c4+=1){
                    S0((i_p),(j_p),(c3),(c4));
                }
            }
            if(c3==i_p+2){
                S0((i_p),(j_p),(i_p+2),(j_p-2));
            }
            if(c3>=i_p+2){
                S0((i_p),(j_p),(c3),(j_p-1));
            }
        }
    }
    #undef S0
    return reduceVar;
}
int fillMatrices1_unafold_Qprime_body3_reduce_1(long N,long MAXLOOP,int i_p,int j_p,int* Qprime_ip,int* Q,int* QBI_SR1,int* QBI_SR1_init,int* QBI_SR1_add,int* Qprime,int* Qprime_body2,int* Qprime_body3,int* QM,int* QBI_SR2){
    int reduceVar = INT_MAX;
#define S0(i,j,d) {int __temp__ = (min(QBI_SR1(i,j,d),QBI_SR2(d)))+(Ebi_sizePenalty(-i+j-d-2)); reduceVar = min(reduceVar,__temp__); }
    {
        //Domain
    //{i,j,d|-j_p+j== 0 && -i_p+i== 0 && -i_p+j_p-10>= 0 && MAXLOOP-8>= 0 && i_p-1>= 0 && N-j_p>= 0 && N-8>= 0 && d-4>= 0 && -i_p+j_p-d-6>= 0 && MAXLOOP+i_p-j_p+d+2>= 0}
        int c3;
        for(c3=MAX(4,j_p-MAXLOOP-i_p-2);c3 <= j_p-i_p-6;c3+=1){
            S0((i_p),(j_p),(c3));
        }
    }
    #undef S0
    return reduceVar;
}
int fillMatrices1_unafold_QM_body_reduce_1(long N,long MAXLOOP,int i_p,int j_p,int* Qprime_ip,int* Q,int* QBI_SR1,int* QBI_SR1_init,int* QBI_SR1_add,int* Qprime,int* Qprime_body2,int* Qprime_body3,int* QM,int* QBI_SR2){
    int reduceVar = INT_MAX;
#define S0(i,j,k) {int __temp__ = (Q(i,k-1))+(Q(k,j)); reduceVar = min(reduceVar,__temp__); }
    {
        //Domain
    //{i,j,k|-j_p+j== 0 && -i_p+i== 0 && N-j_p>= 0 && MAXLOOP-8>= 0 && i_p-1>= 0 && -i_p+j_p-9>= 0 && N-8>= 0 && k-3>= 0 && N-i_p>= 0 && -i_p+k-4>= 0 && j_p-k-5>= 0 && N-k>= 0 && j_p-2>= 0}
        int c3;
        for(c3=i_p+4;c3 <= j_p-5;c3+=1){
            S0((i_p),(j_p),(c3));
        }
    }
    #undef S0
    return reduceVar;
}
int fillMatrices1_unafold_QBI_SR1_init_reduce_1(long N,long MAXLOOP,int i_p,int j_p,int ip_p,int* Qprime_ip,int* Q,int* QBI_SR1,int* QBI_SR1_init,int* QBI_SR1_add,int* Qprime,int* Qprime_body2,int* Qprime_body3,int* QM,int* QBI_SR2){
    int reduceVar = INT_MAX;
#define S0(i,j,ip,jp) {int __temp__ = ((Ebi_stacking(jp,ip))+(Ebi_asymmetry(-i+ip-1,j-jp-1)))+(Qprime(ip,jp)); reduceVar = min(reduceVar,__temp__); }
    {
        //Domain
    //{i,j,ip,jp|-ip_p-ip+jp== 0 && -j_p+j== 0 && -i_p+i== 0 && N-j_p>= 0 && MAXLOOP-8>= 0 && -i_p+j_p-ip_p-6>= 0 && i_p-j_p+ip_p+7>= 0 && ip_p-4>= 0 && i_p-1>= 0 && N-8>= 0 && j_p-ip-8>= 0 && j_p-ip_p-ip-4>= 0 && -i_p+ip-2>= 0 && -i_p+ip_p+ip-6>= 0}
        int c3;
        for(c3=i_p+2;c3 <= j_p-ip_p-4;c3+=1){
            S0((i_p),(j_p),(c3),(c3+ip_p));
        }
    }
    #undef S0
    return reduceVar;
}
int fillMatrices1_unafold_QBI_SR1_add_reduce_1(long N,long MAXLOOP,int i_p,int j_p,int ip_p,int* Qprime_ip,int* Q,int* QBI_SR1,int* QBI_SR1_init,int* QBI_SR1_add,int* Qprime,int* Qprime_body2,int* Qprime_body3,int* QM,int* QBI_SR2){
    int reduceVar = INT_MAX;
#define S0(i,j,ip,jp) {int __temp__ = ((Ebi_stacking(jp,ip))+(Ebi_asymmetry(-i+ip-1,j-jp-1)))+(Qprime(ip,jp)); reduceVar = min(reduceVar,__temp__); }
    {
        //Domain
//{i,j,ip,jp|-i_p-ip_p+jp-2== 0 && -i_p+ip-2== 0 && -j_p+j== 0 && -i_p+i== 0 && N-j_p>= 0 && MAXLOOP-8>= 0 && -i_p+j_p-ip_p-8>= 0 && MAXLOOP+i_p-j_p+ip_p+2>= 0 && ip_p-4>= 0 && i_p-1>= 0} || {i,j,ip,jp|-j_p+jp+4== 0 && -j_p+ip_p+ip+4== 0 && -j_p+j== 0 && -i_p+i== 0 && N-j_p>= 0 && MAXLOOP-8>= 0 && -i_p+j_p-ip_p-8>= 0 && MAXLOOP+i_p-j_p+ip_p+2>= 0 && ip_p-4>= 0 && i_p-1>= 0}
        S0((i_p),(j_p),(i_p+2),(i_p+ip_p+2));
        S0((i_p),(j_p),(j_p-ip_p-4),(j_p-4));
    }
    #undef S0
    return reduceVar;
}
#undef Qprime_ip
#undef Q
#undef QBI_SR1
#undef QBI_SR1_init
#undef QBI_SR1_add
#undef Qprime
#undef Qprime_body2
#undef Qprime_body3
#undef QM
#undef QBI_SR2
//Postamble
#undef MAX
#undef max
#undef MIN
#undef min
#undef CEILD
#undef ceild
#undef FLOORD
#undef floord
#undef CDIV
#undef FDIV
#undef LB_SHIFT
#undef MOD
#undef RADD
#undef RMUL
#undef RMAX
#undef RMIN
