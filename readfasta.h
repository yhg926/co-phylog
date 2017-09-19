#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#define HALFCXT 9
#define CL 2*(HALFCXT)
#define TL (CL)+1
#define BITCL 2*(CL)
#define BITTL 2*(TL)
#define HIBITTL 64-(BITTL)
#define _64MASK 0xffffffffffffffffLLU
#define HIGHMASK (_64MASK) >> (HIBITTL)
#define BIT1MASK 0x0000000000000001LLU
typedef unsigned long long int llong;
char gch; llong tuple;

static inline void prepare(FILE *);
//static inline llong revllong(llong , int);
static inline void ch_work(FILE *infp, FILE *outfp);
