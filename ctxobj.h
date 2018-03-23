#ifndef CTXOBJ_H
#define CTXOBJ_H

#define HALFCTX 9 // 1/2 context length
#define CTXLEN (2*HALFCTX)
#define TUPLEN (CTXLEN+1)
#define BITTL (2*TUPLEN)
//#define HIBITTL (64-BITTL)
#define _64MASK 0xffffffffffffffffLLU
#define TUPMASK (_64MASK >> (64-BITTL))
#define BIT1MASK 0x0000000000000001LLU
#define CTXMASK (3LLU << CTXLEN)
#define OBJMASK (TUPMASK^CTXMASK)
#define LF 0.5             //load factor

#define LINEARPROBE 1
#if LINEARPROBE==0
#define Pr 20000081 //size of hash  //prime around tuple counts/load factor
#define LPr Pr-1
#define h1(k) ((k)%(Pr)) 
#define h2(k) (1+((k)%(LPr)))
#define h(k,i) ((h1(k) + i * h2(k))%Pr)
#endif 
//linear probe: take advantage of cache line; faster than above method
//#define h(k,i) ((k+i)%Pr) 

typedef unsigned long long int llong;

//static int bittl = BITTL;
//static int hibittl = HIBITTL;
//static llong highmask = HIGHMASK;
//static llong objmask =  HIGHMASK^(3LLU << CTXLEN) ;
//static int TL = TUPLEN ;
#endif
