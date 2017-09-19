#ifndef CRVS_H
#define CRVS_H

#include "ctxobj.h"
#define SWAP2  0x3333333333333333ULL
#define SWAP4  0x0F0F0F0F0F0F0F0FULL
#define SWAP8  0x00FF00FF00FF00FFULL
#define SWAP16 0x0000FFFF0000FFFFULL
#define SWAP32 0x00000000FFFFFFFFULL  


static inline llong crvs64bits(llong n) {

	n = ((n >> 2 ) & SWAP2 ) | ((n & SWAP2 ) << 2 );
	n = ((n >> 4 ) & SWAP4 ) | ((n & SWAP4 ) << 4 );
	n = ((n >> 8 ) & SWAP8 ) | ((n & SWAP8 ) << 8 );
	n = ((n >> 16) & SWAP16) | ((n & SWAP16) << 16);
	n = ((n >> 32) & SWAP32) | ((n & SWAP32) << 32);
	return ~n;
}

#endif








