#ifndef BASEMAP_H
#define BASEMAP_H

#include <string.h>
#define DEFAULT -1
int Basemap[128];

static int *basemapset(int map[128]) {
	memset(map, DEFAULT, 128*sizeof(int));
	map['A']=0; map['a']=0;
	map['C']=1; map['c']=1;
	map['G']=2; map['g']=2;
	map['T']=3; map['t']=3;
	return map;
}

char Mapbase[]={'A','C','G','T'};


#endif
