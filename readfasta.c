/*#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#define CL 18
#define TL CL+1 
typedef unsigned long long int llong;
char gch; llong tuple;
*/
#include "readfasta.h"

static inline void prepare(FILE *fp){
	int base=0;tuple=0;
	while((gch=fgetc(fp)) && base < TL){
		if	(gch=='A')		{tuple<<=2 ;base++;}
		else if (gch=='C')	{tuple= (tuple<<2)|1 ;base++;}
		else if (gch=='G')	{tuple= (tuple<<2)|2 ;base++;}
		else if (gch=='T')	{tuple= (tuple<<2)|3 ;base++;}
		else if (gch=='a')  {tuple<<=2   ;base++;}
		else if (gch=='c')  {tuple= (tuple<<2)|1 ;base++;}
		else if (gch=='g')  {tuple= (tuple<<2)|2 ;base++;}
		else if (gch=='t')  {tuple= (tuple<<2)|3 ;base++;}
		else if (gch=='\n') continue;
		else if (gch=='\r') continue;
		else if (isalpha(gch)) {tuple=0;base=0;continue;}
		else if (gch=='>') {while( ((gch=fgetc(fp))!='\n') && gch != EOF ) ;base=0; tuple=0; continue;}
		else if (gch==EOF) break;
		//these cases usally don't happen
		else if (gch<33) continue;// the same with '\n'
		else if ((gch <= 128) && (gch >= 0) ){tuple=0;base=0;continue;} // the same with isalpha
		else {printf("illegal charactor!\n");exit(1);}
	}
}

//llong tmp[151689122];int i;

static inline  void ch_work(FILE *infp,FILE *outfp){
//	char gch='\0'; llong tuple=0;
//	llong out=0;
	 prepare(infp);
	 while ((gch=fgetc(infp))){
		 if  (gch=='A')      {tuple<<=2          ;}
         else if (gch=='C')  {tuple= (tuple<<2)|1;}
         else if (gch=='G')  {tuple= (tuple<<2)|2;}
         else if (gch=='T')  {tuple= (tuple<<2)|3;}
         else if (gch=='a')  {tuple<<=2          ;}
         else if (gch=='c')  {tuple= (tuple<<2)|1;}
         else if (gch=='g')  {tuple= (tuple<<2)|2;}
         else if (gch=='t')  {tuple= (tuple<<2)|3;}
         else if (gch=='\n') continue;
         else if (gch=='\r') continue;
         else if (isalpha(gch)) {prepare(infp);continue;}
         else if (gch=='>') {while( ((gch=fgetc(infp))!='\n') && gch != EOF ) ;prepare(infp); continue;}
         else if (gch==EOF) break;
         //usally don't happen
         else if (gch<33) continue;// the same with '\n'
		 else if ((gch <= 128) && (gch >= 0) ){prepare(infp);continue;} // the same with isalpha
		 else {printf("illegal charactor!\n");exit(1);}
		//	~revllong(tuple,BITTL) & HIGHMASK; 
	     //printf("%llu\t\n",tuple&HIGHMASK);// ~revllong(tuple,BITTL) & HIGHMASK );//) & HIGHMASK) ; 
//		 tmp[i++]=(tuple&HIGHMASK);
//	 	 fwrite(tmp,8,151689122,outfp);
//		 out= ~revllong(tuple,BITTL) & HIGHMASK;
//		 fwrite(&out, 8, 1,outfp);
	 }
//	fwrite(tmp,8,151689122,outfp);

}

/*
static inline llong revllong(llong tuple, int length){
	llong revtuple = 0;
	while(length--){
		 revtuple+=(tuple&BIT1MASK);
		 revtuple<<=1;
		 tuple>>=1;
	}
	return revtuple;
} 

*/








