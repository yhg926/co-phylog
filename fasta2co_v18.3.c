/*  line107: find and correct vital bug in old version 
	if(CO[n]!=unituple) //  if(CO[n]!=tuple) ?? bug find? 20180302
*/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "ctxobj.h"
#include "basemap.h"
#include "comp_rvs.h"
#define HIBITSET1 0x8000000000000000LLU
#ifndef HASHSIZE
#define HASHSIZE 40000003
#endif  
#define h1(k) ((k)%(HASHSIZE)) 
#define h2(k) (1+((k)%(HASHSIZE-1)))
#define h(k,i) ((h1(k) + i * h2(k))%HASHSIZE)

static int comp_bittl = 64-BITTL ;
static int crvsaddmove = CTXLEN*2; //bits the new base should right move before add to crvstuple
static llong tupmask = TUPMASK;
static llong objmask = OBJMASK;
static int TL = TUPLEN ; 
static llong tuple,crvstuple,unituple ;//unituple is unified tuple either tuple or crvstuple according to their value

static inline int fst_pre(FILE *fp){
	int base=1; char ch;
	while( (base < TL)&&((ch=fgetc(fp))!=EOF) ){
		 if(Basemap[(int)ch]!=DEFAULT) {tuple= (tuple<<2)|Basemap[(int)ch];base++;}
		 else if (ch=='\n') continue;
		 else if (ch=='\r') continue;
	     else if (isalpha(ch)) {base=1;continue;}
	     else if (ch=='>') {while( ((ch=fgetc(fp))!='\n') && ch != EOF ); base=1; continue;}
	     else {printf("ignorn illegal charactor%c!\n",ch); base=1; continue;}		
	}
	crvstuple = crvs64bits(tuple)>>comp_bittl ;
	return 1;
}

int main(int argc, char* argv[])
{
	basemapset(Basemap);//set Basemap 
	if(argc!=3){
		printf("USAGE: ./fastaCovert2co <*.fasta/*.fna> <*.co> \n");
		exit(1);
	}

	FILE *infp, *outfp;
	if( (infp=fopen(argv[1],"r"))==NULL ){
		printf("can't open argv[1]!\n");
		exit(1);
	}
	else if ( (outfp=fopen(argv[2],"wb"))==NULL ){
		printf("can't open argv[2]!\n");
		exit(1);
	}
	//allocate memory for CO table
	llong *CO;				//CO is the hashtable
	struct stat statbuf;
	stat(argv[1], &statbuf);
	printf("Size of file in bytes: %d\n", HASHSIZE);

	if ( (CO = calloc(HASHSIZE,sizeof(llong)) ) == NULL ) {
		printf("allocate memory failed!");
		exit(1);
	}
	/*Convert sequence to llong tuple*/
	fst_pre(infp);	
	llong basenum;
	char gch;

    while ( (gch=fgetc(infp))!=EOF ){

        if( (basenum=Basemap[(int)gch])!=DEFAULT ){
			tuple=((tuple<<2)|basenum) & tupmask;
			crvstuple=(crvstuple>>2)+((basenum^3)<<crvsaddmove);
		}			
        else if ((gch=='\n')|(gch=='\r'))
		   	continue;
        else if (isalpha(gch)) {
			fst_pre(infp);
			continue;
		}
        else if (gch=='>') {
			while( ((gch=fgetc(infp))!='\n') && gch != EOF );
			fst_pre(infp);
		   	continue;
		}
        else {
			printf("ignorn illegal charactor %c!\n",gch);
			fst_pre(infp);
			continue;
		}
		//begin hashing by linear probe   //FstTupleInsert(tuple,crvstuple);
		int i,n;
		if ((tuple&objmask) == (crvstuple&objmask))
			continue;
		unituple = tuple > crvstuple ? crvstuple : tuple ; // 2012.10.28 important improvment! only need consider tuple on either strand 
    	for(i=0;i<HASHSIZE;i++){
        	//n=( (tuple&objmask)+i )%hashsize; //tuple&objmask: get context 
      		n=h((unituple&objmask),i);
			if (CO[n]==0) {
                CO[n]=unituple;
                break;
        	}
        	else if ( (CO[n]&objmask)==(unituple&objmask) ) {
            	if(CO[n]!=unituple) //  if(CO[n]!=tuple) ?? bug find? 20180302
					CO[n]|=HIBITSET1;
            	break;          
        	}
    	}
		//end hashing
    }
	fclose(infp);

	int count;
	for(count=0;count<HASHSIZE;count++)
		if(CO[count]!=0 && CO[count]< HIBITSET1)
			fwrite(CO+count,8,1,outfp);
	;
	fclose(outfp);
	return 1;
}






















