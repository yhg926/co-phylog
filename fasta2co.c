/******************************************************************************
 
 Co-phylog: an assembly-free phylogenomic approach for closely related organisms
 
 Copyright(c) 2013 Huiguang Yi (yhg926@gmail.com)
 
	This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>
 
 *************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "ctxobj.h"
#include "basemap.h"
#include "comp_rvs.h"
#define HIBITSET1 0x8000000000000000LL
#define HASHSIZE 20000081 //320001293 //40000003 
#define h1(k) ((k)%(HASHSIZE)) 
#define h2(k) (1+((k)%(HASHSIZE-1)))
#define h(k,i) ((h1(k) + i * h2(k))%HASHSIZE)

static int comp_bittl = 64-BITTL ;
static int crvsaddmove = CTXLEN*2; //bits the new base should right move before add to crvstuple
static llong tupmask = TUPMASK;
static llong objmask = OBJMASK;
static int TL = TUPLEN ; 
static llong tuple,crvstuple;

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
//	int hashsize = statbuf.st_size/LF >HASHSIZE?statbuf.st_size/LF:HASHSIZE ;
	int hashsize =HASHSIZE;
	printf("Size of file in bytes: %d\n", hashsize);

	if ( (CO = calloc(hashsize, sizeof(llong)) ) == NULL ) {
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
		int i,j,n,crn ;

    	for(i=0;i<hashsize;i++){
        	//n=( (tuple&objmask)+i )%hashsize; //tuple&objmask: get context 
      		n=h((tuple&objmask),i);
			if (CO[n]==0) {
            	for(j=0;j<hashsize;j++) {
                	//crn=( (crvstuple&objmask)+j )%hashsize;
               		crn=h((crvstuple&objmask),j);
					if (CO[crn]==0){
                    	CO[n]=tuple;
                   		break;
                	}
               		else if ( (CO[crn]&objmask)==(crvstuple&objmask) ){
                    	if(CO[crn]!=crvstuple)
							CO[crn]|=HIBITSET1;

                     	break;                  
               		}
            	}
            	break;
        	}
        	else if ( (CO[n]&objmask)==(tuple&objmask) ) {
            	if(CO[n]!=tuple)
					CO[n]|=HIBITSET1;

            	break;          
        	}
    	}
		//end hashing
    }
	fclose(infp);

	int count;
	for(count=0;count<hashsize;count++)
		if(CO[count]!=0 && CO[count]<0x8000000000000000LLU)
			fwrite(CO+count,8,1,outfp);
	;
	fclose(outfp);
	/*
	FILE *fp;
	if((fp=fopen("./testfile","rb"))==NULL)
		printf("can't open testfile\n");
	llong read[20];
	if(fread(read,8,20,fp))
		for(count=0;count<20;count++)
			printf("%llx\n",read[count]);
	*/

	return 1;
}






















