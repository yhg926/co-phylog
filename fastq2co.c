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
#include <ctype.h>
#include <stdlib.h>
#include "ctxobj.h"
#include "basemap.h"
#include "comp_rvs.h"
//fastq macro
#define LEN 4096 //line length
#define SRA_LOWSTD  52 //62 
#define SRA_HIGHSTD 82 //72 //'H'
#define ILL_LOWSTD 84
#define ILL_HIGHSTD 104
#define HASHSIZE 20000081 //320001293*2 // 40000003 //20000081*2 //default hashsize
#define HLF 0.7 //up limit of loadfacor
#define CONFLITBIT 0x8000000000000000LLU     //HIBITSET1 0x8000000000000000LL
#define QBITSET1   0x4000000000000000LLU //set qualified tuple this bit 1

#define h1(k) ((k)%(HASHSIZE)) 
#define h2(k) (1+((k)%(HASHSIZE-1)))
#define h(k,i) ((h1(k) + i * h2(k))%HASHSIZE)

static int comp_bittl = 64-BITTL ;
static int crvsaddmove = CTXLEN*2; //bits the new base should right move before add to crvstuple
static llong tupmask = TUPMASK;
static llong objmask = OBJMASK;
static llong ctxmask = TUPMASK^OBJMASK;
static int TL = TUPLEN ;
static llong tuple,crvstuple;
//char buf[Titlelen];
char seq[LEN]; char qual[LEN]; //static int i,j;
static int bottom,top;
llong *CO; 
static int hashsize=HASHSIZE; //hashtable and size
static int keycountlimit=(int)HASHSIZE*HLF ;
static int keycount; 

static inline void fq2num_pre(void){	
	int tl = TL-1;
	int basecode;
	while(tl--){
		if( (basecode=Basemap[(int)seq[bottom++]])!=DEFAULT ) {
        	tuple=((tuple<<2)|basecode);
		}
		else if (top - bottom > TL)
			tl=TL-1;
		else {
			bottom=top;
			break;
		}	 
	}
	crvstuple = crvs64bits(tuple)>>comp_bittl;
}

static inline int FqhashInsert(){
    int i,j,n,crn ;
    for(i=0;i<hashsize;i++){
        n=((tuple&objmask)+i)%hashsize; //tuple&objmask: get context 
//		n=h((tuple&objmask),i);
        if (CO[n]==0) {
            for(j=0;j<hashsize;j++) {				
                crn=((crvstuple&objmask)+j)%hashsize;
				//crn=h((crvstuple&objmask),j);
                if (CO[crn]==0){
                    CO[n]= tuple|( QBITSET1>>((tuple&ctxmask)>>CTXLEN) ); 
                    keycount++;
					break;
                }
                else if ((CO[crn]&objmask)==(crvstuple&objmask)){
					if(CO[crn]&CONFLITBIT)
						break;
					llong crnflag=(QBITSET1>>((CO[crn]&ctxmask)>>CTXLEN));
					if( (CO[crn]&tupmask)==crvstuple ){
						CO[crn]&=~crnflag;
						break;
					}
					llong crvstupleflag=(QBITSET1>>((crvstuple&ctxmask)>>CTXLEN));
					CO[crn]^=crvstupleflag;
					if((CO[crn]&crvstupleflag)==0) {
						if((CO[crn]&crnflag)==0)
							CO[crn]|=CONFLITBIT;
						else
							CO[crn]=(CO[crn]&objmask)|crvstuple;						
					}
					break;
				} 					
            }
            break;
        }
        else if ( (CO[n]&objmask)==(tuple&objmask) ) {
                if (CO[n]&CONFLITBIT)
					break;
				llong nflag      =(QBITSET1>>((CO[n]&ctxmask)>>CTXLEN));
				if( (CO[n]&tupmask)==tuple ){
                    CO[n]&=~nflag;
                    break;
                }
				llong tupleflag=(QBITSET1>>((tuple&ctxmask)>>CTXLEN));
				CO[n]^=tupleflag;

				if((CO[n]&tupleflag)==0) {
                    if((CO[n]&nflag)==0)
                        CO[n]|=CONFLITBIT;
                    else
                        CO[n]=(CO[n]&objmask)|tuple;
                }
                break;
        }
   }
   return 1;
                        //end hashing
}


static inline int HQ_FqhashInsert(){
    int i,j,n,crn ;
    for(i=0;i<hashsize;i++){
        n=((tuple&objmask)+i)%hashsize; //tuple&objmask: get context 
        //n=h((tuple&objmask),i);
		if (CO[n]==0) {
            for(j=0;j<hashsize;j++) {
                crn=((crvstuple&objmask)+j)%hashsize;
                //crn=h((crvstuple&objmask),j);
				if (CO[crn]==0){
                    CO[n]= tuple;//|( QBITSET1>>((tuple&ctxmask)>>CTXLEN) ); 
                    keycount++;
					break;
                }
                else if ((CO[crn]&objmask)==(crvstuple&objmask)){
                    if(CO[crn]&CONFLITBIT)
                        break;
                    llong crnflag=(QBITSET1>>((CO[crn]&ctxmask)>>CTXLEN));
                    if( (CO[crn]&tupmask)==crvstuple ){
                        CO[crn]&=~crnflag;
                        break;
                    }
                    llong crvstupleflag=(QBITSET1>>((crvstuple&ctxmask)>>CTXLEN));
                    //CO[crn]^=crvstupleflag;
                    if((CO[crn]&crnflag)==0) 
                            CO[crn]|=CONFLITBIT;
                    else
                        CO[crn]= ( (CO[crn]&objmask)|crvstuple )&~crvstupleflag  ;                    
					break;
                }
            }
            break;
        }
        else if ( (CO[n]&objmask)==(tuple&objmask) ) {
                if (CO[n]&CONFLITBIT)
                    break;
                llong nflag      =(QBITSET1>>((CO[n]&ctxmask)>>CTXLEN));
                if( (CO[n]&tupmask)==tuple ){
                    CO[n]&=~nflag;
                    break;
                }
                llong tupleflag=(QBITSET1>>((tuple&ctxmask)>>CTXLEN));
                //CO[n]^=tupleflag;
                if((CO[n]&nflag)==0)
                    CO[n]|=CONFLITBIT;
                else
                    CO[n]= ( (CO[n]&objmask)|tuple )&~tupleflag  ;                
                break;
        }
   }
   return 1;
                        //end HQhashing
}

int main(int argc, char* argv[])
{
	if(argc != 3){
	    printf("USEAGE: ./a.out <*.fq> <*.co>\n");
		exit(1);
    }
    FILE *infp, *outfp;
	if((infp=fopen(argv[1],"r"))==NULL){
		printf("can't open argv[1]\n");
		exit(1);
	}
	else if((outfp=fopen(argv[2],"wb"))==NULL){
		printf("can't open argv[2]\n");
		exit(1);
	}
	//allocate memory for CO table
	hashsize=HASHSIZE;// 20000081;
	if ( (CO = calloc(hashsize, sizeof(llong)) ) == NULL ) {
        printf("allocate memory failed!");
        exit(1);
    }
	basemapset(Basemap);//set Basemap 
	/*Convert sequence to llong tuple*/
	llong basenum; char *buf;
	while(fgets(seq,LEN,infp)!=NULL) {
		buf=fgets(seq,LEN,infp);	
		buf=fgets(qual,LEN,infp);
		buf=fgets(qual,LEN,infp);

	    top=bottom=0; 
		while(qual[top++]!='\n'){
			if (qual[top] < SRA_LOWSTD) {
				if (top - bottom >= TL){
                    /*-------------------------------*/
					fq2num_pre();
    				while (bottom < top) {
						if( (basenum=Basemap[(int)seq[bottom++]])!=DEFAULT ) {
							tuple=((tuple<<2)|basenum) & tupmask;
							crvstuple=(crvstuple>>2)+((basenum^3)<<crvsaddmove);
							//hashing
							FqhashInsert();
						}
						else if(top - bottom > TL)
							fq2num_pre();
						else break;
					}
                    /*------------------------------*/
				}
				bottom=top+1;	
			}
			else if (qual[top] > SRA_HIGHSTD) {
				int hqst=top;//high qual start
				while (qual[++top] > SRA_HIGHSTD) ;

				if (top - hqst >= TL){
					int hqtst=hqst+TL-2;//high qual tuple start 
					fq2num_pre();
		            while (bottom < top) {
                        if( (basenum=Basemap[(int)seq[bottom++]])!=DEFAULT ) {
                            tuple=((tuple<<2)|basenum) & tupmask;
                            crvstuple=(crvstuple>>2)+((basenum^3)<<crvsaddmove);

							 if(bottom > hqtst)
								HQ_FqhashInsert();
							 else
								FqhashInsert(); 
                        }
                        else if (top - bottom > TL)
                            fq2num_pre();
                        else break;
                    }
				}
			}			
		}
		//resize hash
		if(keycount >keycountlimit){	
			int oldhashsize=hashsize;
			hashsize*=4;
			keycountlimit = hashsize*HLF ;
			printf("old hash size %d is not enough, procedure began to resize %d memory",oldhashsize,hashsize);
			llong *newCO;
			if ( (newCO = calloc(hashsize, sizeof(llong)) ) == NULL ) {
				printf("allocate memory failed!");
				exit(1);
			}
			int i,j,n;
			for(i=0;i<oldhashsize;i++){
				if (CO[i]!=0){
					for(j=0;j<hashsize;j++){
						n=((CO[i]&objmask)+j)%hashsize;
						if (newCO[n]==0){
							newCO[n]=CO[i];
							break;
						}
					}
				}
			}
			free(CO);
			CO=newCO;
		}
		// resize end
	}
	fclose(infp);
	//put out
	int i;
	for(i=0;i<hashsize;i++)
		if((CO[i]!=0) && (CO[i]< CONFLITBIT) && ((CO[i]&(QBITSET1>>((CO[i]&ctxmask)>>CTXLEN)))==0))
			fwrite(CO+i,8,1,outfp);
	;

	fclose(outfp);
	//	printf("%d\t%d\n",keycount,keycountlimit);	
	return 1;
}



/*
inline int FqhashInsert(llong qmask){
	int i,j,n,crn ;
	for(i=0;i<hashsize;i++){
		n=( (tuple&objmask)+i )%hashsize; //tuple&objmask: get context 
		if (CO[n]==0) {
			for(j=0;j<hashsize;j++) {
				crn=( (crvstuple&objmask)+j )%hashsize;
				if (CO[crn]==0){
					CO[n]=tuple|qmask; //qmask set the 2nd bit 1 if qualified 0 if not already qualified
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
   return 1;	
                        //end hashing
}
*/














