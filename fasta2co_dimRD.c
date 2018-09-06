#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <math.h>
#include <string.h>
#include "ctxobj.h"
#include "basemap.h"
#include "comp_rvs.h"
#include "dim_rdc.h"
#define HIBITSET1 0x8000000000000000LLU
#define HASHSIZE 4000003
#define h1(k) ((k)%(HASHSIZE))
#define h2(k) (1+((k)%(HASHSIZE-1)))
#define h(k,i) ((h1(k) + i * h2(k))%HASHSIZE)

static int comp_bittl = 64-BITTL ;
static int crvsaddmove = CTXLEN*2; //bits the new base should right move before add to crvstuple
static llong tupmask = TUPMASK;
static llong objmask = OBJMASK;
static int TL = TUPLEN ;
static llong tuple,crvstuple,unituple,unituple2,drtuple,pfilter;
//MASK for Dimension Reduction
static llong domask = ( (1LLU << (DO*2)) - 1 ) << (CTXLEN+2);  
static llong undomask =  ( (1LLU << ((HALFCTX-DO)*2)) - 1 ) << (CTXLEN+2+DO*2); 

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


int main(int argc, char* argv[]){
    if(argc!=3){
        printf("USAGE: ./fastaCovert2co <*.fasta/*.fna> <*.co> \n");
        exit(1);
    }
    if(DRR!=PD/PSTDR){
        printf("NOT correct dimension reduction rate \n"); 
        exit(1);
    }else if(PD!=pow(16,DO)){
        printf("NOT correct dimension order\n");
        exit(1);
    }else if (1<<DRbits != PSTDR){
        printf("NOT correct bits\n");
        exit(1);
    }
    
   
    if( access("./cophylog.seed",F_OK) == -1 ){
        //file not exists
        FILE *seedout;
        seedout=fopen("./cophylog.seed","wb");
        short pr_dm[PSTDR]; char meet[PD]={0}; int tmp;
        for (int i=0;i < PSTDR; i++ ){
            while( meet[tmp=(int)rand()%PD]!=0 ){};
            pr_dm[i]=tmp;
            meet[tmp]=1;
        }
        short dmp=PSTDR,dmi=PD;
         
        fwrite(&dmp,sizeof(short),1,seedout);
        
        fwrite(&dmi,sizeof(short),1,seedout);
        fwrite(pr_dm,sizeof(short),PSTDR,seedout);
        fclose(seedout);
        printf("%d fold Dimension reduction from %d to %d\n",DRR, PD,PSTDR); 
    }
    
    FILE *seedfp;
    seedfp=fopen("./cophylog.seed","rb");
    short dmp_r,dmi_r;
    fread(&dmp_r,sizeof(short),1,seedfp);
    fread(&dmi_r,sizeof(short),1,seedfp);
    if(dmp_r != PSTDR){
        printf("dimension from cophylog.seed %d is inconsist with %d \n",dmp_r,PSTDR); 
        exit(1);
    }else if(dmi_r != PD){
        printf("dimension from cophylog.seed %d is inconsist with %d \n",dmi_r,PD); 
        exit(1);
    }
    printf("infile: %d fold Dimension reduction from %d to %d\n",DRR, PD,PSTDR); 
    short pr_dm[PSTDR], rd_dm[PD];
    //set defaut value as -1
    memset(rd_dm,-1,sizeof(rd_dm));
    fread(pr_dm,sizeof(short),PSTDR,seedfp);
    for (int i=0;i<PSTDR;i++)
        rd_dm[pr_dm[i]]=i; 
    //for(int i=0 ; i< PD;i++){
    //   printf("%d %d\n",i,rd_dm[i]);
    //}
//-------------//
    basemapset(Basemap);
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
    llong *CO;              //CO is the hashtable
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
        //the lexi order will not be determined by object
        if ((tuple&objmask) == (crvstuple&objmask))
            continue;
        //project ctx to a uniq number
        if(tuple > crvstuple ){
            unituple=crvstuple;
            unituple2=tuple;
        }else{
            unituple=tuple;
            unituple2=crvstuple;
        }
        //only for 64bits Kmer storage && length(obj)=1, 
        short dim_tup = ((unituple & domask) >> (CTXLEN+2 - DO*2 )) + ( (unituple2 & domask) >> (CTXLEN+2));
        //DEAULT is -1 , which will be filtered
   /* 
        for(int i=0;i<TUPLEN;i++)
				printf("%c", Mapbase[(unituple<<(64-2*(TUPLEN-i)))>>62]) ;
			printf("\n");
            for(int i=0;i<TUPLEN;i++)
				printf("%c", Mapbase[(unituple2<<(64-2*(TUPLEN-i)))>>62]) ;
			printf("\n");

        printf("value:%d\n",rd_dm[dim_tup]);
 */
        if(rd_dm[dim_tup]==-1)
            continue;
        pfilter=rd_dm[dim_tup];//post filter
        drtuple = (((unituple & undomask) + ((unituple2 & undomask) >> ((HALFCTX-DO)*2))) >> (DO*4-DRbits))
        + (pfilter << 2) +  ((unituple & CTXMASK) >> CTXLEN);
 
      /*  for(int i=0;i<TUPLEN;i++)
			printf("%c", Mapbase[(drtuple<<(64-2*(TUPLEN-i)))>>62]) ;
		printf("\n");
*/
        int i,n;
        for(i=0;i<HASHSIZE;i++){
            //n=( (tuple&objmask)+i )%hashsize; //tuple&objmask: get context
            n=h((drtuple>>2),i);
            if (CO[n]==0) {
                CO[n]=drtuple;
                break;
            }
            else if ( (CO[n]>>2)==(drtuple>>2) ) {
                if(CO[n]!=drtuple) //  if(CO[n]!=tuple) ?? bug find? 20180302
                    CO[n]|=HIBITSET1;
                break;
            }
            else if(CO[n] > HIBITSET1)
                break;
        }
        //end hashing
    }
    fclose(infp);
    int count;
    for(count=0;count<HASHSIZE;count++)
        if(CO[count]!=0 && CO[count]< HIBITSET1)
            fwrite(CO+count,sizeof(llong),1,outfp);
    ;
    printf("%d sites writed\n",count);
    fclose(outfp);
    return 1;
}
