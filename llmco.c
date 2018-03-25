/* *************************massive-co.v18_323***************************************
// highly parallelable pairwise co-distances computing method using linked list 
// (need only < 1G memory per hashtable for about 64 genomes)
//                  contact: yhg926@gmail.com  
//                   beta version: 2018-3-24
// ********************************************************************************/
#include "comp_rvs.h"
#include "ctxobj.h"
#include <dirent.h>
#include <err.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
//set remoted related taxa as unrelated cluster: defined by co-distance >ULB, or #common ctx < SCLB
//for which a unified distance UNRD is set 
//co-dist < ULB && (#shared ctx > SCOL)&& jaccard index >JCLB or co-distance = URD
#define URD 0.5 // unrelated taxa distance: defined distace to URD
#define SUBC 0.25 //subcluster distance set to SUBC
#define ULB 0.09  //unrelated taxa boundary: based on C9,9O1 co-distance 
#define SCLB 1000 // for 3M # share ctx should >150K, for extreme small genome,need impose a JCLB
#define JCLB 0.08  //shared ctx/smaller genome or jaccard index lower boundary
//the optimized hash methods and parametors were chose here based on a test on 358 Salmonella genomes(20180323)
#define HASHSIZE 400000009 //a prime 400000009
#define NLF 0.5 
#define h1(k) ((k)%(HASHSIZE)) 
#define h2(k) (1+((k)%(HASHSIZE-1)))
#define h(k,i) ((h1(k) + i * h2(k))%HASHSIZE)//test: much faster than h(k,i)=(h1(k)+i(i+1)/2)%HASHSIZE
#define BTS 1  //bytes store SP.count+obj(2bits)
#define GCB 6 //genome code bits GCB=BTS*8-2
#define HKL 64 //maximum #sp. allowed in a hashtable,HKL = 2^GCB 
#define OBJAR 16 //obj array size
#define SHIFT_GC 52 //high 64-HSPB bits for genome code
#define SHIFT_MOBJ 40 //for count of other genmes match this ctx

static llong objmask = OBJMASK;
//static llong tupmask =TUPMASK;
static llong ctxmask = CTXMASK;
typedef struct objs{char gcobj[OBJAR]; struct objs *next;} objs_t; //gcobj[i]:genomecode:6bits,obj:2bits
struct entry {llong key; objs_t *new;} entry;
//key bits: genome_code|ctxcount|kmer:  12|12|40
void help(int exit_code)
{
	static const char str[] = "Usage: llmcodist DIR\n";
	fprintf(exit_code == 0 ? stdout : stderr, str);
	exit(exit_code);
}

int main(int argc, char *argv[])
{
#ifdef _OMP
	THREADS = omp_get_num_procs();
#endif

	int index;
	while ((index = getopt(argc, argv, "h"))) {
		if (index == -1) break;
		else if (index == 'h') help(0);
		else help(1);
	}

	argc -= optind;
	argv += optind;

	const char *dirpath = argv[0];
	if (argc != 1) {
		help(1);
	}
	DIR *dh;
	if ((dh = opendir(dirpath)) == NULL) {
		err(errno, "can't open %s\n", dirpath);
	}
	char fullname[512] = {0};
	char cofile[HKL][512]={""};
	llong fsize[HKL]={0};
	struct stat statbuf;
	struct dirent *dirent;
	int n = 0;
	llong size=0;
	while ((dirent = readdir(dh)) != NULL && (n < HKL)) {
		char *substr = strstr(dirent->d_name, ".co\0"); //end with .co not match .co
		if (substr) {
			snprintf(fullname, 512, "%s/%s", dirpath, dirent->d_name);
			stat(fullname, &statbuf);
			size += statbuf.st_size/8 ;
			// if(size < HASHSIZE * NLF ){	
			strncpy(cofile[n], dirent->d_name, 512);
			fsize[n] = statbuf.st_size/8;
			n++;
			//}else		break;
		}
	}
	closedir(dh);
	printf("loading less than %d .co files:\n",n);
	struct entry *CO = malloc(HASHSIZE*sizeof(*CO));
    //initial CO
    for(unsigned int i=0;i<HASHSIZE;i++){
        CO[i].key=0LLU;
        //CO[i].gc=0; //genome code
    }
	if (!CO) err(errno, "oom");
	double *matrix = malloc(n * n * sizeof(*matrix));
	if (!matrix) err(errno, "oom");
#define MAT(X, Y) (matrix[(X)*n + (Y)])
   // llong gcu = 1LLU << HSPB; // genome count unit
	llong hcl=HASHSIZE, keycount=0 ;// hash collision limited to HASHSIZE
	int j;
	for (j = 0; (j < n) && (keycount < HASHSIZE*NLF); j++) { //n should == HKL
		MAT(j, j) = 0.0;
		//tmp arr for ctx,obj;
		unsigned int *CTX=calloc(j,sizeof(*CTX));
		unsigned int *OBJ=calloc(j,sizeof(*OBJ));
		//unsigned int OBJ[j]={0};
		llong tuple=0LLU;
		llong mobju=1LLU<<SHIFT_MOBJ; //one obj unit; 
		llong gc=(llong)j << SHIFT_GC; // 52 default 
		char fullname[512] = {0};
		snprintf(fullname, 512, "%s/%s", dirpath, cofile[j]);
		FILE *fp = fopen(fullname, "rb");
		if (!fp) {
			err(errno, "can't open file %s", fullname);
		}
        printf("%dth\t",j);
		printf("file:%s\n",cofile[j]);
//hash method: when collision to same ctx(different sp.)--->linear probe,i from 1
//when collision to different ctx--->quadratic probe, i continue increase from last interupt quadratic probe until find empty bucket  
//quadratic proble i from 1 
		while(!feof(fp)&&fread(&tuple,8,1,fp)){
			unsigned int qpi=0,ind;// qudratic probe i //qpi=qw ;
			//ind = h1((tuple&objmask)),reprobe=0;
            for(qpi=0;qpi<hcl;qpi++){
                ind=h((tuple&objmask),qpi);
                if(CO[ind].key==0){
                    CO[ind].key=gc|tuple;
                    break;
                }
                else if ( (CO[ind].key&objmask) == (tuple&objmask) ) {
                    CTX[CO[ind].key>>SHIFT_GC]++;
                    if( CO[ind].key != tuple ) //make sure key is 40bits
                        OBJ[CO[ind].key>>SHIFT_GC]++;
                    CO[ind].key=+mobju; //increase 1,but alway the total copies of the ctx - 1
                    unsigned int moreobjc= (CO[ind].key<<12) >>SHIFT_MOBJ, //make sure 64-SHIFT_GC==12
                    mod=moreobjc%OBJAR,
                    rd=(unsigned int)(moreobjc/OBJAR);
                    unsigned char nobj = (tuple&ctxmask >> CTXLEN)%4; //char type only for HKL<64 otherwise use short
                    objs_t *tmp ;
                    if (mod==1){
                        tmp = malloc(sizeof(objs_t));
                        if(!tmp) err(errno,"oom");
                        tmp->next=CO[ind].new;
                        CO[ind].new=tmp;
                        CO[ind].new->gcobj[0]=(j<<2)+nobj;                       
                    }
                    unsigned char it;   //
                    for(unsigned int i=0;i<mod-1;i++){
                        it=CO[ind].new->gcobj[i];
                        CTX[it>>2]++;
                        if(it%4!=nobj) OBJ[it>>2]++;
                    }
                    tmp=CO[ind].new->next;
                    CO[ind].new->gcobj[mod-1]=(j<<2)+nobj;
                    for(unsigned int blk=0;blk<rd;blk++){
                        for(int i=0;i<OBJAR;i++){
                            it=tmp->gcobj[i];
                            CTX[it>>2]++;
                            if(it%4!=nobj) OBJ[it>>2]++;
                        }
                        tmp=tmp->next;
                    }
                    break;
                }
            }
		}
		for(int i=0;i<j;i++){
	//		if ( (CTX[i] >SCLB)&& ((double)CTX[i]/(fsize[i]+fsize[j]-CTX[i]) > JCLB ) && ((double)OBJ[i] / (double)CTX[i] < ULB) )
				MAT(j, i) = MAT(i, j) = (double)OBJ[i] / (double)CTX[i];
	//		else
	//			MAT(j, i) = MAT(i, j) = (double)URD;
		}
		fclose(fp);
	}
	printf("total %d files loaded\n",j);
        //free CO
        objs_t* tmp;
        for (unsigned int i=0;i<HASHSIZE;i++){
            while(CO[i].new){
                tmp=CO[i].new;
                CO[i].new=tmp->next;
                free(tmp);
            }
        }
		free(CO);

	printf("%i\n", n);
	for (int j = 0; j < n; j++) {
		char *ptr = strchr(cofile[j], '.');
		if (ptr) {
			*ptr = '\0';
		}
		printf("%-9s", cofile[j]);
		for (int k = 0; k < n; k++) {
			printf(" %lf", MAT(j, k));
		}
		printf("\n");
	}
	return 0;
}
