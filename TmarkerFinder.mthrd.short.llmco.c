/* *************************linked list massive-co for>64 genomes***************************************
// parallelable pairwise co-distances computing method using linked list
// highly compact ctx storage, highly immune to clustering problem   
//                  contact: yhg926@gmail.com  
//                   beta version: 2018-3-24
// ********************************************************************************/
#include "comp_rvs.h"
#include "ctxobj.h"
#include "basemap.h"
#include <dirent.h>
#include <err.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#ifdef _OMP
#include <omp.h>
#endif
#define URD 0.5 // unrelated taxa distance: defined distace to URD
#define SUBC 0.25 //subcluster distance set to SUBC
#define ULB 0.09  //unrelated taxa boundary: based on C9,9O1 co-distance 
#define SCLB 1000 // for 3M # share ctx should >150K, for extreme small genome,need impose a JCLB
#define JCLB 0.08  //shared ctx/smaller genome or jaccard index lower boundary
//the optimized hash methods and parametors were chose here based on a test on 358 Salmonella genomes(20180323)
#define HASHSIZE 1000000009 //a prime 400000009
#define NLF 0.6 
#define h1(k) ((k)%(HASHSIZE)) 
#define h2(k) (1+((k)%(HASHSIZE-1)))
#define h(k,i) ((h1(k) + i * h2(k))%HASHSIZE)//test: much faster than h(k,i)=(h1(k)+i(i+1)/2)%HASHSIZE
//#define BTS 2  //bytes store SP.count+obj(2bits)
#define GCB 12 //genome code bits GCB=BTS*8-2
#define HKL 800 //4096 //maximum #sp. allowed in a hashtable,HKL = 2^GCB 
#define OBJAR 32 //16 obj array size
#define SHIFT_GC 52 //high 64-HSPB bits for genome code
#define SHIFT_MOBJ 40 //for count of other genmes match this ctx
#define ALEL_TH 0.01 // allel frequency lower than ALEL_TH can not consider as a real allel ,might be error?  

static llong objmask = OBJMASK;
static llong tupmask =TUPMASK;
//static llong ctxmask = CTXMASK;
typedef struct objs{unsigned short gcobj[OBJAR]; struct objs *next;} objs_t; //gcobj[i]:genomecode:6bits,obj:2bits
struct entry {llong key; objs_t *new;} entry;

int THREADS = 0;

void help(int exit_code)
{
	static const char str[] = "Usage: llmcodist DIR\n";
	fprintf(exit_code == 0 ? stdout : stderr, str);
	exit(exit_code);
}
long long timeInMilliseconds(void) {
    struct timeval tv;

    gettimeofday(&tv,NULL);
    return (((long long)tv.tv_sec)*1000)+(tv.tv_usec/1000);
}
long long timeInMicroseconds(void) {
    struct timeval tv;

    gettimeofday(&tv,NULL);
    return (((long long)tv.tv_sec)*1000000)+(tv.tv_usec);
}
struct timespec diff(struct timespec start, struct timespec end)
{
        struct timespec temp;
        if ((end.tv_nsec - start.tv_nsec) < 0) 
        {
                temp.tv_sec = end.tv_sec - start.tv_sec - 1;
                temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
        } 
        else 
        {
                temp.tv_sec = end.tv_sec - start.tv_sec;
                temp.tv_nsec = end.tv_nsec - start.tv_nsec;
        }
        return temp;
}

int main(int argc, char *argv[])
{
#ifdef _OMP
	THREADS = omp_get_num_procs();
#endif

	int index;
	while ((index = getopt(argc, argv, "h"))) {
		if (index == -1) break;
		if (index == 't'){
#ifdef _OMP
			// parse
			errno = 0;
			unsigned long threads = strtoul(optarg);
			if (errno || end == optarg || *end != '\0') {
				warnx("Expected a number for -t argument, but '%s' was "
					  "given. Ignoring -t argument.",
					  optarg);
				break;
			}

			if (threads > omp_get_num_procs()) {
				warnx(
					"The number of threads to be used, is greater than the "
					"number of available processors; Ignoring -t %lu "
					"argument.",
					threads);
				break;
			}
#else
			warnx("This version of co2dist_p was built without OpenMP and "
				  "thus does not support multi threading. Ignoring -t "
				  "argument.");
			break;
#endif
		} else if (index == 'h') help(0);
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
	unsigned int n = 0;
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
	FILE *taxon;
	taxon=fopen("./taxoncode.out","w");
	for(unsigned int i=0;i<n;i++)
		fprintf(taxon,"%d\t%s\n",i,cofile[i]);

	fclose(taxon);

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
	llong hcl=HASHSIZE, keycount=0 ;// hash collision limited to HASHSIZE
	llong *fStore=malloc(20000000*8);
        unsigned int CTX[4096];//={0};//calloc(j,sizeof(*CTX));
        unsigned int OBJ[4096];//={0};//calloc(j,sizeof(*OBJ));
        llong mobju=1LLU<<SHIFT_MOBJ,gc;FILE *fp;llong start,end;
#pragma omp parallel num_threads(20) //private(gc,fullname,fp) 
{
	for (unsigned short j = 0; j < n; j++) { //n should == HKL

#pragma omp single 
{
		MAT(j, j) = 0.0;
		memset(CTX, 0, sizeof(CTX));		
		memset(OBJ, 0, sizeof(OBJ));
		gc=(llong)j << SHIFT_GC; // 52 default 
		snprintf(fullname, 512, "%s/%s", dirpath, cofile[j]);
		fp = fopen(fullname, "rb");
		if (!fp) {
			err(errno, "can't open file %s", fullname);
		}
		fread(fStore,8,fsize[j],fp);
		start=timeInMilliseconds();
  //     printf("%dth\t",j);
    //    printf("file:%s\n",cofile[j]);		
}
#pragma omp for reduction(+:CTX,OBJ) schedule(guided)  // need gcc version >6.0 for reduction array
		for(unsigned int tupit=0;tupit<fsize[j];tupit++){//tupit<fsize[j]; //replace while loop so can apply #pragma omp parallel for 
			llong tuple=fStore[tupit];
			unsigned int qpi,ind;// private
			for(qpi=0;qpi<hcl;qpi++){
                ind=h((tuple&objmask),qpi);
     	           if(__sync_bool_compare_and_swap(&CO[ind].key,0LLU,gc|tuple)){ //,gc|tuple
                    break;
                }                
                else if ( (CO[ind].key&objmask) == (tuple&objmask) ) {
                    CTX[(int)(CO[ind].key >>SHIFT_GC)]++;
                    if( (CO[ind].key&tupmask) != tuple ) 
                        OBJ[(int)(CO[ind].key>>SHIFT_GC)]++;
					CO[ind].key+=mobju;//fix CO[ind].key=+mobju; //increase 1,but alway the total copies of the ctx - 1
					unsigned int moreobjc= (CO[ind].key<<12) >>52, //make sure 64-SHIFT_GC==12,SHIFT_OBJM+12==52
					rd=(unsigned int)(moreobjc/OBJAR);
					if( moreobjc%OBJAR==0 ) rd--;//rd:extra chunks of objs;
					int mod=moreobjc%OBJAR; //mod should be signed to make sense mod-1 in line 160, otherwise it become UNIT_MAX so loop not end;  
                    unsigned char nobj = (tuple >> CTXLEN)%4; //char type only for HKL<64 otherwise use short
					objs_t *tmp ;
                    if (mod==1){
                        tmp = malloc(sizeof(objs_t));
                        if(!tmp) err(errno,"oom");
                        tmp->next=CO[ind].new;
                        CO[ind].new=tmp;
                    }
                    unsigned short it;   //
		    		int rmod = (OBJAR+mod-1)%OBJAR;// correct mod if mod -1 < 0; mod=0,1,2..15,rmod=15,0,1..14;
                    for(int i=0;i< rmod ;i++){ // i take 0..14
                        it=CO[ind].new->gcobj[i];
                        CTX[it>>2]++;
                        if(it%4!=nobj) 
							OBJ[it>>2]++;
                    }
                    tmp=CO[ind].new->next;
                    CO[ind].new->gcobj[rmod]=(unsigned short)(j<<2)+nobj; // when rd>0,mod==0, mod-1 != (OBJAR+mod-1)%OBJAR,
		    		for(unsigned int blk=0;blk<rd;blk++){
                        for(unsigned int i=0;i<OBJAR;i++){
                            it=tmp->gcobj[i];
                            CTX[it>>2]++;
                            if(it%4!=nobj)
								 OBJ[it>>2]++;
                        }
                        tmp=tmp->next;
                    }
					break;
                }//else if end				
				//else break;
            }//hash a key
		}//hash a file
#pragma omp single nowait
{
		end=timeInMilliseconds();
//		printf("malloc took %llums\n", end-start);
	fclose(fp); //tool 0ms
}
#pragma omp for //took 0ms
		for(unsigned short i=0;i<j;i++){
	//		if ( (CTX[i] >SCLB)&& ((double)CTX[i]/(fsize[i]+fsize[j]-CTX[i]) > JCLB ) && ((double)OBJ[i] / (double)CTX[i] < ULB) )
				MAT(j, i) = MAT(i, j) = (double)OBJ[i] / CTX[i];
	//		else
	  //   		MAT(j, i) = MAT(i, j) = (double)URD;
		}
	}//hash all files
}
        //free CiO
FILE *fout;
fout= fopen("./markerMembership.out", "w");

basemapset(Basemap);
objs_t* tmp;
int nonsingleton=0,ctx_count=0;
unsigned int allel[4],allel_max,allel_min;
char allel_ct;
unsigned short tmparr[512][4]; //keep tmp genome code; this size should > #total genomes(this testdata have 272 genomes); 
for (unsigned int i=0;i<HASHSIZE;i++){
	if(CO[i].key==0)
		continue;
	else{
		ctx_count++;
		if(CO[i].new==NULL)
			continue;
	}
	nonsingleton++;
	memset(allel,0,4*sizeof(unsigned int));
	allel_ct=0,allel_max=0,allel_min=10000;
	for(int k=0;k<TUPLEN;k++)
                printf("%c", Mapbase[(CO[i].key<<(64-2*(TUPLEN-k)))>>62]) ;
        unsigned int obj_count = (CO[i].key<<12) >>52;	
	printf("\t%d\t",obj_count);	
	tmp=CO[i].new; 
	for(int k = obj_count - 1 ; k >=0 ; k--){
		unsigned short it = tmp->gcobj[k%OBJAR];
	//	printf("\t%d|%c",it>>2, Mapbase[it%4]);
		tmparr[allel[it%4]][it%4] = it>>2;
		allel[it%4]++;
		if(k%OBJAR==0){
			CO[i].new=tmp->next;
			free(tmp);
			tmp=CO[i].new;
		}
	}
	for(int k=0;k<4;k++){
		if(allel[k]>0){ 
			printf("%c:%d|",Mapbase[k],allel[k]);
			allel_ct++;
			if (allel[k] > allel_max)
				allel_max = allel[k] ;
			if (allel[k] < allel_min)
                                allel_min = allel[k] ;
		}	
	}
	printf("\t%d\t%f\n",allel_ct, log((double)allel_max/allel_min));
	//print to file:
	if(obj_count>150 && allel_ct==2){		
		for(int k=0;k<4;k++){	
			if(allel[k]>0){
			 	fprintf(fout,"%d\t%010llx", allel[k],(CO[i].key&OBJMASK)|((llong)k<<CTXLEN) ); 
				for(int j=allel[k]-1; j>=0;j-- )
					fprintf(fout,"\t%d",tmparr[j][k]);
				fprintf(fout,"\n");
			}	
		}
	}	
}
fclose(fout);
//printf("ctx_count:%d\tnonsingleton:%d\n",ctx_count,nonsingleton);
free(CO);
	
/*
	for (unsigned short i = 0; i < n; i++) {
		char *ptr = strchr(cofile[i], '.');
		if (ptr) {
			*ptr = '\0';
		}
		printf("%-9s", cofile[i]);
		for (unsigned short k = 0; k < n; k++) {
			printf(" %lf", MAT(i, k));
		}
		printf("\n");
	}
*/
	return 0;
}

