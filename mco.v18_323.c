/* *************************massive-co.v18_323***************************************
// faster pairwise co-distances computing method for >300 genome in one core 
// (parallism is support as well) by create a genomes-pooled hash tables 
// (need 32G memory per hashtable) and using new hash method hybriding linear probing 
// and quadratic probing
//                  contact: yhg926@gmail.com  
//                   beta version: 2018-3-23
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
#define HASHSIZE 4000000007//test: use a prime(combine with the chosed hash fun.) is faster than use power of 2 of similar size,e.g. 2^32 
#define NLF 0.4 //test: when load factor NLF >0.4,there is no gain of efficiency any more(compare to running smaller batches serially)  
#define h1(k) ((k)%(HASHSIZE)) 
#define h2(k) (1+((k)%(HASHSIZE-1)))
#define h(k,i) ((h1(k) + i * h2(k))%HASHSIZE)//test: much faster than h(k,i)=(h1(k)+i(i+1)/2)%HASHSIZE
#define HKL 2048 //maximum #sp. allowed in a hashtable 

static llong objmask = OBJMASK;
static llong tupmask =TUPMASK;

void help(int exit_code)
{
	static const char str[] = "Usage: mcodist DIR\n";
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
	while ((dirent = readdir(dh)) != NULL) {
		char *substr = strstr(dirent->d_name, ".co\0"); //end with .co not match .co
		if (substr) {
			snprintf(fullname, 512, "%s/%s", dirpath, dirent->d_name);
			stat(fullname, &statbuf);
			size += statbuf.st_size/8 ;
			 if(size < HASHSIZE * NLF ){	
				strncpy(cofile[n], dirent->d_name, 512);
				fsize[n] = statbuf.st_size/8;
				n++;
			}else	
				break;
		}
	}
	closedir(dh);
	printf("loading %d .co files:\n",n);
	llong *CO ;
	CO = calloc(HASHSIZE, sizeof(*CO));
	if (!CO) err(errno, "oom");
	double *matrix = malloc(n * n * sizeof(*matrix));
	if (!matrix) err(errno, "oom");
#define MAT(X, Y) (matrix[(X)*n + (Y)])
#define HSPB 52 //high 64-HSPB bits for species code
	llong hcl=HASHSIZE ;// hash collision limited to HASHSIZE
	for (int j = 0; j < n ; j++) {
		MAT(j, j) = 0.0;
		//tmp arr for ctx,obj;
		unsigned int *CTX=calloc(j,sizeof(*CTX));
		unsigned int *OBJ=calloc(j,sizeof(*OBJ));
			//unsigned int OBJ[j]={0};
		llong tuple=0LLU;
		llong spc=(llong)j << HSPB; // 52 default 
		char fullname[512] = {0};
		snprintf(fullname, 512, "%s/%s", dirpath, cofile[j]);
		FILE *fp = fopen(fullname, "rb");
		if (!fp) {
			err(errno, "can't open file %s", fullname);
		}
		printf("%s\n",cofile[j]);
//hash method: when collision to same ctx(different sp.)--->linear probe,i from 1
//when collision to different ctx--->quadratic probe, i continue increase from last interupt quadratic probe until find empty bucket  
//quadratic proble i from 1 
		while(!feof(fp)&&fread(&tuple,8,1,fp)){
			unsigned int qpi=0,// qudratic probe i //qpi=qw
			ind = h1((tuple&objmask)),reprobe=0;//coltype:last collission type 0,1,2
			while(CO[ind]!=0){
				if((CO[ind]&objmask) != (tuple&objmask)){
					qpi++;
					ind=h((tuple&objmask),qpi) ;//hs(tuple&objmask,qpi); //quadratic probe
				}
				else{ //CO[ind]&obmask) == (tuple&objmask)
					CTX[CO[ind]>>HSPB]++;
					if((CO[ind]&tupmask) != tuple)
						OBJ[CO[ind]>>HSPB]++;
					ind++;
					ind=(ind)%HASHSIZE ;// //linear probe
				}
				reprobe++;
				if (reprobe > hcl)
					err(errno,"#reprobe > hcl");
			}
			CO[ind]= spc|tuple ;
		}
		for(int i=0;i<j;i++){
	//		if ( (CTX[i] >SCLB)&& ((double)CTX[i]/(fsize[i]+fsize[j]-CTX[i]) > JCLB ) && ((double)OBJ[i] / (double)CTX[i] < ULB) )
				MAT(j, i) = MAT(i, j) = (double)OBJ[i] / (double)CTX[i];
	//		else
	//			MAT(j, i) = MAT(i, j) = (double)URD;
		}
		fclose(fp);
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
