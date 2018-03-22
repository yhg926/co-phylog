//fast compute pairwise distance for >500 genome in one core by create a genomes pool hash table   
//need 32G memory for one hashtable
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
#define SCLB 10000 // for 3M # share ctx should >150K, for extreme small genome,need impose a JCLB
#define JCLB 0.08  //shared ctx/smaller genome or jaccard index lower boundary

#define HASHSIZE (BIT1MASK<<30) //(BIT1MASK<<32) // limit for largest hash size
//hash fuc for collision of different ctx
#define HKN 200 // HKN is the number of sp. in the hash table
#define HKL 2048 //up limit of species in a hashtable(2^11)
#define h1(k) ((k)%(HASHSIZE)) 
#define h2(i) (i*(i+1)/2) // i > HKN^(1/2) ,i=21
#define hs(k,i) ((h1(k) + h2(i)) % HASHSIZE)
//hash fuc for collision of same ctx in different sp.(linear proble)
#define hg(k,i) ((h1(k) + i) % HASHSIZE) // i start for 0
/*
		llong tuple = 0LLU;
		while (!feof(fp) && fread(&tuple, 8, 1, fp)) {
			for (unsigned int offset = 0; offset < HASHSIZE; offset++) {
				unsigned int ind = hs((tuple & objmask), offset);
				if (CO[ind] == 0) {
					CO[ind] = tuple;
					break;
				}
			}
		}
		fclose(fp);

*/


static llong objmask = OBJMASK;
static llong tupmask =TUPMASK;

void help(int exit_code)
{
	static const char str[] = "Usage: 4gmco2dist DIR\n";
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

	//	argv[1]="cofile";
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
	int pw=(int)pow(HKN,0.5);
	//printf("pw=%d\n",pw);
	while ((dirent = readdir(dh)) != NULL) {
		char *substr = strstr(dirent->d_name, ".co\0"); //end with .co not match .co
		if (substr) {
			snprintf(fullname, 512, "%s/%s", dirpath, dirent->d_name);
			stat(fullname, &statbuf);
			size += statbuf.st_size ;
			 if(size < HASHSIZE * LF ){	
				strncpy(cofile[n], dirent->d_name, 512);
				fsize[n] = statbuf.st_size/8;
				n++;
			}else	
				break;
		}
	}
	closedir(dh);

	llong *CO ;
	CO = calloc(HASHSIZE, sizeof(*CO));
	if (!CO) err(errno, "oom");
	//double *jcM= malloc(n * n * sizeof(*jcM));
	//if (!jcM) err(errno, "oom");
	double *matrix = malloc(n * n * sizeof(*matrix));
	if (!matrix) err(errno, "oom");
	//unsigned int *ctxM=malloc(n * n * sizeof(*ctxM));
	//if (!ctxM) err(errno, "oom");
	//unsigned int *objM=malloc(n * n * sizeof(*objM));
	//if (!objM) err(errno, "oom");
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
		llong spc=(llong)j << HSPB;
		char fullname[512] = {0};
		snprintf(fullname, 512, "%s/%s", dirpath, cofile[j]);
		FILE *fp = fopen(fullname, "rb");
		if (!fp) {
			err(errno, "can't open file %s", fullname);
		}
//hash method: when collision to same ctx(different sp.)--->linear probe,i from 1
//when collision to different ctx--->quadratic probe, i continue increase from last interupt quadratic probe until find empty bucket  
//quadratic proble i from 1 or pw?
		while(!feof(fp)&&fread(&tuple,8,1,fp)){
			unsigned int qpi=0,// qudratic probe i //qpi=qw
			ind = h1((tuple&objmask)),reprobe=0;//coltype:last collission type 0,1,2
			//0:no collision,1:different ctx(use quadratic probe),2:same ctx,different sp.(use linear probe).
			while(CO[ind]!=0){
				if((CO[ind]&objmask) != (tuple&objmask)){
					qpi++;
					ind=(ind+h2(qpi))%HASHSIZE ;//hs(tuple&objmask,qpi); //quadratic probe
				}
				else{ //CO[ind]&obmask) == (tuple&objmask)
					CTX[CO[ind]>>HSPB]++;
					if((CO[ind]&tupmask) != tuple)
						OBJ[CO[ind]>>HSPB]++;

					ind++;
					ind=(ind)%HASHSIZE ;// hg(tuple&objmask,offset); //linear probe
				}
				reprobe++;
				if (reprobe > hcl)
					err(errno,"#reprobe > hcl");
			}
			CO[ind]= spc|tuple ;

		}
		for(int i=0;i<j;i++){
			if ( (CTX[i] >SCLB)&& (CTX[i]/(fsize[i]+fsize[j]-CTX[i]) > JCLB ) && ((double)OBJ[i] / (double)CTX[i] < ULB) )
				MAT(j, i) = MAT(i, j) = (double)OBJ[i] / (double)CTX[i];
			else
				MAT(j, i) = MAT(i, j) = (double)URD;
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
