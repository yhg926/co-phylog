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
//set remoted related taxa as unrelated cluster: defined by co-distance >ULB, or #common ctx < SCLB
//for which a unified distance UNRD is set 
//co-dist < ULB && (#shared ctx > SCOL)&& jaccard index >JCLB or co-distance = URD
#define URD 0.5 // unrelated taxa distance: defined distace to URD
#define SUBC 0.25 //subcluster distance set to SUBC
#define ULB 0.09  //unrelated taxa boundary: based on C9,9O1 co-distance 
#define SCLB 10000 // for 3M # share ctx should >150K, for extreme small genome,need impose a JCLB
#define JCLB 0.08  //shared ctx/smaller genome or jaccard index lower boundary
//========/these boundary are empirical only===============
/*#define GEB 0.05 // genus boundary 
#define SPB 0.01 // ?species boundary 
#define POPB 0.002 //within POPULATION distance; co-distance < POPB defined as the same population
*/

#ifdef _OMP
#include <omp.h>
#endif

#define LF 0.5
#define HSL 40000000 // limit for largest hash size
#define hs(fsize) ((llong)((fsize)/8/LF) ) // get hashsize from filesize
#define h1(k) ((k)%(hashsize)) 
#define h2(k) (1+((k)%(hashsize-1)))
#define h(k,i) ((h1(k) + i * h2(k))%hashsize)

static llong objmask = OBJMASK;

int THREADS = 0;

void help(int exit_code)
{
	static const char str[] = "Usage: co2dist_p [-t THREADS] DIR\n";
	fprintf(exit_code == 0 ? stdout : stderr, str);
	exit(exit_code);
}

int main(int argc, char *argv[])
{
#ifdef _OMP
	THREADS = omp_get_num_procs();
#endif

	int index;
	while ((index = getopt(argc, argv, "t:h"))) {
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

	//	argv[1]="cofile";
	if (argc != 1) {
		help(1);
	}

	DIR *dh;
	if ((dh = opendir(dirpath)) == NULL) {
		err(errno, "can't open %s\n", dirpath);
	}
	char fullname[512] = {0};
	struct stat statbuf;
	unsigned int hashsize,size=0;
	struct dirent *dirent;
	int count = 0;

	while ((dirent = readdir(dh)) != NULL) {
		char *substr = strstr(dirent->d_name, ".co\0");
		if (substr) {
			count++;
			snprintf(fullname, 512, "%s/%s", dirpath, dirent->d_name);
			stat(fullname, &statbuf);
			 if(statbuf.st_size > size)
				size = statbuf.st_size ;
		}
	}
	hashsize = hs(size);
	if (hashsize > HSL)
		err(errno, "hash size larger than HSL!");
	char cofile[count][512];
	memset(cofile, 0, 512 * count);
	rewinddir(dh);

	int n = 0;
	while ((dirent = readdir(dh)) != NULL) {
		char *substr = strstr(dirent->d_name, ".co\0");
		if (!substr) continue;
		strncpy(cofile[n], dirent->d_name, 512);
		n++;
	}

	closedir(dh);

	double *matrix = malloc(n * n * sizeof(*matrix));
	if (!matrix) err(errno, "oom");
#define MAT(X, Y) (matrix[(X)*n + (Y)])

#pragma omp parallel for num_threads(THREADS)
	for (int j = 0; j < n - 1; j++) {

		llong *CO = malloc(hashsize * sizeof(*CO));
		if (!CO) err(errno, "oom");
		memset(CO, 0, sizeof(*CO) * hashsize);

		char fullname[512] = {0};
		snprintf(fullname, 512, "%s/%s", dirpath, cofile[j]);
		FILE *fp = fopen(fullname, "rb");
		if (!fp) {
			err(errno, "can't open file %s", fullname);
		}
		struct stat statbuf;
		stat(fullname, &statbuf);
		//kmer count
		llong kmcr = statbuf.st_size/8 ; 
		llong tuple = 0LLU;
		while (!feof(fp) && fread(&tuple, 8, 1, fp)) {
			for (unsigned int offset = 0; offset < hashsize; offset++) {
				unsigned int ind = h((tuple & objmask), offset);
				if (CO[ind] == 0) {
					CO[ind] = tuple;
					break;
				}
			}
		}
		fclose(fp);

		MAT(j, j) = 0.0;

		for (int k = j + 1; k < n; k++) {
			char fullname[512];
			snprintf(fullname, 512, "%s/%s", dirpath, cofile[k]);
			fp = fopen(fullname, "rb");

			if (!fp) {
				err(errno, "can't open file %s", fullname);
			}
			struct stat statbuf;
			stat(fullname, &statbuf);
			//kmer count
			llong kmcq = statbuf.st_size/8 ; 
			//llong kmc = kmcr < kmcq?kmcr:kmcq;
			unsigned int ctx = 0, obj = 0;
			while(!feof(fp)&&fread(&tuple,8,1,fp)){
				for(unsigned int offset=0;offset<hashsize;offset++){
					unsigned int ind=h((tuple&objmask),offset);
					if((CO[ind]&objmask)==(tuple&objmask)){
						ctx++; 
						if(CO[ind]!=tuple)
							obj++;
						break;
					}
					else if(CO[ind]==0)
						break;
				}
			}
			fclose(fp);
			if ( (ctx >SCLB)&& (ctx/(kmcr+kmcq-ctx) > JCLB ) && ((double)obj / (double)ctx < ULB) )
				MAT(j, k) = MAT(k, j) = (double)obj / (double)ctx;
			else
				MAT(j, k) = MAT(k, j) = (double)URD;

		}

		free(CO);
	}

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
