//the paralelism is contributed by @kloetzl
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

#ifdef _OMP
#include <omp.h>
#endif

#define HASHSIZE 20000081 // 40000161 //20000081*2 //default hashsize
#define h1(k) ((k) % (HASHSIZE))
#define h2(k) (1 + ((k) % (HASHSIZE - 1)))
#define h(k, i) ((h1(k) + i * h2(k)) % HASHSIZE)

static int comp_bittl = 64 - BITTL;
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

	int hashsize = HASHSIZE;
	struct dirent *dirent;
	int count = 0;

	while ((dirent = readdir(dh)) != NULL) {
		char *substr = strstr(dirent->d_name, ".co");
		if (substr) {
			count++;
		}
	}

	char cofile[count][512];
	memset(cofile, 0, 512 * count);
	rewinddir(dh);

	int n = 0;
	while ((dirent = readdir(dh)) != NULL) {
		char *substr = strstr(dirent->d_name, ".co");
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

		llong tuple = 0LLU;
		while (!feof(fp) && fread(&tuple, 8, 1, fp)) {
			for (int offset = 0; offset < hashsize; offset++) {
				int ind = h((tuple & objmask), offset);
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

			int cxt = 0, obj = 0;
			while (!feof(fp) && fread(&tuple, 8, 1, fp)) {
				for (int offset = 0; offset < hashsize; offset++) {
					int ind = h((tuple & objmask), offset);
					if ((CO[ind] & objmask) == (tuple & objmask)) {
						cxt++;
						if (CO[ind] != tuple) obj++;
						break;
					} else if (CO[ind] == 0) {
						llong crvstuple = crvs64bits(tuple) >> comp_bittl;
						for (offset = 0; offset < hashsize; offset++) {
							ind = h((crvstuple & objmask), offset);
							if (CO[ind] == 0)
								break;
							else if ((CO[ind] & objmask) ==
									 (crvstuple & objmask)) {
								cxt++;
								if (CO[ind] != crvstuple) obj++;
								break;
							}
						}
						break;
					}
				}
			}
			fclose(fp);
			MAT(j, k) = MAT(k, j) = (double)obj / (double)cxt;
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
