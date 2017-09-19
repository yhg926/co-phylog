#include <stdio.h>
#include <stdlib.h>
#include "ctxobj.h"
#include "basemap.h" 
int main(int argc, char* argv[]){
	if(argc<2)
		printf("USAGE:./a.out <*.co> <string/hex>");
	FILE *fp;
	if((fp=fopen(argv[1],"rb"))==NULL){
		printf("can't open %s", argv[1]);
		exit(1);
	}
	basemapset(Basemap);
	llong tuple=0LLU;
	while(!feof(fp)){
		fread(&tuple,8,1,fp);
		if(argv[2]){
			int i;
			for(i=0;i<TUPLEN;i++)
				printf("%c", Mapbase[(tuple<<(64-2*(TUPLEN-i)))>>62]) ;
			printf("\n");
		//	exit(1);		
		}
		else
			printf("%010llx\n",tuple);	
	}
	fclose(fp);
	return 1;

}
