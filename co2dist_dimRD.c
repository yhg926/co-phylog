#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include "ctxobj.h"
#include "comp_rvs.h"
#include "dim_rdc.h"
#define HASHSIZE 4000081 //40000161 //20000081*2 //default hashsize
#define h1(k) ((k)%(HASHSIZE)) 
#define h2(k) (1+((k)%(HASHSIZE-1)))
#define h(k,i) ((h1(k) + i * h2(k))%HASHSIZE)

//static int comp_bittl = 64-BITTL ;
//static llong objmask = OBJMASK;
time_t ti,newti;
//int collision,coll;
int main(int argc, char*argv[])
{	
//	argv[1]="cofile";
	if(argc!=2){
		printf("USAGE: ./co2dist <co dir>\n");
		exit(1);
	}
	DIR *dh;
	if((dh=opendir(argv[1]))==NULL){
		printf("can't open %s\n",argv[1]);
		exit(1);
	}

	struct stat statbuf;
	int hashsize = HASHSIZE,size=0;//statbuf.st_size/LF;
	struct dirent * dirent; 
	int count=0;
    char fullname[512];
	while ((dirent=readdir(dh))!=NULL){
	    memset(fullname,0,512);
		char *substr="";
		 if( ((substr=strstr((*dirent).d_name, ".co"))!=NULL) && (strcmp(substr, ".co")==0) ){
			 count++;
			 strcat(strcat( strcat(fullname,argv[1]),"/"),(*dirent).d_name);
			 stat(fullname, &statbuf);
			 if(statbuf.st_size > size)
				size= statbuf.st_size ;
			 //printf("%d\t%f\t%d\n",size,LF,(int)(size/8/LF));
		}
		 
	}
	//get hash size
	if((int)(size/8/LF)*2 > HASHSIZE)
		hashsize=(int)(size/8/LF)*2;
	printf("%d\n",hashsize);
	char cofile[count][512];
	memset(cofile,0,512*count);	
	rewinddir(dh);
	int i,n=0;
	while ((dirent=readdir(dh))!=NULL){
        char *substr="";
        if( ((substr=strstr((*dirent).d_name, ".co"))!=NULL) && (strcmp(substr, ".co")==0) ){
			for(i=0;(cofile[n][i]=((*dirent).d_name)[i])!='\0';i++) ;
            n++;
		}
    }

	closedir(dh);
	llong tuple=0LLU;
	llong *CO;
	CO=malloc(hashsize*sizeof(llong));
	FILE *fp;
	int j,k,ind,offset;
	for(j=0;j< n-1;j++){
		memset(CO,0,sizeof(llong)*hashsize);
		memset(fullname,0,512);
		strcat(strcat( strcat(fullname,argv[1]),"/"),cofile[j]);
		if((fp=fopen(fullname,"rb"))==NULL){			
            printf("can't open file %s",fullname);
            exit(1);
        }
		//printf("%s\t",cofile[j]);
		while(!feof(fp)&&fread(&tuple,8,1,fp)){
			for(offset=0;offset<hashsize;offset++){
				ind=h((tuple>>2),offset);
				if(CO[ind]==0){
					CO[ind]=tuple;
					break;
				}
			}
		}
		fclose(fp);
		for(k=j+1;k<n;k++){
			//time(&ti);
			memset(fullname,0,512);
			strcat(strcat( strcat(fullname,argv[1]),"/"),cofile[k]);
			if((fp=fopen(fullname,"rb"))==NULL){
				printf("can't open file %s",fullname);
				exit(1);
			}
			//printf("%s\t",cofile[k]);
			int cxt=0,obj=0;
			while(!feof(fp)&&fread(&tuple,8,1,fp)){
               // printf("%llu\n",tuple);
				for(offset=0;offset<hashsize;offset++){
					ind=h((tuple>>2),offset);
            //        printf("%llx\t%llx\n",CO[ind],tuple);
					if((CO[ind]>>2)==(tuple>>2)){
						cxt++; 
                       // printf("%llx\t%llx\n",CO[ind],tuple);
						if(CO[ind]!=tuple)
							obj++;
						break;
					}
					else if(CO[ind]==0)
						break;					
				}
			}
          //  exit(1);
			fclose(fp);
			//printf("%f\t%d\t%d\t",(double)obj/(double)cxt,obj,cxt);
			printf("%s\t%s\t%f\t%d\t%d\t\n",cofile[j],cofile[k],(double)obj/(double)cxt,obj,cxt);
			//time(&newti);
		//	printf("%ld\n",newti-ti);
        }
	}
return 1;
}

















