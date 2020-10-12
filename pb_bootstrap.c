/****************************************************************************
 ****************************************************************************
 *  Copyright (C) 2013-2023                                                 *
 *  Author: Huiguang Yi (yhg926@gmail.com)            	                    *
 *  All rights reserved.                   		                    *
 ****************************************************************************
 ****************************************************************************/

/******bootstrap.c******/
/***limit to 255 co files in a co dir***/
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include "ctxobj.h"  
#define LF 0.5
#define COFN 255 
#define LMARK 0x8000000000000000LLU  
#define h1(k) ((k)%(hashsize)) 
#define h2(k) (1+((k)%(hashsize-1)))
#define h(k,i) ((h1(k) + i * h2(k))%hashsize)

static llong objmask = OBJMASK;
static llong ctxmask = TUPMASK^OBJMASK;
typedef struct copos {llong co; unsigned int pos;} copos;
typedef struct entry {struct entry* next;unsigned int site; unsigned char fnum; unsigned char nt;} entry;
typedef struct primentry {struct entry* next;unsigned int site; unsigned char nt;} primentry;

int main(int argc, char*argv[])
{
    if(argc!=3){
        printf("USAGE: ./bootsrap <co dir> <resample times N>\n");
        exit(1);
    }
    DIR *dh;
    if((dh=opendir(argv[1]))==NULL){
        printf("can't open %s\n",argv[1]);
        exit(1);
    }

    struct stat statbuf;
    struct dirent * dirent; 
	unsigned int hashsize, fcount=0; 
    char fullname[512];
	FILE *fp;	
	llong *filemap;
	struct copos *coarray[COFN];
	struct primentry *entryarray[COFN];
	unsigned int arz[COFN]; 
	unsigned int ulz[COFN]; 
	unsigned int haz[COFN]; 
    char *filename[COFN];	
    while ((dirent=readdir(dh))!=NULL){
        char *substr="";
        if( ((substr=strstr((*dirent).d_name, ".co"))==NULL) || (strcmp(substr, ".co")!=0) )
			continue;
		filename[fcount]=calloc(256,sizeof(char *));
		strcpy(filename[fcount],(*dirent).d_name);
		memset(fullname,0,512);
        strcat(strcat( strcat(fullname,argv[1]),"/"),(*dirent).d_name);
        if((fp=fopen(fullname,"rb"))==NULL){            
            printf("can't open file %s",fullname);
            exit(1);
        }
		stat(fullname, &statbuf);
		if((filemap = malloc(statbuf.st_size))==NULL){
			 printf("filemap allocate memory failed!\n");
			 exit(1);
		};
		fread(filemap,statbuf.st_size,1,fp);
		fclose(fp);		

		unsigned int i,j,offset,ind,linknum=0, nco = arz[fcount] = statbuf.st_size/sizeof(llong) ;
		llong ctx;
		if(fcount >0 ){
			for (i=0;i < nco;i++){
				ctx=filemap[i]&objmask;
				for(j=0;j<fcount;j++){
					hashsize = haz[j];
					for(offset=0;offset < hashsize ;offset++){
						ind=h(ctx,offset);
						if((coarray[j][ind].co & objmask) == ctx){
						   struct entry *elm;
               			   if((elm = calloc(1,sizeof(struct entry *)))==NULL){
                     			printf("elm allocate memory failed!\n");
                     			exit(1);
                			}; 
                			(*elm).fnum=fcount;
                			(*elm).site=i;              
                			(*elm).nt=(filemap[i]&ctxmask)>>CTXLEN;	
							(*elm).next = entryarray[j][coarray[j][ind].pos].next;
							entryarray[j][coarray[j][ind].pos].next = elm;
							filemap[i]|=LMARK; 
							linknum++;  	
							goto LOOP;
						}
						else if(coarray[j][ind].co == 0)
							break;
						
					}
				}
				LOOP: ;
			}
		}
		unsigned int unlinknum = ulz[fcount] = nco - linknum ;
		if((entryarray[fcount]  = calloc(unlinknum,sizeof(struct primentry)))==NULL){
            printf("primarray %d allocate memory failed!\n",fcount);
            exit(1);
        }; 
		hashsize = haz[fcount] = (unsigned int)(unlinknum/LF);
		if((coarray[fcount] = calloc(hashsize,sizeof(struct copos)))==NULL){
            printf("coarray %d allocate memory failed!\n",fcount);
            exit(1);
        };
		j=0;
		for (i=0;i < nco ;i++){
			if (filemap[i] < LMARK) {
            	entryarray[fcount][j].nt=(filemap[i]&ctxmask)>>CTXLEN;
               	entryarray[fcount][j].next=NULL;
				entryarray[fcount][j].site = i	;
              	ctx=filemap[i]&objmask;
                for(offset=0;offset < hashsize;offset++){
                    ind=h(ctx,offset);
                    if(coarray[fcount][ind].co == 0){                                                                                           						coarray[fcount][ind].co=filemap[i]; 
                       	coarray[fcount][ind].pos = j;
                      	break;
                    }
                }
				j++; 
			}
        } 
		if (j!=unlinknum){  
            printf("no LMARK non equal to unlinknum\n");
			exit(1);	 
		}                                                                                                                                       
		free(filemap);
		fcount++;           
	}
	closedir(dh);
	int i;
	for	(i=0;i< fcount;i++)
		free(coarray[i]);
	printf("/************linking finished, bootstraping begins***********************/\n");
	unsigned int context[fcount][fcount];
	unsigned int dobj[fcount][fcount];
	unsigned int *bootstrap[fcount];
	for (i=0;i< fcount;i++){
		if(arz[i] > RAND_MAX){
			printf("arz %d >RAND_MAX!\n",i);
			exit(1);
		}
		if((bootstrap[i] = calloc(arz[i],sizeof(unsigned int)))==NULL){
			printf("bootstrap %d allocate memory failed!\n",i);
			exit(1);
		};
	}
	unsigned int j,k,pow,run = atoi(argv[2]);
	struct entry *node,*node2;
	srand( (unsigned int)time(NULL)) ;
	for (i=0;i<run;i++){
		memset(context,0,fcount*fcount*sizeof(unsigned int));
		memset(dobj,0,fcount*fcount*sizeof(unsigned int));
		for(j=0;j<fcount;j++){
			for (k = 0; k < arz[j]; k++){
				bootstrap[j][(unsigned int)rand() % arz[j]]++ ;
			}
		}
		printf("%dth run bootstraping finished!\t",i+1);
		for (j=0;j<fcount-1;j++){ 
			for(k=0;k< ulz[j];k++){
				for(node=entryarray[j][k].next; node!=NULL;node=(*node).next){
					if(bootstrap[j][entryarray[j][k].site] > bootstrap[(*node).fnum][(*node).site])
						pow=bootstrap[(*node).fnum][(*node).site];
					else 
						pow=bootstrap[j][entryarray[j][k].site];
					context[j][(*node).fnum]+=pow;
					context[(*node).fnum][j]+=pow;
					if (entryarray[j][k].nt != (*node).nt){
						dobj[j][(*node).fnum]+=pow;
						dobj[(*node).fnum][j]+=pow;
					}

					for(node2=(*node).next;node2!=NULL;node2=(*node2).next){
						if(bootstrap[(*node).fnum][(*node).site] > bootstrap[(*node2).fnum][(*node2).site])
                        	pow=bootstrap[(*node2).fnum][(*node2).site];
                    	else 
							pow=bootstrap[(*node).fnum][(*node).site];
                    	context[(*node).fnum][(*node2).fnum]+=pow;
						context[(*node2).fnum][(*node).fnum]+=pow;
						if( (*node).nt != (*node2).nt ){	
							dobj[(*node).fnum][(*node2).fnum]+=pow;
							dobj[(*node2).fnum][(*node).fnum]+=pow;
						}
					}
					
				}

			}	
		}
		char str[50];
		memset(str,0,50);
		sprintf(str,"./bootstap.%u.dist",i);
		if((fp=fopen(str,"w"))==NULL){            
            printf("can't open file %s",str);
            exit(1);
        }
		for (j=0;j<fcount;j++){
			for(k=j+1;k<fcount;k++){
				fprintf(fp,"%s\t%s\t%f\n",filename[j],filename[k],(double)dobj[j][k]/context[j][k]) ;
			}
			memset(bootstrap[j],0,arz[j]*sizeof(unsigned int));
		}
		printf("%dth run pairwise disances computed.\n",i+1);
		fclose(fp);
	}
	


	return 1;
}




































