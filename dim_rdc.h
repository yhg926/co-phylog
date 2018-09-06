#ifndef DIM_RDC_H
#define DIM_RDC_H
#define DO 3 //dimension order: one side flanking nt length taken for reduction
#define PD 4096 //(pow(16,DO)) // core k-mer dimension before reduction
#define DRR 16 //dimension reduction rate
#define DRbits 8 //bits for after dimention reduction
#define PSTDR 256 // must be 1<<DRbits, (int)(PD/DRR) //dimension post reduction
#endif
