#co-phylog tools set

Please cite:
“Co-phylog: an assembly-free phylogenomic approach for closely related organisms
H Yi, L Jin
Nucleic acids research 41 (7), e75-e75”

======INSTALL AND USAGE========

1.	Compile:

gcc fasta2co.v18.3.c -O6 -Wall -o fasta2co ;

gcc fastq2co.c -O6 -Wall -o fastq2co ;

gcc co2dist2.c -O6 -Wall -o co2dist  ;

gcc readco.c  -Wall -o readco ;

2.	Usage:

fasta2co : Convert fasta file to co file ;
Usage : ./fasta2co <*.fasta> <*.co>

fastq2co: Convert fastq file to co file ;
Usage : ./fasta2co <*.fastq> <*.co>

readco: read binary co file 
Usage : ./readco <string/hex>

co2dist: generated pairwise distance from co files 
Usage: ./co2dist <co file dir> > *.dist

dis_matr_ge.pl: Generated standard distance matrix and orangnism name code from pairwise distances
Usage : perl ./Pl/dis_matr_ge.pl [pairwise distances file] [organism number] >[output file]

nwkrename.pl
	When using 'neighbor' in Phylip package to generated tree file. if you need substitute
	the output code using origin organism name,you could use this program
Usage: perl ./Pl/nwkrename.pl [outtree/outfile] > [output file]

EXAMPLE:

	tar xzvf example_brucella.tar.gz
	cp example_brucella/* algorithm/Genomes/
    cd algorithm
	
	1) Build cofile files using complete genomes

		for i in `ls Genomes`;
			do ./C/fasta2co ./Genomes/$i ./CO-index/$i.co ;
		done;

    You can also bulid co file using fastq NGS data by:
	   for i in `ls Genomes`;
	       do ./C/fasta2co ./Genomes/$i ./CO-index/$i.co ;
		done;

	2) Generate pairwise distances using index files
	./C/co2dist CO-index >dist

	3) format the pairwise distances to standard distance matrix
	        perl ./Pl/dis_matr_ge.pl dist 15 >dis_matrix

	4) If the 'neighbor' contained in Phylip package has been installed then
			neighbor <<END
						dis_matrix
						Y
   			        END

	5) If want to replace the namecode by organisms' name								    
		perl ./Pl/nwkrename.pl outtree
