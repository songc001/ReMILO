### LATEST NEWS
The ReMILO manuscript is submitted to ISMB'17! 

### Overview
ReMILO is software that detects misassembly errors in the contigs of a species  using  the corresponding short reads and the reference genome of a closely related species, as well as the corresponding long reads.

### Copy right
ReMILO is under the [Artistic License 2.0](http://opensource.org/licenses/Artistic-2.0).

### Short manual
1. System requirements

   ReMILO is suitable for 32-bit or 64-bit machines with Linux operating systems. At least 4GB of system memory is recommended for larger data sets.

2. Installation

   Aligners [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [BWA-MEM](http://bio-bwa.sourceforge.net/) are required to run ReMILO.  
   * The source files in 'src' folder can be compiled to generate a 'bin' folder by running Makefile: `make all`. 
   * Put Bowtie2, BWA-MEM and the 'bin' folder to $PATH: `export PATH=PATH2Bowtie2:$PATH`, `export PATH=PATH2BWA-MEM:$PATH` and `export PATH=PATH2bin:$PATH`, respectively.  


3. Inputs
   * Contigs of a species assembled by any genome assembler (ABySS, ALLPATHS-LG, SOAPdenovo2, Velvet, etc.).
   * Reference genome from a closely related species.
   * Corresponding paired-end short reads in FASTA format.
   * Corresponding long reads in FASTA format (optional).

4. Using ReMILO

   ```
   runReMILO.py contigs.fa reference_genome.fa short_reads1.fa short_reads2.fa [-options | -options]
   ```

   Options (default value):  
   -i/-insert n (500)  
   Insert length of short reads.  
   -k/-kmer n (19)  
   Size of k bases.  
   -d/-distance n (85)  
   Minimum alignment distance of adjacent contig positions to detect a misassembly error (using either reference genome or long reads).  
   -l/-longread long_reads.fa (yes)  
   Corresponding long reads.  
   -c/-coverage n (1)  
   Minimum number of long reads required to detect a misassembly error.

5. Outputs
   * Misassembly locations. Each line of the file records the misassembly locations of one contig and has the following format : `initial contig name : offset of first misassembly error; offset of second misassembly error ; ...`.
   * Split contigs at misassembly locations. Specification of each split contig has the following format: `initial contig name : split contig ID`, where the split contig ID starts from 0.
   * Remaining contigs without detected misassembly locations. 
