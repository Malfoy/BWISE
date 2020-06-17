BWISE
=====

de Bruijn Workflow using Integral information of Short paired-End reads

Work in progress -- so far, this assembler works well  with \> 100X of
2\*250 bp Illumina reads (recommended fragment library size \~500-600 bp).


[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

[![Build Status](https://travis-ci.org/Malfoy/BWISE.svg?branch=master)](https://travis-ci.org/Malfoy/BWISE)



INSTALLATION from git
------------
(requires GCC\>=4.9, GIT, MAKE and CMAKE3).

`git clone https://github.com/Malfoy/BWISE --depth 1  `

`cd BWISE`

`./install_git.sh` 

will download and compile the latest version of the needeed programs 
Possible compilation options:

\-f absolute path of folder to put the binaries

\-t use multiple threads for compilation (default is 8)

You can test your installation as follows:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ bash
cd data
./test_install.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


INSTALLATION from release file (experimental)
------------
(requires GCC\>=4.9, GIT, MAKE and CMAKE3).

Get the latest release:

`wget https://github.com/Malfoy/BWISE/releases/download/V0.1/Bwise01.tgz`

`tar -xvf Bwise01.tgz`

`cd BWISE`


`./install_src.sh` 

will compile the programs of the downloaded release
Possible compilation options:

\-f absolute path of folder to put the binaries

\-t use multiple threads for compilation (default is 8)



UPDATE
------------

`./update.sh` will update the components of Bwise

Possible compilation options:

\-f absolute path of folder to put the binaries

\-t use multiple threads for compilation (default is 8)

You can test your installation as follows:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ bash
cd data
./test.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




RUN
---

`python3 Bwise.py -x examplePairedReads.fa -o workingFolder`

Options and default values:

-h, --help            show this help message and exit

  -x PAIRED_READFILE    input fasta/fastq (fa, fa.gz, fq or fq.gz) paired-end read files. Only one input file of this type can be used: if you have for instance two paired-end librairies (e.g. libA and libB), first interleave the forward and reverse reads for each library. This may be done using the provided script `src/two_fastq_to_interleaved_fasta.py`:

```bash
python src/two_fastq_to_interleaved_fasta.py libA_1.fq libA_2.fq > interleavedlibA.fa
```

and 

```bash
python src/two_fastq_to_interleaved_fasta.py libB_1.fq libB_2.fq > interleavedlibB.fa
```

and then concatenate the two resulting files  in order to generate the PAIRED_READFILE input fasta): 

```cat interleavedlibA.fa interleavedlibB.fa > interleavedlibsA+B.fa```

 -u SINGLE_READFILE    input fasta/fastq (fa, fa.gz, fq or fq.gz) single-read files. Only one input file of this type can be used: if you have several single-read librairies, concatenate them to generate the SINGLE_READFILE input fasta).

  -s KMER_SOLIDITY      an integer, k-mers present strictly less than this
                        number of times in the dataset will be discarded
                        (default 2)

  -S KMER_COVERAGE      an integer, minimal unitig coverage for first cleaning
                        (default 5)

  -p SR_SOLIDITY        an integer, super-reads present strictly less than
                        this number of times will be discarded (default 3)

  -P SR_COVERAGE        an integer X, unitigs with less than X reads
                        mapped is filtred (default 3)

  -k K_MIN              an integer, smallest k-mer size (default 63)

  -K K_MAX              an integer, largest k-mer size (default 201). Advice: put a -K near slightly below your read size (241 for 250 bp reads for example)

  -e MAPPING_EFFORT     number of anchors to test for mapping (default max)

  -a ANCHOR_SIZE        size of the anchors (default 31)

  -m MISSMATCH_ALLOWED  number of missmatchs allowed in mapping (default 0)

  -t NB_CORES           number of cores used (default 1)

  -o OUT_DIR            path to store the results (default = current
                        directory)

  --version             show program's version number and exit
