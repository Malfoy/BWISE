INSTALATION

./install.sh

Will download and compile all the needeed software


UPDATE

./install.sh

Will pull the different git and recompile them
(if this do not work, you can try to erase the folder itself and do the part proposed in install.sh for a clean build)

RUN

./assemble.sh -x examplePairedReads.fa -u exampleUnpairedReads.fa

Will run the complete pipeline in the directory folderAssembly
-x option to give paired reads
-u option to give unpaired reads

