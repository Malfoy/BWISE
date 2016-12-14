#include <iostream>
#include <chrono>
#include <string>
#include <unistd.h>
#include <vector>
#include <stdlib.h>



#define XSTR(x) #x
#define STR(x) XSTR(x)



using namespace std;


void help(){
	cout<<"BWISE"<<endl
		<<"-x for paired read file"<<endl
		<<"-u for unpaired read file"<<endl
		<<"-o for working folder"<<endl
		<<"-s for kmer solidity threshold"<<endl
		<<"-k for largest kmer size"<<endl
		<<"-p for superReads cleaning threshold"<<endl
		<<"-c for correction step"<<endl
		;
}


int main(int argc, char *argv[]) {
	if(argc<2){
		help();
	}
	string pairedFile(""),unPairedFile(""),workingFolder("."),prefixCommand(""),folderStr(STR(folder));
	uint kMax(220),solidity(2),superReadsCleaning(3),correctionStep(3);
	if(folderStr!=""){
		prefixCommand=folderStr+"/";
	}
	char c;
	while ((c = getopt (argc, argv, "u:x:o:s:k:p:c:")) != -1){
	switch(c){
		case 'u':
			unPairedFile=optarg;
			break;
		case 'x':
			pairedFile=realpath(optarg,NULL);
			break;
		case 'o':
			workingFolder=(optarg);
			break;
		case 's':
			solidity=stoi(optarg);
			break;
		case 'k':
			kMax=stoi(optarg);
			break;
		case 'p':
			superReadsCleaning=stoi(optarg);
			break;
		case 'c':
			correctionStep=stoi(optarg);
			break;
		}
	}
	if(pairedFile==""){
		help();
		exit(0);
	}
	c=chdir(workingFolder.c_str());
	c=system("mkdir logs");
	//TODO  unpaired file
	cout<<"Reads Correction"<<endl;
	string fileToCorrect(pairedFile);
	vector<string> kmerSizeCorrection={"31","63","127"};
	vector<string> bloocooversion={"32","64","128"};
	for(uint i(0);i<min(correctionStep,(uint)kmerSizeCorrection.size());++i){
		c=system((prefixCommand+"Bloocoo"+bloocooversion[i]+" -file "+fileToCorrect+" -kmer-size "+kmerSizeCorrection[i]+" -abundance-min 5 -out reads_corrected"+to_string(i)+".fa >>logs/logBloocoo 2>>logs/logBloocoo").c_str());
		fileToCorrect="reads_corrected"+to_string(i)+".fa";
	}
	if(fileToCorrect==pairedFile){
		 c=system(("cp "+fileToCorrect+" reads_corrected.fa").c_str());
	}else{
		c=system(("mv "+fileToCorrect+" reads_corrected.fa").c_str());
	}

	//cout<<"Reads Correction ended"<<endl;

	//TODO better kmerlist
	vector<string> kmerList{"50","100","150","200",to_string(kMax)};
	string fileBcalm(pairedFile),kmerSize;

	for(uint i(0);i<kmerList.size();++i){
		kmerSize=kmerList[i];
		cout<<"Graph construction "+to_string(i)<<endl;
		c=system((prefixCommand+"bcalm -in "+fileBcalm+" -kmer-size "+kmerSize+" -abundance-min "+to_string(solidity)+" -out out >>logs/logBcalm 2>>logs/logBcalm").c_str());
		c=system((prefixCommand+"kMILL out.unitigs.fa $(("+kmerSize+"-1)) $(("+kmerSize+"-2))>>logs/logBcalm 2>>logs/logBcalm").c_str());
		c=system(("mv out_out.unitigs.fa.fa dbg"+to_string(i)+".fa").c_str());
		//cout<<"Graph construction "+to_string(i)+"  ended"<<endl;

		cout<<"Read mapping on the graph "+to_string(i)<<endl;
		c=system((prefixCommand+"bgreat -k "+kmerSize+" -x reads_corrected.fa  -g dbg"+to_string(i)+".fa -t 20  -c -m 0 -e 1 >>logs/logBgreat 2>>logs/logBgreat").c_str());
		fileBcalm="paths";
	}

	cout<<"SuperReads Cleaning"<<endl;
	c=system((prefixCommand+"sequencesCleaner paths "+to_string(superReadsCleaning)+" >>logs/logBready 2>>logs/logBready").c_str());
	c=system(("cat dbg"+to_string(kmerList.size()-1)+".fa >> noduplicate.fa").c_str());
	c=system((prefixCommand+"dsk -file noduplicate.fa -kmer-size 63 -abundance-min 1 -out out_dsk >>logs/logBready 2>>logs/logBready").c_str());
	c=system(("echo noduplicate.fa > bankBready"));
	c=system((prefixCommand+"BREADY -graph out_dsk -bank bankBready -query bankBready -out maximalSuperReads.fa -kmer_threshold 1 -fingerprint_size 60 -core 0 -gamma 10 >>logs/logBready 2>>logs/logBready").c_str());
	//cout<<"SuperReads Cleaning ended"<<endl;

	cout<<"SuperReads Compaction"<<endl;
	c=system((prefixCommand+"kMILL maximalSuperReads.fa >>logs/logkmill 2>>logs/logkmill").c_str());
	c=system(("mv out_maximalSuperReads.fa.fa contigs.fa >>logs/logkmill 2>>logs/logkmill"));
	//cout<<"SuperReads Cleaning ended"<<endl;
	c=system(("rm -rf trashme* *.h5 out.unitigs.fa notAligned.fa bankBready >>logs/logkmill 2>>logs/logkmill"));
	cout<<"The end"<<endl;

    return 0;
}
