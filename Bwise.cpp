#include <iostream>
#include <fstream>
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
	<<"Options and default values: "<<endl
	<<"-x for paired read file"<<endl
	<<"-u for unpaired read file"<<endl
	<<"-o for working folder (.)"<<endl
	<<"-s for kmer solidity threshold (5)"<<endl
	<<"-k for largest kmer size (220)"<<endl
	<<"-p for superReads cleaning threshold (2)"<<endl
	<<"-c for correction step (max)"<<endl
	<<"-t for core used (max)"<<endl
	;
}



int main(int argc, char *argv[]) {
	if(argc<2){
		help();
		exit(0);
	}
	string pairedFile(""),unPairedFile(""),workingFolder("."),prefixCommand(""),folderStr(STR(folder)),bgreatArg,bloocooArg;
	uint kMax(220),solidity(5),superReadsCleaning(2),correctionStep(3),coreUsed(0);
	if(folderStr!=""){
		prefixCommand=folderStr+"/";
	}
	char c;
	while ((c = getopt (argc, argv, "u:x:o:s:k:p:c:t:")) != -1){
	switch(c){
		case 'u':
			unPairedFile=realpath(optarg,NULL);
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
		case 't':
			coreUsed=stoi(optarg);
			break;
		}
	}
	if(pairedFile=="" and unPairedFile==""){
		help();
		exit(0);
	}

	c=system(("mkdir "+workingFolder).c_str());
	c=chdir(workingFolder.c_str());
	c=system("mkdir logs");
	ofstream param("ParametersUsed.txt");
	ofstream bankBcalm("bankBcalm.txt");
	param<<"kmax: "<<kMax<<" solidity: "<<solidity<<" SRcleaning: "<<superReadsCleaning<<" correction steps: "<<correctionStep<<endl;
	uint filesCase(0);
	if(pairedFile==""){
		filesCase=1;
		bloocooArg=unPairedFile;
		bgreatArg=" -u reads_corrected.fa ";
		bankBcalm<<"reads_corrected.fa"<<endl;
	}else{
		if(unPairedFile==""){
			filesCase=2;
			bloocooArg=pairedFile;
			bgreatArg=" -x reads_corrected.fa ";
			bankBcalm<<"reads_corrected.fa"<<endl;
		}else{
			filesCase=3;
			bloocooArg=pairedFile+","+unPairedFile;
			bankBcalm<<"reads_corrected1.fa"<<endl<<"reads_corrected2.fa"<<endl;
			bgreatArg=" -x reads_corrected1.fa  -u reads_corrected2.fa ";
		}
	}
	cout<<filesCase<<endl;


	cout<<"Reads Correction"<<endl;
	string fileToCorrect(pairedFile);
	vector<string> kmerSizeCorrection={"31","63","127"};
	vector<string> bloocooversion={"32","64","128"};
	uint indiceCorrection(0);
	for(;indiceCorrection<min(correctionStep,(uint)kmerSizeCorrection.size());++indiceCorrection){
		c=system((prefixCommand+"Bloocoo"+bloocooversion[indiceCorrection]+" -file "+bloocooArg+" -kmer-size "+kmerSizeCorrection[indiceCorrection]+" -abundance-min "+to_string(solidity)+" -out reads_corrected"+to_string(indiceCorrection)+".fa -nb-cores "+to_string(coreUsed)+"  >>logs/logBloocoo 2>>logs/logBloocoo").c_str());
		if(filesCase==3){
			c=system(("mv reads_corrected"+to_string(indiceCorrection)+"_0_.fasta reads_corrected"+to_string(indiceCorrection)+"1.fa").c_str());
			c=system(("mv reads_corrected"+to_string(indiceCorrection)+"_1_.fasta reads_corrected"+to_string(indiceCorrection)+"2.fa").c_str());
			bloocooArg="reads_corrected"+to_string(indiceCorrection)+"1.fa,reads_corrected"+to_string(indiceCorrection)+"2.fa";
		}else{
			bloocooArg="reads_corrected"+to_string(indiceCorrection)+".fa ";
		}
	}
	if(correctionStep==0){
		if(filesCase==3){
			c=system(("ln -s "+pairedFile+" reads_corrected1.fa").c_str());
			c=system(("ln -s "+unPairedFile+" reads_corrected2.fa").c_str());
		}else{
			c=system(("ln -s "+bloocooArg+" reads_corrected.fa").c_str());
		}
	}else{
		if(filesCase==3){
			c=system(("ln -s reads_corrected"+to_string(indiceCorrection-1)+"1.fa  reads_corrected1.fa").c_str());
			c=system(("ln -s reads_corrected"+to_string(indiceCorrection-1)+"2.fa  reads_corrected2.fa").c_str());
		}else{
			c=system(("mv "+bloocooArg+" reads_corrected.fa").c_str());
		}
	}


	//TODO better kmerlist
	vector<string> kmerList{"50","100","150","200","210","220","230","240","250","260","270","280","290","300","310"};
	string fileBcalm("bankBcalm.txt"),kmerSize;
	uint indiceGraph(0);
	for(;indiceGraph<kmerList.size() and stoi(kmerList[indiceGraph])<=kMax;++indiceGraph){
		kmerSize=kmerList[indiceGraph];
		cout<<"Graph construction "+to_string(indiceGraph)<<endl;
		string kmerSizeTip(to_string(min((stoi(kmerSize)+50),250)));
		c=system((prefixCommand+"bcalm -in "+fileBcalm+" -kmer-size "+kmerSize+" -abundance-min "+to_string(solidity)+" -out out  -nb-cores "+to_string(coreUsed)+"  >>logs/logBcalm 2>>logs/logBcalm").c_str());
		c=system((prefixCommand+"kMILL out.unitigs.fa $(("+kmerSize+"-1)) $(("+kmerSize+"-2)) >>logs/logBcalm 2>>logs/logBcalm").c_str());
		c=system((prefixCommand+"tipCleaner out_out.unitigs.fa.fa $(("+kmerSize+"-1)) "+kmerSizeTip+" >>logs/logTip 2>>logs/logTip").c_str());
		//~ c=system(("mv tiped.fa swag"));
		//~ c=system((prefixCommand+"tipCleaner swag $(("+kmerSize+"-1)) $(("+kmerSize+"+50)) >>logs/logTip 2>>logs/logTip").c_str());
		//~ c=system((prefixCommand+"tipCleaner out_out.unitigs.fa.fa $(("+kmerSize+"-1)) 1 >>logs/logTip 2>>logs/logTip").c_str());
		c=system((prefixCommand+"kMILL tiped.fa $(("+kmerSize+"-1)) $(("+kmerSize+"-2)) >>logs/logBcalm 2>>logs/logBcalm").c_str());
		c=system(("mv out_tiped.fa.fa dbg"+to_string(indiceGraph)+".fa").c_str());

		cout<<"Read mapping on the graph "+to_string(indiceGraph)<<endl;
		c=system((prefixCommand+"bgreat -k "+kmerSize+" "+bgreatArg+" -g dbg"+to_string(indiceGraph)+".fa -t "+to_string((coreUsed==0)?10:coreUsed) +"  -c -m 0 -e 1 >>logs/logBgreat 2>>logs/logBgreat").c_str());
		fileBcalm="paths";
		solidity=2;
		if(stoi(kmerSize)>=200){
			solidity=1;
		}
	}


	cout<<"SuperReads Cleaning"<<endl;
	//~ cout<<"sequencesCleaner paths "+to_string(superReadsCleaning)<<endl;
	c=system((prefixCommand+"sequencesCleaner paths "+to_string(superReadsCleaning)+" >>logs/logBready 2>>logs/logBready").c_str());
	c=system(("cat dbg"+to_string(indiceGraph-1)+".fa >> noduplicate.fa").c_str());
	c=system((prefixCommand+"dsk -file noduplicate.fa -kmer-size 31 -abundance-min 1 -out out_dsk -nb-cores "+to_string(coreUsed)+"  >>logs/logBready 2>>logs/logBready").c_str());
	c=system(("echo noduplicate.fa > bankBready"));
	c=system((prefixCommand+"BREADY -graph out_dsk -bank bankBready -query bankBready -out maximalSuperReads.fa -kmer_threshold 1 -fingerprint_size 8 -core "+to_string(coreUsed)+"  -gamma 10 >>logs/logBready 2>>logs/logBready").c_str());
	//cout<<"SuperReads Cleaning ended"<<endl;

	cout<<"SuperReads Compaction"<<endl;
	c=system((prefixCommand+"kMILL maximalSuperReads.fa >>logs/logkmill 2>>logs/logkmill").c_str());
	c=system(("mv out_maximalSuperReads.fa.fa contigs.fa >>logs/logkmill 2>>logs/logkmill"));
	//cout<<"SuperReads Cleaning ended"<<endl;
	c=system(("rm -rf trashme* *.h5 out.unitigs.fa notAligned.fa bankBready bankBcalm >>logs/logkmill 2>>logs/logkmill"));
	cout<<"The end"<<endl;

    return 0;
}
