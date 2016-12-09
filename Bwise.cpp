#include <iostream>
#include <chrono>
#include <string>
#include <unistd.h>
#include <vector>



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
		;
}


int main(int argc, char *argv[]) {
	if(argc<2){

	}
	string pairedFile(""),unPairedFile(""),workingFolder(""),prefixCommand(""),folderStr(STR(folder));
	uint kMax(220),solidity(2),superReadsCleaning(3);
	if(folderStr!=""){
		prefixCommand=folderStr+"/";
	}
	//TODO working folder
	//TODO cleaning of working folder
	char c;
	while ((c = getopt (argc, argv, "u:x:o:s:k:p:")) != -1){
	switch(c){
		case 'u':
			unPairedFile=optarg;
			break;
		case 'x':
			pairedFile=optarg;
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
		}
	}
	//~ cout<<folderStr<<endl;
	//~ cin.get();
	if(pairedFile==""){
		help();
		exit(0);
	}
	//TODO parameter reads correction number, unpaired file
	cout<<"Reads Correction"<<endl;
	c=system((prefixCommand+"Bloocoo -file "+pairedFile+" -kmer-size 31 -abundance-min 5 -out reads_corrected1.fa >>log 2>>log").c_str());
	c=system((prefixCommand+"Bloocoo -file reads_corrected1.fa  -kmer-size 63 -abundance-min 5 -out reads_corrected.fa >>log 2>>log").c_str());
	//~ c=system((prefixCommand+"Bloocoo -file reads_corrected2.fa  -kmer-size 127 -abundance-min 5 -out reads_corrected.fa >>log 2>>log").c_str());
	cout<<"Reads Correction ended"<<endl;

	//TODO better kmerlist
	vector<string> kmerList{"50","100","150","200",to_string(kMax)};
	string fileBcalm(pairedFile),kmerSize;

	for(uint i(0);i<kmerList.size();++i){
		kmerSize=kmerList[i];
		cout<<"Graph construction "+to_string(i)<<endl;
		c=system((prefixCommand+"bcalm -in "+fileBcalm+" -kmer-size "+kmerSize+" -abundance-min "+to_string(solidity)+" -out out >>log 2>>log").c_str());
		c=system((prefixCommand+"kMILL out.unitigs.fa $(("+kmerSize+"-1)) $(("+kmerSize+"-2))>>log 2>>log").c_str());
		c=system(("mv out_out.unitigs.fa.fa out"+to_string(i)+".fa").c_str());
		cout<<"Graph construction "+to_string(i)+"  ended"<<endl;

		cout<<"Read mapping on the graph "+to_string(i)<<endl;
		c=system((prefixCommand+"bgreat -k "+kmerSize+" "+pairedFile+" -g out.fa -t 20  -c -m 0 -e 1 >>log 2>>log").c_str());
		cout<<"Read mapping on the graph "+to_string(i)+" ended "<<endl;
		fileBcalm="paths";
	}

	cout<<"SuperReads Cleaning"<<endl;
	c=system((prefixCommand+"pathsCleaner paths "+to_string(superReadsCleaning)+" >>log 2>>log").c_str());
	c=system(("cat out"+to_string(kmerList.size()-1)+".fa >> noduplicate.fa").c_str());
	c=system((prefixCommand+"dsk -file noduplicate.fa -kmer-size 63 -abundance-min 1 -out out_dsk >>log 2>>log").c_str());
	c=system(("echo noduplicate.fa > bankBready"));
	c=system((prefixCommand+"BREADY -graph out_dsk -bank bankBready -query bankBready -out maximalSuperReads.fa -kmer_threshold 1 -fingerprint_size 60 -core 0 -gamma 10 >>log 2>>log").c_str());
	cout<<"SuperReads Cleaning ended"<<endl;

	cout<<"SuperReads Compaction"<<endl;
	c=system((prefixCommand+"kMILL maximalSuperReads.fa >>log 2>>log").c_str());
	c=system(("mv out_maximalSuperReads.fa.fa contigs.fa >>log 2>>log").c_str());
	cout<<"SuperReads Cleaning ended"<<endl;

	cout<<"The end"<<endl;

    return 0;
}
