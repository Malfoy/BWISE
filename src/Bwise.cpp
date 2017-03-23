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
using namespace chrono;



void help(){
	cout<<"BWISE"<<endl
	<<"Options and default values: "<<endl
	<<"-x for paired read file"<<endl
	<<"-u for unpaired read file"<<endl
	<<"-s for kmer solidity threshold (2)"<<endl
	<<"-S for unitig solidity threshold (2)"<<endl
	<<"-o for working folder (.)"<<endl
	<<"-k for largest kmer size (301)"<<endl
	<<"-p for superReads cleaning threshold (2)"<<endl
	<<"-c for correction step (max)"<<endl
	<<"-t for core used (max)"<<endl
	;
}



bool exists_test (const string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}



int main(int argc, char *argv[]) {
	if(argc<2){
		help();
		exit(0);
	}
	string pairedFile(""),unPairedFile(""),workingFolder("."),prefixCommand(""),folderStr(STR(folder)),bgreatArg,bloocooArg,slowParameter(" -slow ");
	uint kMax(301),solidity(2),superReadsCleaning(2),correctionStep(4),coreUsed(0),unitigFilter(2);
	if(folderStr!=""){
		prefixCommand=folderStr+"/";
	}
	char c;
	while ((c = getopt (argc, argv, "u:x:o:s:k:p:c:t:S:")) != -1){
	switch(c){
		case 'u':
			if(not exists_test(optarg)){
				cout<<optarg<<" not found..."<<endl;
				exit(0);
			}
			unPairedFile=realpath(optarg,NULL);
			break;
		case 'x':
			if(not exists_test(optarg)){
				cout<<optarg<<" not found..."<<endl;
				exit(0);
			}
			pairedFile=realpath(optarg,NULL);
			break;
		case 'o':
			workingFolder=(optarg);
			break;
		case 's':
			solidity=stoi(optarg);
			break;
		case 'S':
			unitigFilter=stoi(optarg);
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
	param<<"kmax: "<<kMax<<" solidity: "<<solidity<<" unitig solidity: "<<unitigFilter<<" SRcleaning: "<<superReadsCleaning<<" correction steps: "<<correctionStep<<endl;
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



	//CORRECTION
	cout<<"Reads Correction... "<<flush;
	auto start=system_clock::now();
	auto realStart=start;
	string fileToCorrect(pairedFile);
	vector<string> kmerSizeCorrection={"31","63","95","127"};
	vector<string> bloocooversion={"32","64","128","128"};
	uint indiceCorrection(0);
	for(;indiceCorrection<min(correctionStep,(uint)kmerSizeCorrection.size());++indiceCorrection){
		c=system((prefixCommand+"Bloocoo"+bloocooversion[indiceCorrection]+" -file "+bloocooArg+" "+slowParameter+"   -kmer-size "+kmerSizeCorrection[indiceCorrection]+" -nbits-bloom 24  -out reads_corrected"+to_string(indiceCorrection)+".fa -nb-cores "+to_string(coreUsed)+"  >>logs/logBloocoo 2>>logs/logBloocoo").c_str());
		//~ c=system((prefixCommand+ "h5dump -y -d histogram/histogram  reads_corrected"+to_string(indiceCorrection)+".fa.h5  > logs/histocorr"+to_string(indiceCorrection)).c_str());
		c=system(("rm  reads_corrected"+to_string(indiceCorrection-1)+"* 2>> logs/histocorr"+to_string(indiceCorrection)).c_str());

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
	auto end=system_clock::now();
	cout<<"Correction took "<<duration_cast<minutes>(end-start).count()<<" minutes"<<endl;



	//TODO better kmerlist
	vector<string> kmerList{"0","51","101","151","201","251","301","351","401","451","501"};
	string fileBcalm("bankBcalm.txt"),kmerSize;
	uint indiceGraph(1);
	for(;indiceGraph<kmerList.size() and (uint)stoi(kmerList[indiceGraph])<=kMax;++indiceGraph){
		kmerSize=kmerList[indiceGraph];
		cout<<"Graph construction "+to_string(indiceGraph)<<"... "<<flush;
		start=system_clock::now();
		string kmerSizeTip(to_string((2*stoi(kmerSize))));
		//GRAPH CONSTRUCTION
		c=system((prefixCommand+"bcalm -in "+fileBcalm+" -kmer-size "+kmerSize+" -abundance-min "+to_string(solidity)+" -out out  -nb-cores "+to_string(coreUsed)+"  >>logs/logBcalm 2>>logs/logBcalm").c_str());
		c=system((prefixCommand+ "h5dump -y -d histogram/histogram  out.h5  > logs/histodbg"+(kmerSize)).c_str());
		//GRAPH CLEANING
		c=system((prefixCommand+"kMILL out.unitigs.fa $(("+kmerSize+"-1)) $(("+kmerSize+"-2)) >>logs/logBcalm 2>>logs/logBcalm").c_str());
		c=system((prefixCommand+"tipCleaner out_out.unitigs.fa.fa $(("+kmerSize+"-1)) "+kmerSizeTip+" >>logs/logTip 2>>logs/logTip").c_str());
		c=system((prefixCommand+"kMILL tiped.fa $(("+kmerSize+"-1)) $(("+kmerSize+"-2)) >>logs/logBcalm 2>>logs/logBcalm").c_str());
		c=system(("mv out_tiped.fa.fa dbg"+to_string(indiceGraph)+".fa").c_str());
		cout<<"Read mapping on the graph "+to_string(indiceGraph)<<"... "<<flush;
		//READ MAPPING
		c=system((prefixCommand+"bgreat -k "+kmerSize+" -M -i 10 "+bgreatArg+" -g dbg"+to_string(indiceGraph)+".fa -t "+to_string((coreUsed==0)?20:coreUsed) +" -a 63  -m 0 -e 100 >>logs/logBgreat 2>>logs/logBgreat").c_str());
		if(indiceGraph==1){
			c=system((prefixCommand+"numbersFilter paths "+to_string(unitigFilter)+" "+to_string(superReadsCleaning)+" dbg"+to_string(indiceGraph)+".fa   $(("+kmerSize+"))  > cleanedPaths 2>>logs/logBgreat").c_str());
			c=system((prefixCommand+"numbersToSequences dbg"+to_string(indiceGraph)+".fa  cleanedPaths  $(("+kmerSize+"-1)) >newPaths 2>>logs/logBgreat").c_str());
		}else{
			if((uint)stoi(kmerList[indiceGraph])<kMax){
				c=system((prefixCommand+"numbersFilter paths "+to_string(unitigFilter)+" "+to_string(superReadsCleaning)+" dbg"+to_string(indiceGraph)+".fa   $(("+kmerSize+"))  > cleanedPaths 2>>logs/logBgreat").c_str());
				c=system(("python3 "+prefixCommand+"K2000.py cleanedPaths dbg"+to_string(indiceGraph)+".fa  $(("+kmerSize+")) -e > compacted_unitigs"+to_string(indiceGraph)+".txt  2>>logs/logK2000").c_str());
				c=system((prefixCommand+"numbersToSequences dbg"+to_string(indiceGraph)+".fa  compacted_unitigs"+to_string(indiceGraph)+".txt  $(("+kmerSize+"-1)) >newPaths 2>>logs/logBgreat").c_str());
			}else{
				c=system((prefixCommand+"numbersFilter paths "+to_string(unitigFilter)+" "+to_string(superReadsCleaning)+" dbg"+to_string(indiceGraph)+".fa   $(("+kmerSize+"))  > cleanedPaths 2>>logs/logBgreat").c_str());
				//~ c=system((prefixCommand+"numbersToSequences dbg"+to_string(indiceGraph)+".fa  cleanedPaths  $(("+kmerSize+"-1)) >newPaths 2>>logs/logBgreat").c_str());
			}
		}
		fileBcalm="newPaths";
		solidity=1;
		end=system_clock::now();
		cout<<"Step "+to_string(indiceGraph)+" took "<<duration_cast<minutes>(end-start).count()<<" minutes"<<endl;
	}


	//~ //SUPERREADS COMPACTION
	cout<<"SuperReads Compactions ..."<<flush;
	start=system_clock::now();
	//K2000
	c=system(("python3 "+prefixCommand+"K2000.py cleanedPaths dbg"+to_string(indiceGraph-1)+".fa  $(("+kmerSize+")) -e > compacted_unitigs.txt  2>>logs/logK2000").c_str());
	c=system((prefixCommand+"numbersToSequences dbg"+to_string(indiceGraph-1)+".fa  compacted_unitigs.txt $(("+kmerSize+"-1)) > contigs.fa 2>>logs/logSRC").c_str());
	//~ c=system(("python3 "+prefixCommand+"K2000_msr_to_gfa.py compacted_unitigs.txt  dbg"+to_string(indiceGraph-1)+".fa  "+(kmerSize)+" > outk2000.gfa").c_str());
	//~ c=system(("python3 "+prefixCommand+"K2000_gfa_to_fasta.py outk2000.gfa  > contigsk2000.fa").c_str());
	//BREADY and KMILL
	//~ c=system((prefixCommand+"dsk -file newPaths -kmer-size 31 -abundance-min 1 -out out_dsk -nb-cores "+to_string(coreUsed)+"  >>logs/logBready 2>>logs/logBready").c_str());
	//~ c=system(("echo newPaths > bankBready"));
	//~ c=system((prefixCommand+"BREADY -graph out_dsk -bank bankBready -query bankBready -out maximalSuperReads.fa -kmer_threshold 1 -fingerprint_size 8 -core "+to_string(coreUsed)+"  -gamma 10 >>logs/logBready 2>>logs/logBready").c_str());
	//~ c=system((prefixCommand+"kMILL maximalSuperReads.fa >>logs/logkmill 2>>logs/logkmill").c_str());
	//~ c=system(("mv out_maximalSuperReads.fa.fa contigs.fa >>logs/logkmill 2>>logs/logkmill"));
	c=system(("rm -rf trashme* *.h5 out.unitigs.fa notAligned.fa bankBready bankBcalm.txt maximalSuperReads.fa newPaths out_out.unitigs.fa.fa tiped.fa paths >>logs/logkmill 2>>logs/logkmill"));
	end=system_clock::now();
	cout<<"SuperReads compaction took "<<duration_cast<minutes>(end-start).count()<<" minutes"<<endl;
	cout<<"The end !"<<endl;
	cout<<"BWISE assembly took "<<duration_cast<minutes>(end-realStart).count()<<" minutes"<<endl;

    return 0;
}
