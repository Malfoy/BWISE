#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <unistd.h>
#include <vector>
#include <stdlib.h>
#include <unordered_map>
#include <algorithm>
#include <atomic>



using namespace std;



char revCompChar(char c) {
	switch (c) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
	}
	return 'A';
}



string revComp(const string& s){
	string rc(s.size(),0);
	for (int i((int)s.length() - 1); i >= 0; i--){
		rc[s.size()-1-i] = revCompChar(s[i]);
	}
	return rc;
}



string getCanonical(const string& str){
	return (min(str,revComp(str)));
}



void canonicalVector(vector<int64_t>& V){
	vector<int64_t> RC;
	for(uint i(0);i<V.size();++i){
		RC.push_back(-V[V.size()-i-1]);
	}
	if(V>RC){
		V=RC;
	}
}


vector<int64_t> reverseVector(const vector<int64_t>& V){
	vector<int64_t> RC;
	for(uint i(0);i<V.size();++i){
		RC.push_back(-V[V.size()-i-1]);
	}
	return RC;
}


string compactionEndNoRC(const string& seq1,const string& seq2, uint k){
	uint s1(seq1.size()),s2(seq2.size());
	if(s1==0 or s2==0){return "";}
	string end1(seq1.substr(s1-k,k)), beg2(seq2.substr(0,k));
	if(end1==beg2){return seq1+(seq2.substr(k));}
	//~ string begrc2(rc2.substr(0,k));
	//~ if(end1==begrc2){return seq1+(rc2.substr(k));}
	cout<<"fail"<<endl;
	return "";
}



bool isPrefix(const vector<int64_t>& v1, const vector<int64_t>& v2){
	if(v2.size()<v1.size()){return false;}
	for(uint i(0);i<v1.size();++i){
		if(v1[i]!=v2[i]){
			return false;
		}
	}
	return true;
}


bool isSuffix(const vector<int64_t>& v1, const vector<int64_t>& v2){
	if(v2.size()<v1.size()){return false;}
	for(uint64_t i(0);i<v1.size();++i){
		if(v1[i]!=v2[i+v2.size()-v1.size()]){
			return false;
		}
	}
	return true;
}


bool isInclued(const vector<int64_t>& v1, const vector<int64_t>& v2){
	if(v2.size()<v1.size()){return false;}
	for(uint j(0);j<v2.size();++j){
		bool success(true);
		for(uint i(0);i<v1.size();++i){
			if(i+j>v2.size()){return false;}
			if(v1[i]!=v2[i+j]){
				success=false;
				break;
			}
		}
		if(success){
			return true;
		}
	}
	return false;
}



void help(){
	cout<<"./numbersFilter  numbers.txt threshold outputFile [core to use] [superReadsThreshold] [unitig.fa] [k] [affine]"<<endl;
}


//THIS CAN WORK ON BUCKETS
int main(int argc, char *argv[]) {
	if(argc<4){
		help();
		exit(0);
	}

	vector<vector<int64_t>> lines;
	vector<int64_t> coucouch;
	vector<uint> sizeUnitig;
	string seqFile(argv[1]),unitigFile;
	uint threshold(stoi(argv[2])),superThreshold(0),kmerSize,afineThreshold(20),coreUsed(4);
	bool headerNeed(false);
	if(argc>4){
		coreUsed=(stoi(argv[4]));
	}
	if(argc>5){
		superThreshold=(stoi(argv[5]));
	}
	if(argc>7){
		unitigFile=((argv[6]));
		kmerSize=(stoi(argv[7]));
	}
	//~ if(argc>7){
		//~ headerNeed=true;
	//~ }
	if(argc>8){
		afineThreshold=(stoi(argv[8]));
	}
	int64_t uNumber;
	string line,useless,msp,number;
	ifstream numStream(seqFile);
	vector<uint64_t> count;
	vector<vector<uint64_t>> unitigsToReads;
	cout<<"Loading unitigs"<<endl;
	if(unitigFile!=""){
		//~ sizeUnitig.push_back(0);
		ifstream unitigStream(unitigFile);
		while(not unitigStream.eof()){
			getline(unitigStream,useless);
			getline(unitigStream,line);
			int size(0);
			size+=line.size();
			size-=1*(kmerSize-1);
			if(size>0){
				sizeUnitig.push_back((uint)size/afineThreshold);
			}else{
				sizeUnitig.push_back(0);
			}
		}
	}

	//LOADING and Counting
	//~ TODO WHY WE COUNT TWO TIMES ?
	while(not numStream.eof()){
		getline(numStream,useless);
		getline(numStream,line);
		coucouch={};
		if(line.size()>1){
			uint64_t i(1),lasti(0);
			while(i<line.size()){
				if(line[i]==';'){
					number=line.substr(lasti,i-lasti);
					lasti=i+1;
					uNumber=stoi(number);
					coucouch.push_back(uNumber);
					uint64_t uUNumber(uNumber>0?uNumber:-uNumber);
					if(uUNumber>count.size()){
						count.resize(uUNumber,0);
					}
					count[uUNumber-1]++;
				}
				++i;
			}
			if(coucouch.size()!=0){
				lines.push_back(coucouch);
			}
		}
	}

	//CLEANING
	cout<<"Cleaning low coverage unitigs"<<endl;
	for(uint64_t i(0);i<lines.size();++i){
		for(uint64_t j(0);j<lines[i].size();++j){
			uNumber=(lines[i][j]);
			uint64_t uUNumber(uNumber>0?uNumber:-uNumber);
			if(unitigFile!=""){
				if(count[uUNumber-1]<threshold+(sizeUnitig[uUNumber-1])){
					lines[i]={};
				}else{
				}
			}else{
				if(count[uUNumber-1]<threshold){
					lines[i]={};
				}
			}
		}
		canonicalVector(lines[i]);
	}

	//DEDUPLICATING
	cout<<"Removing Duplicate"<<endl;
	sort(lines.begin(),lines.end());
	uint64_t pred(0),counter(1);
	for(uint64_t i(1);i<lines.size();++i){
		if(lines[i]==lines[pred]){
			++counter;
			lines[i]={};
		}else{
			if(counter<superThreshold){
				lines[pred]={};
			}else{
				if(isPrefix(lines[pred],lines[i]) or isPrefix(reverseVector(lines[pred]),lines[i])){
					lines[pred]={};
				}
			}
			pred=i;
			counter=1;
		}
	}
	sort(lines.begin(),lines.end());
	lines.erase( unique( lines.begin(), lines.end() ), lines.end() );
	unitigsToReads.resize(count.size()+1,{});

	cout<<"Computing MSR"<<endl;
	//FILLING
	for(uint64_t i(0);i<lines.size();++i){
		for(uint64_t j(0);j<lines[i].size();++j){
			unitigsToReads[abs(lines[i][j])].push_back(i);
		}
	}

	//MSR COMPUTATION
	ofstream outputFile;
    outputFile.open(argv[3]);
	atomic<uint64_t> counterMSR(0),counterSR(0);

	#pragma omp parallel for num_threads(coreUsed)
	for(uint64_t i=(0);i<lines.size();++i){
		bool toPrint(true);
		unordered_map<uint64_t,int64_t> readScore;
		unordered_map<uint64_t,int64_t> stillCandidate;
		for(uint64_t j(0);j<lines[i].size() and toPrint and (j==0 or readScore.size()>0) ;++j){
			stillCandidate=readScore;
			for(uint64_t ii(0); ii<unitigsToReads[abs(lines[i][j])].size() and toPrint; ++ii){
				uint64_t friendRead=unitigsToReads[abs(lines[i][j])][ii];
				if(friendRead!=i and (j==0 or readScore.count(friendRead)!=0) ){
					++readScore[friendRead];
					stillCandidate[friendRead]=-1;
				}
				if(readScore[friendRead]>=(int)lines[i].size()){
					if( isInclued( lines[i], lines[friendRead] ) or  isInclued( reverseVector(lines[i]), lines[friendRead] )  ){
						toPrint=false;
					}
				}
			}
			for (auto it : stillCandidate){
				if(it.second>0){
					readScore.erase(it.first);
				}
			}
		}
		if(++counterSR%1000==0){
			cout<<100*counterSR/lines.size()<<"% done"<<endl;
		}
		if(toPrint){
			if(lines[i].size()>=1){
				#pragma omp critical(dataupdate)
				{

					if(headerNeed){
						outputFile<<">"+to_string(counterMSR++)<<endl;
					}

					for(uint j(0);j<lines[i].size();++j){
						outputFile<<lines[i][j]<<";";
					}
					outputFile<<"\n";
				}
			}
		}
	}

	//~ cout<<"Output"<<endl;
    //~ ofstream outputFile;
    //~ outputFile.open(argv[3]);
	//~ //OUTPUT
	//~ counter=(0);
	//~ for(uint i(0);i<lines.size();++i){
		//~ if(lines[i].size()>=1){
			//~ if(headerNeed){
				//~ outputFile<<">"+to_string(counter++)<<endl;
			//~ }
			//~ for(uint j(0);j<lines[i].size();++j){
				//~ outputFile<<lines[i][j]<<";";
			//~ }
			//~ outputFile<<endl;
		//~ }
	//~ }
    outputFile.close();
    cout<<"End"<<endl;
    return 0;
}
