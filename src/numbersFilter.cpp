#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <unistd.h>
#include <vector>
#include <stdlib.h>
#include <unordered_map>
#include <algorithm>



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



void canonicalVector(vector<int>& V){
	vector<int> RC;
	for(uint i(0);i<V.size();++i){
		RC.push_back(-V[V.size()-i-1]);
	}
	if(V>RC){
		V=RC;
	}
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



bool isInclued(const vector<int>& v1, const vector<int>& v2){
	if(v2.size()<v1.size()){return false;}
	for(uint i(0);i<v1.size();++i){
		if(v1[i]!=v2[i]){
			return false;
		}
	}
	return true;
}



void help(){
	cout<<"./numbersFilter  numbers.txt threshold [superReadsThreshold] [unitig.fa] [k] [header]"<<endl;
}


int main(int argc, char *argv[]) {
	if(argc<3){
		help();
		exit(0);
	}

	vector<vector<int>> lines;
	vector<int> coucouch;
	vector<uint> sizeUnitig;
	string seqFile(argv[1]),unitigFile;
	uint threshold(stoi(argv[2])),superThreshold(0),kmerSize;
	bool headerNeed(false);
	if(argc>3){
		superThreshold=(stoi(argv[3]));
	}
	if(argc>5){
		unitigFile=((argv[4]));
		kmerSize=(stoi(argv[5]));
	}
	if(argc>6){
		headerNeed=true;
	}
	int uNumber;
	string line,useless,msp,number;
	ifstream numStream(seqFile);
	vector<uint> count;
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
				sizeUnitig.push_back((uint)size/10);
			}else{
				sizeUnitig.push_back(0);
			}
		}
	}

	//LOADING and Counting
	while(not numStream.eof()){
		getline(numStream,useless);
		getline(numStream,line);
		coucouch={};
		if(line.size()>1){
			uint i(1),lasti(0);
			while(i<line.size()){
				if(line[i]==';'){
					number=line.substr(lasti,i-lasti);
					lasti=i+1;
					uNumber=stoi(number);
					coucouch.push_back(uNumber);
					uint uUNumber(uNumber>0?uNumber:-uNumber);
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
	for(uint i(0);i<lines.size();++i){
		for(uint j(0);j<lines[i].size();++j){
			uNumber=(lines[i][j]);
			uint uUNumber(uNumber>0?uNumber:-uNumber);
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
	sort(lines.begin(),lines.end());
	if(superThreshold==0){
		lines.erase( unique( lines.begin(), lines.end() ), lines.end() );
	}else{
		uint pred(0),count(1);
		for(uint i(1);i<lines.size();++i){
			if(lines[i]==lines[pred]){
				++count;
				lines[i]={};
			}else{
				if(count<superThreshold){
					lines[pred]={};
				}
				pred=i;
				count=1;
			}
		}
		if(not lines.empty()){
			if(count<superThreshold){
				lines[pred]={};
			}
		}
	}
	sort(lines.begin(),lines.end());
	lines.erase( unique( lines.begin(), lines.end() ), lines.end() );
	//~ uint pred(0);
	for(uint i(1);i<lines.size();++i){
		if(isInclued(lines[i-1],lines[i])){
			lines[i-1]={};
		}
	}

	//OUTPUT
	uint counter(0);
	for(uint i(0);i<lines.size();++i){
		if(lines[i].size()>=1){
			if(headerNeed){
				cout<<">"+to_string(counter++)<<endl;
			}
			for(uint j(0);j<lines[i].size();++j){
				cout<<lines[i][j]<<";";
			}
			cout<<endl;
		}
	}
    return 0;
}
