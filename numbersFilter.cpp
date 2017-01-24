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



void help(){
	cout<<"./numbersToSequence  numbers.txt threshold"<<endl;
}


int main(int argc, char *argv[]) {
	if(argc<3){
		help();
		exit(0);
	}
	vector<vector<int>> lines;
	vector<int> coucouch;
	string seqFile(argv[1]);
	uint threshold(stoi(argv[2]));
	int uNumber;
	string line,useless,msp,number;
	ifstream numStream(seqFile);
	vector<uint> count;

	//LOADING and Counting
	while(not numStream.eof()){
		getline(numStream,line);
		coucouch={};
		if(line.size()>3){
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
			if(coucouch.size()){
				lines.push_back(coucouch);
			}
		}
	}

	//CLEANING
	for(uint i(0);i<lines.size();++i){
		for(uint j(0);j<lines[i].size();++j){
			if(count[lines[i][j]]<threshold){
				lines[i]={};
			}
		}
		if(lines[i][0]<0){
			for(uint j(0);j<lines[i].size();++j){
				lines[i][j]=-lines[i][j];
			}
			reverse(lines[i].begin(),lines[i].end());
		}
	}

	//DEDUPLICATING
	sort(lines.begin(),lines.end());
	//IS THIS THE GOOD SOLUTION ?
	lines.erase( unique( lines.begin(), lines.end() ), lines.end() );

	//OUTPUT
	for(uint i(0);i<lines.size();++i){
		for(uint j(0);j<lines.size();++j){
			cout<<lines[i][j]<<";";
		}
		cout<<endl;
	}

    return 0;
}
