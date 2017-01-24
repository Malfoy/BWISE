#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <unistd.h>
#include <vector>
#include <stdlib.h>
#include <unordered_map>



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
	cout<<"./numbersToSequence unitigs.fa numbers.txt kmersize"<<endl;
}


int main(int argc, char *argv[]) {
	if(argc<4){
		help();
		exit(0);
	}
	vector<string> unitigs;
	string unitigFile(argv[1]);
	string seqFile(argv[2]);
	uint  k(stoi(argv[3]));
	string unitig,useless,msp;
	ifstream unitigStream(unitigFile);
	ifstream numStream(seqFile);
	while(not unitigStream.eof()){
		getline(unitigStream,useless);
		getline(unitigStream,unitig);
		unitigs.push_back(unitig);
	}
	string number,contig;
	int count(0),uNumber;
	while(not numStream.eof()){
		contig="";
		getline(numStream,msp);
		uint i(1),lasti(0);
		while(i<msp.size()){
			if(msp[i]==';'){
				number=msp.substr(lasti,i-lasti);
				lasti=i+1;
				uNumber=stoi(number);
				if(uNumber>0){
					unitig=unitigs[uNumber];
				}else{
					unitig=revComp(unitigs[-uNumber]);
				}
				if(contig!=""){
					contig=compactionEndNoRC(contig,unitig,k);
				}else{
					contig=unitig;
				}
			}
			++i;
		}
		cout<<">"+count<<endl<<contig<<endl;
	}

    return 0;
}
