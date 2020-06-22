#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <unistd.h>
#include <vector>
#include <stdlib.h>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <atomic>
#include "zstr.hpp"




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
	cout<<"./get_max  counted_path k  cores out_file  "<<endl;
}


bool almost_equal(const vector<int64_t>& V1,const vector<int64_t>& V2){
	if(V1.size()!=V2.size()){
		return false;
	}
	for(uint i(0);i<V1.size()-1;++i){
		if(V1[i]!=V2[i]){
			return false;
		}
	}
	return V2[V2.size()-1]==0;
}

void print_vector(vector<int64_t>& path){
	for (auto i = path.begin(); i != path.end(); ++i){
		cout << *i << ';';
	}
	cout<<endl;
}


//THIS CAN WORK ON BUCKETS
int main(int argc, char *argv[]) {
	if(argc<4){
		help();
		exit(0);
	}

	vector<vector<int64_t>> lines;
	vector<pair<vector<int64_t>,uint32_t>> kmers;
	vector<int64_t> coucouch;
	vector<uint> sizeUnitig;
	uint64_t MaxUnitigNumber(0);
	string seqFile(argv[1]),compactedFile("");
	uint kmer_size(stoi(argv[2])),coreUsed(stoi(argv[3]));

	string out_file=argv[4];
	bool headerNeed(false);
	if(coreUsed==0){
		cout<<"lets go ZERO"<<endl;
	}
	cout<<seqFile<<endl;
	cout<<out_file<<endl;
	uint64_t uNumber;
	string line,useless,msp,number;
	zstr::ifstream numStream(seqFile);
	cout<<"open1"<<endl;
	vector<vector<uint64_t>> unitigsToReads;
	zstr::ofstream outputFile(out_file,ofstream::app);
	cout<<"I am path_to_kmer"<<endl;
	//LOADING and Counting the countedpathWEAK
	while(not numStream.eof()){
		getline(numStream,line);
		if(line.empty()){
			break;
		}
		uint abundanceSR(stoi(line));

		getline(numStream,line);
		coucouch={};
		if(line.size()>=1){
			uint64_t i(1),lasti(0);
			while(i<line.size()){
				if(line[i]==';'){
					number=line.substr(lasti,i-lasti);
					lasti=i+1;
					uNumber=stoi(number);
					coucouch.push_back(uNumber);
				}
				++i;
			}
			if(coucouch.size()!=0){
				//WE SPLIT THE SR INTO KMERS
				for(uint j(0);j+kmer_size<=coucouch.size();++j){
					vector<int64_t> kmer(&coucouch[j],&coucouch[j+kmer_size]);
					canonicalVector(kmer);
					kmers.push_back({kmer,abundanceSR});
				}
			}
		}
	}


	sort(kmers.begin(),kmers.end());

	uint count(0);
	for(uint i(0);i<kmers.size();++i){
		count+=kmers[i].second;
		if(kmers[i].first!=kmers[i+1].first){
			if(kmers[i].first.size()>0){
				outputFile<<""+to_string(count)<<"\n";
				for(uint j(0);j<kmers[i].first.size();++j){
					outputFile<<kmers[i].first[j]<<";";
				}
				outputFile<<"\n";
				count=0;
			}
		}
	}

	outputFile<<flush;
	cout<<"End"<<endl;
	return 0;
}
