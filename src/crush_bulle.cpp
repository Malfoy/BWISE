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
	cout<<"./crush_bulle  numbers.txt  outputFile core to use "<<endl;
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
	string seqFile(argv[1]);
	uint coreUsed(stoi(argv[3]));
	bool headerNeed(false);

	int64_t uNumber;
	string line,useless,msp,number;
	ifstream numStream(seqFile);
	vector<uint64_t> count;
	vector<vector<uint64_t>> unitigsToReads;

	//LOADING and Counting
	//~ TODO WHY WE COUNT TWO TIMES ?
	while(not numStream.eof()){
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
				canonicalVector(coucouch);
				lines.push_back(coucouch);
			}
		}
	}
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
    outputFile.open(argv[2]);
	atomic<uint64_t> counterMSR(0),counterSR(0);
	cout<<"go"<<endl;
	#pragma omp parallel num_threads(coreUsed)
	{
		#pragma omp for
		for(uint64_t i=(0);i<lines.size();++i){
			bool toPrint(true);
			for(uint64_t ii(0); ii<unitigsToReads[abs(lines[i][0])].size() and toPrint; ++ii){
				uint64_t friendRead=unitigsToReads[abs(lines[i][0])][ii];
				if(friendRead!=i and lines[friendRead][0]==lines[i][0] and lines[friendRead][lines[friendRead].size()-1]==lines[i][lines[i].size()-1]
				 //~ and lines[friendRead].size()==lines[i].size()
				 ){
					if(lines[i]>lines[friendRead]){
						toPrint=false;
						//~ for(uint j(0);j<lines[i].size();++j){
							//~ cout<<lines[i][j]<<";";
						//~ }
						//~ cout<<endl;
						//~ cout<<"because"<<endl;
						//~ for(uint j(0);j<lines[friendRead].size();++j){
							//~ cout<<lines[friendRead][j]<<";";
						//~ }
						//~ cout<<endl<<endl;;
					}
				}
			}


			if(++counterSR%10000==0){
				cout<<100*counterSR/lines.size()<<"% done"<<endl;
			}
			if(toPrint){
				if(lines[i].size()>=1){
					#pragma omp critical(dataupdate)
					{
						if(headerNeed){
							outputFile<<">"+to_string(counterMSR++)<<"\n";
						}

						for(uint j(0);j<lines[i].size();++j){
							outputFile<<lines[i][j]<<";";
						}
						outputFile<<"\n";
					}
				}
			}
		}
	}

    cout<<"End"<<endl;
    return 0;
}
