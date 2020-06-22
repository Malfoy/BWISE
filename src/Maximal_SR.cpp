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
	cout<<"./get_max  counted_path threshold  outputfile core to use  compacted_path "<<endl;
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
		cout << *i << ' ';
	}
	cout<<endl;
}


//THIS CAN WORK ON BUCKETS
int main(int argc, char *argv[]) {
	if(argc<4){
		help();
		exit(0);
	}
	cout<<"I am get_max"<<endl;

	vector<vector<int64_t>> lines;
	vector<int64_t> coucouch;
	vector<uint> sizeUnitig;
	uint64_t MaxUnitigNumber(0);
	string seqFile(argv[1]),compactedFile("");
	uint superThreshold(stoi(argv[2])),coreUsed(stoi(argv[4]));
	bool headerNeed(false);
	if(argc>5){
		compactedFile=argv[5];
	}


	int64_t uNumber;
	string line,useless,msp,number;
	zstr::ifstream numStream(seqFile);
	vector<vector<uint64_t>> unitigsToReads;



	//LOADING and Counting the countedpath
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
					uint64_t uUNumber(uNumber>0?uNumber:-uNumber);
					if(uUNumber>MaxUnitigNumber){
						MaxUnitigNumber=uUNumber;
					}
				}
				++i;
			}
			if(coucouch.size()!=0){
				canonicalVector(coucouch);
				coucouch.push_back(abundanceSR);
				lines.push_back(coucouch);
			}
		}
	}

	if(not compactedFile.empty()){
		zstr::ifstream compactedfile(compactedFile);
		if(compactedfile.good()){
			cout<<"READING COMPACT"<<endl;
			while(not compactedfile.eof()){
				getline(compactedfile,line);
				coucouch={};
				if(line.size()>=1){
					uint64_t i(1),lasti(0);
					while(i<line.size()){
						if(line[i]==';'){
							number=line.substr(lasti,i-lasti);
							lasti=i+1;
							uNumber=stoi(number);
							coucouch.push_back(uNumber);
							uint64_t uUNumber(uNumber>0?uNumber:-uNumber);
							if(uUNumber>MaxUnitigNumber){
								MaxUnitigNumber=uUNumber;
							}
						}
						++i;
					}
					if(coucouch.size()!=0){
						canonicalVector(coucouch);
						coucouch.push_back(0);
						lines.push_back(coucouch);
					}
				}
			}
		}
	}

	unitigsToReads.resize(MaxUnitigNumber+1,{});
	sort(lines.begin(),lines.end());
	uint64_t pred(0),counter(1);
	for(uint64_t i(1);i<lines.size();++i){
		if(almost_equal(lines[i],lines[i-1])){
			lines[i-1]={};
		}else{
			if(not lines[i-1].empty()){
				lines[i-1].pop_back();
			}
		}
		if(lines[i][lines[i].size()-1]!=0){
			if(lines[i][lines[i].size()-1]<superThreshold){
				lines[i]={};
			}
		}
	}
	if(not lines[lines.size()-1].empty()){
		lines[lines.size()-1].pop_back();
	}
	sort(lines.begin(),lines.end());
	lines.erase( unique( lines.begin(), lines.end() ), lines.end() );


	//REMOVING NONMAXIMAL
	cout<<"Remove NONMAX "<<endl;

	pred=(0);
	for(uint64_t i(1);i<lines.size();++i){
		if(isPrefix(lines[pred],lines[i]) or isPrefix(reverseVector(lines[pred]),lines[i])){
			lines[pred]={};
		}
	}


	cout<<"Removing trivial msr"<<endl;
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
	cout<<"Computing MSR"<<endl;
	#pragma omp parallel num_threads(coreUsed)
	{
		unordered_set<uint64_t> unitig_set;
		vector<uint64_t> candidates;
		vector<int64_t> Line;
		vector<int64_t> ami;
		#pragma omp for
		for(uint64_t i=(0);i<lines.size();++i){//FOR EACH SR
			if(lines[i].size()>=1){
				Line=lines[i];
				unitig_set.clear();
				for(uint64_t ii(0); ii<Line.size(); ++ii){
					unitig_set.insert(abs(Line[ii]));
				}
				bool toPrint(true);
				candidates=unitigsToReads[abs(Line[0])];
				for(uint64_t ii(0); ii<candidates.size(); ++ii){
					if(candidates[ii]==i){continue;}
					ami=(lines[candidates[ii]]);
					uint score(0);
					for(uint ami_i(0);ami_i<ami.size();++ami_i){
						if(unitig_set.count(abs(ami[ami_i]))==1){
							score++;
						}
					}
					if(score>=Line.size()){
						if( isInclued( Line, ami ) or  isInclued( reverseVector(Line), ami )  ){
							toPrint=false;
							break;
						}
					}
				}

				if(++counterSR%100000==0){cout<<100*counterSR/lines.size()<<"% done"<<endl;}
				if(toPrint){
					#pragma omp critical(dataupdate)
					{
						for(uint j(0);j<Line.size();++j){
							outputFile<<Line[j]<<";";
						}
						outputFile<<"\n";
					}
				}
			}
		}
	}
	outputFile<<flush;
    cout<<"End"<<endl;
    return 0;
}
