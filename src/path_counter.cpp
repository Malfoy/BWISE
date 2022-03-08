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
#include "robin_hood.h"
#include "bf.h"



using namespace std;


typedef  int32_t UN;


int64_t hash64shift(int64_t key) {
	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return key;
}



size_t hashV(const vector<UN>& v){
	uint64_t res(0);
	for(uint i(0);i<v.size();++i){
		res^= hash64shift(v[i]);
	}
    return res;
}



struct MyHashV{
  size_t operator()(vector<UN> m) const {
    return hashV(m);
  }
};


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



void canonicalVector(vector<UN>& V){
    vector<UN> RC;
    for(uint i(0);i<V.size();++i){
        RC.push_back(-V[V.size()-i-1]);
    }
    if(V>RC){
        V=RC;
    }
}



vector<UN> reverseVector(const vector<UN>& V){
    vector<UN> RC;
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



bool isPrefix(const vector<UN>& v1, const vector<UN>& v2){
    if(v2.size()<v1.size()){return false;}
    for(uint i(0);i<v1.size();++i){
        if(v1[i]!=v2[i]){
            return false;
        }
    }
    return true;
}


bool isSuffix(const vector<UN>& v1, const vector<UN>& v2){
    if(v2.size()<v1.size()){return false;}
    for(uint64_t i(0);i<v1.size();++i){
        if(v1[i]!=v2[i+v2.size()-v1.size()]){
            return false;
        }
    }
    return true;
}


bool isInclued(const vector<UN>& v1, const vector<UN>& v2){
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


string get_canon_string(const string& line){
	string result;
	vector<UN> coucouch;
	if(line.size()>1){
		uint64_t i(1),lasti(0);
		string number;
		while(i<line.size()){
			if(line[i]==';'){
				number=line.substr(lasti,i-lasti);
				lasti=i+1;
				UN uNumber=stoi(number);
				coucouch.push_back(uNumber);
			}
			++i;
		}
		if(coucouch.size()!=0){
			canonicalVector(coucouch);
			for(uint i(0);i<coucouch.size();++i){
				result+=to_string(coucouch[i])+";";
			}
		}
	}

	return result;
}



void help(){
    cout<<"./pathCounter  path threshold_unitig outputFile [core to use] [superReads_Threshold] [unitig.fa] [k] [affine]"<<endl;
}


//THIS CAN WORK ON BUCKETS
int main(int argc, char *argv[]) {
    if(argc<4){
        help();
        exit(0);
    }
    cout<<"I am pathcounter"<<endl;
	bloom_parameters parameters;
	parameters.projected_element_count = 3000000000;
	parameters.false_positive_probability = 0.01;
	parameters.random_seed = 0xA5A5A5A5;
	parameters.maximum_number_of_hashes = 3;
	parameters.compute_optimal_parameters();
	//THIS GIVE ROUHLY THAN 4GB
	bloom_filter filter(parameters);

    // vector<vector<int64_t>> lines;
	robin_hood::unordered_node_map<string,uint32_t> lines;
    vector<UN> coucouch;
    string seqFile(argv[1]),unitigFile;
	// TODO COUL REUSE threshold_unitig
    uint threshold_unitig(stoi(argv[2])),superThreshold(0),kmerSize,coreUsed(8);
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



    string line,useless,msp,number;
    auto numStream=new zstr::ifstream(seqFile);

    


    //LOADING and Counting
	#pragma omp parallel num_threads(coreUsed)
	{
		string cstr,linep;
	    while(not numStream->eof()){
            linep.clear();
			#pragma omp critical (file)
			{
				getline(*numStream,linep);
			}
			cstr=(get_canon_string(linep));
			if(not cstr.empty()){
				if(filter.contains(cstr)){
					#pragma omp critical (hash)
					{
						lines[cstr]++;
					}
				}else{
					#pragma omp critical (bloom)
					{
						filter.insert(cstr);
					}
				}
			}
	    }
	}

	string outputFileName(argv[3]);
    zstr::ofstream outputFileSolid(outputFileName+"Solid"),outputFileWeak(outputFileName+"Weak");
	uint64_t singleton(0),doubleton(0);

    // for(uint i(0);i<lines.size();++i){
	for (auto element : lines){
        if(element.first.size()>0){
			uint32_t abundance(element.second);
			if(abundance<superThreshold){
				if(abundance==1){
					singleton++;
				}
				if(abundance==2){
					doubleton++;
				}
				outputFileWeak<<""+to_string(abundance)<<"\n";
				outputFileWeak<<element.first<<"\n";
				outputFileWeak<<"\n";
			}else{
				outputFileSolid<<""+to_string(abundance)<<"\n";
				outputFileSolid<<element.first<<"\n";
			}
        }
    }
	outputFileSolid.flush();
	outputFileWeak.flush();
	cout<<"singleton"<<singleton<<"doubleton"<<doubleton<<endl;
    return 0;
}
