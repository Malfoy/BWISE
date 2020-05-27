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



using namespace std;




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

size_t hashV(const vector<int64_t>& v){
	uint64_t res(0);
	for(uint i(0);i<v.size() and i<2;++i){
		res^= hash64shift(v[i]);
	}
    return res;
}


struct MyHashV{
  size_t operator()(vector<int64_t> m) const {
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
    cout<<"./pathCounter  path threshold_unitig outputFile [core to use] [superReads_Threshold] [unitig.fa] [k] [affine]"<<endl;
}


//THIS CAN WORK ON BUCKETS
int main(int argc, char *argv[]) {
    if(argc<4){
        help();
        exit(0);
    }
    cout<<"I am pathcounter"<<endl;

    // vector<vector<int64_t>> lines;
	robin_hood::unordered_node_map<vector<int64_t>,uint32_t,MyHashV> lines;
    vector<int64_t> coucouch;
    vector<uint> sizeUnitig;
    string seqFile(argv[1]),unitigFile;
    uint threshold_unitig(stoi(argv[2])),superThreshold(0),kmerSize,afineThreshold(0),coreUsed(4);
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
    if(argc>8){
        afineThreshold=(stoi(argv[8]));
    }


    int64_t uNumber(0);
    string line,useless,msp,number;
    auto numStream=new zstr::ifstream(seqFile);
    vector<uint64_t> count,abundance_sr;
    cout<<"Loading unitigs"<<endl;
    if(unitigFile!=""){
        sizeUnitig.push_back(0);
        ifstream unitigStream(unitigFile);
        while(not unitigStream.eof()){
            getline(unitigStream,useless);
            getline(unitigStream,line);
            int size(0);
            size+=line.size();
            size-=1*(kmerSize-1);
            if(size>0){
                sizeUnitig.push_back((uint)size);
            }
        }
    }

    count.resize(sizeUnitig.size(),0);

    //LOADING and Counting
    while(not numStream->eof()){
        getline(*numStream,line);
        coucouch={};
        if(line.size()>1){
            uint64_t i(1),lasti(0);
            while(i<line.size()){
                if(line[i]==';'){
                    number=line.substr(lasti,i-lasti);
                    lasti=i+1;
                    uNumber=stoi(number);
                    coucouch.push_back(uNumber);
                    count[abs(uNumber)]++;
                }
                ++i;
            }
            if(coucouch.size()!=0){
                // lines.push_back(coucouch);
				canonicalVector(coucouch);
				lines[coucouch]++;
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
				for(uint j(0);j<element.first.size();++j){
					outputFileWeak<<element.first[j]<<";";
				}
				outputFileWeak<<"\n";
			}else{
				outputFileSolid<<""+to_string(abundance)<<"\n";
				for(uint j(0);j<element.first.size();++j){
					outputFileSolid<<element.first[j]<<";";
				}
				outputFileSolid<<"\n";
			}
        }
    }
	outputFileSolid.flush();
	outputFileWeak.flush();
	cout<<"singleton"<<singleton<<"doubleton"<<doubleton<<endl;
    return 0;
}
