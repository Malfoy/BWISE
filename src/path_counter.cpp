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
    cout<<"./pathCounter  path threshold_unitig outputFile [core to use] [superReads_Threshold] [unitig.fa] [k] [affine]"<<endl;
}


//THIS CAN WORK ON BUCKETS
int main(int argc, char *argv[]) {
    if(argc<4){
        help();
        exit(0);
    }
    cout<<"I am pathcounter"<<endl;

    vector<vector<int64_t>> lines;
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
    vector<vector<uint64_t>> unitigsToReads;
    cout<<"Loading unitigs"<<endl;
    if(unitigFile!=""){
        sizeUnitig.push_back(0);
        ifstream unitigStream(unitigFile);
        while(not unitigStream.eof()){
            getline(unitigStream,useless);
            getline(unitigStream,line);
            int size(0);
            size+=line.size();
            //~ cout<<line.size()<<endl;
            size-=1*(kmerSize-1);
            //~ cout<<size<<endl;
            if(size>0){
                sizeUnitig.push_back((uint)size);
                //~ cin.get();
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
                lines.push_back(coucouch);
            }
        }
    }

    //CLEANING
    cout<<"Cleaning low coverage unitigs"<<endl;
    for(uint64_t i(0);i<lines.size();++i){
        for(uint64_t j(0);j<lines[i].size();++j){
            uNumber=(lines[i][j]);
            if(unitigFile!=""){
                if(afineThreshold!=0){
                    if(count[abs(uNumber)] < (sizeUnitig[abs(uNumber)]/afineThreshold)){
                        lines[i]={};
                    }
                }
                if(count[abs(uNumber)]<(threshold_unitig)){
                    lines[i]={};
                }
            }
        }
        canonicalVector(lines[i]);
    }

    //DEDUPLICATING
    cout<<"Removing Duplicate"<<endl;
    //TODO HUGE MEMORY PROBLEME HERE ?
    sort(lines.begin(),lines.end());
    uint64_t pred(0),counter(1);
    for(uint64_t i(1);i<lines.size();++i){
        if(lines[i].empty()){
            continue;
        }
        if(lines[i]==lines[pred]){
            ++counter;
            lines[i]={};
        }else{
            if(counter<superThreshold and false){
                lines[pred]={};
            }else{
                //~ abundance_sr.push_back(counter);
                lines[pred].push_back(counter);
            }
            pred=i;
            counter=1;
        }
    }
    if(counter<superThreshold){
        lines[pred]={};
    }else{
        lines[pred].push_back(counter);
    }
    counter=1;

    sort(lines.begin(),lines.end());
    lines.erase( unique( lines.begin(), lines.end() ), lines.end() );

	string outputFileName(argv[3]);
    ofstream outputFileSolid,outputFileWeak;
    outputFileSolid.open(outputFileName+"Solid");
    outputFileWeak.open(outputFileName+"Weak");
    for(uint i(0);i<lines.size();++i){
        if(lines[i].size()>0){
			uint32_t abundance(lines[i][lines[i].size()-1]);
			if(abundance<superThreshold){
				outputFileWeak<<""+to_string(abundance)<<"\n";
				for(uint j(0);j<lines[i].size()-1;++j){
					outputFileWeak<<lines[i][j]<<";";
				}
				outputFileWeak<<"\n";
			}else{
				outputFileSolid<<""+to_string(abundance)<<"\n";
				for(uint j(0);j<lines[i].size()-1;++j){
					outputFileSolid<<lines[i][j]<<";";
				}
				outputFileSolid<<"\n";
			}
        }
    }
    return 0;
}


