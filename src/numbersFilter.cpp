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


vector<int> reverseVector(const vector<int>& V){
	vector<int> RC;
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



bool isPrefix(const vector<int>& v1, const vector<int>& v2){
	if(v2.size()<v1.size()){return false;}
	for(uint i(0);i<v1.size();++i){
		if(v1[i]!=v2[i]){
			return false;
		}
	}
	return true;
}


bool isSuffix(const vector<int>& v1, const vector<int>& v2){
	if(v2.size()<v1.size()){return false;}
	for(uint i(0);i<v1.size();++i){
		if(v1[i]!=v2[i+v2.size()-v1.size()]){
			return false;
		}
	}
	return true;
}


bool isInclued(const vector<int>& v1, const vector<int>& v2){
	if(v2.size()<v1.size()){return false;}
	for(uint j(0);j<v2.size();++j){
		bool success(true);
		for(uint i(0);i<v1.size() and i+j<v2.size();++i){
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
	cout<<"./numbersFilter  numbers.txt threshold outputFile [superReadsThreshold] [unitig.fa] [k] [affine]"<<endl;
}


//THIS CAN WORK ON BUCKETS
int main(int argc, char *argv[]) {
	if(argc<4){
		help();
		exit(0);
	}

	vector<vector<int>> lines;
	vector<int> coucouch;
	vector<uint> sizeUnitig;
	string seqFile(argv[1]),unitigFile;
	uint threshold(stoi(argv[2])),superThreshold(0),kmerSize,afineThreshold(20);
	bool headerNeed(false);
	if(argc>4){
		superThreshold=(stoi(argv[4]));
	}
	if(argc>6){
		unitigFile=((argv[5]));
		kmerSize=(stoi(argv[6]));
	}
	//~ if(argc>7){
		//~ headerNeed=true;
	//~ }
	if(argc>7){
		afineThreshold=(stoi(argv[7]));
	}
	int uNumber;
	string line,useless,msp,number;
	ifstream numStream(seqFile);
	vector<uint> count;
	vector<vector<uint>> unitigsToReads;
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
	sort(lines.begin(),lines.end(),[](const vector<int> a, const vector<int> b) {return a < b;});
	uint pred(0),counter(1);
	for(uint i(1);i<lines.size();++i){
		if(lines[i]==lines[pred]){
			++counter;
			lines[i]={};
		}else{
			if(counter<superThreshold){
				lines[pred]={};
			}else{
				if(isPrefix(lines[pred],lines[i])){
					lines[pred]={};
				}
			}
			pred=i;
			counter=1;
		}
	}
	sort(lines.begin(),lines.end(), [](const vector<int> a, const vector<int> b) {return a < b;});
	lines.erase( unique( lines.begin(), lines.end() ), lines.end() );
	unitigsToReads.resize(count.size()+1,{});
	//FILLING
	for(uint i(0);i<lines.size();++i){
		for(uint j(0);j<lines[i].size();++j){
			unitigsToReads[abs(lines[i][j])].push_back(i);
			//~ cout<<abs(lines[i][j])<<" "<<i<<endl;
		}
	}
	for(uint i(0);i<lines.size();++i){
		//~ cout<<"go"<<i<<endl;
		unordered_map<uint,uint> readScore;
		for(uint j(0);j<lines[i].size();++j){
			//~ cout<<"GO"<<j<<endl;
			//~ cout<<abs(lines[i][j])<<endl;
			for(uint ii(0); ii<unitigsToReads[abs(lines[i][j])].size(); ++ii){
				uint friendRead=unitigsToReads[abs(lines[i][j])][ii];
				if(friendRead!=i){
					++readScore[friendRead];
				}
				//~ cout<<readScore[unitigsToReads[abs(lines[i][j])][ii]]<<endl;
				if(readScore[friendRead]>=lines[i].size()){
					if( isInclued( lines[i], lines[friendRead] ) or  isInclued( reverseVector(lines[i]), lines[friendRead] )  ){
						lines[i]={};
					}
					//~ cout<<"DESTC56TG"<<endl;
				}
			}
		}
	}
    ofstream outputFile;
    outputFile.open(argv[3]);
	//OUTPUT
	counter=(0);
	for(uint i(0);i<lines.size();++i){
		if(lines[i].size()>=1){
			if(headerNeed){
				outputFile<<">"+to_string(counter++)<<endl;
			}
			for(uint j(0);j<lines[i].size();++j){
				outputFile<<lines[i][j]<<";";
			}
			outputFile<<endl;
		}
	}
    outputFile.close();
    return 0;
}
