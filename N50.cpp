#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>



using namespace std;



int main(int argc, char ** argv){
	if(argc<2){
		cout<<"[Fasta file]"<<endl;
		exit(0);
	}
	string input(argv[1]);
	srand (time(NULL));
	string ref, useless;
	ifstream in(input);
	vector<uint> lengths;
	uint size(0);
	while(not in.eof()){
		getline(in,useless);
		getline(in,ref);
		if(not ref.empty() and not useless.empty()){
			lengths.push_back(ref.size());
			size+=ref.size();
		}
	}
	sort(lengths.begin(),lengths.end(),greater<uint>());
	uint total(0),i(0);
	while(total<size*0.5){
		total+=lengths[i];
		++i;
	}
	cout<<"N50: "<<lengths[i-1]<<endl;
	cout<<"L50: "<<i<<endl;
	total=i=0;
	while(total<size*0.80){
		total+=lengths[i];
		++i;
	}
	cout<<"N80: "<<lengths[i-1]<<endl;
	cout<<"L80: "<<i<<endl;
	total=i=0;
	while(total<size*0.90){
		total+=lengths[i];
		++i;
	}
	cout<<"N90: "<<lengths[i-1]<<endl;
	cout<<"L90: "<<i<<endl;
}
