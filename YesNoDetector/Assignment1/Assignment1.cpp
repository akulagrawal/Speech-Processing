// Assignment1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stdio.h"
#include "iostream"
#include "fstream"
#include "string"
#include "vector"

using namespace std;
#define f first
#define s second
#define ll long long
#define ld long double
#define mp make_pair
#define MAX 1000006
#define mod 1000000007
#define pb push_back
#define INF 1e18
#define pii pair<int,int>

// Read CoolEdit data from file "name" and store it in a vector which is
// returned by the function. Heading comprising of top 5 lines is stored
// in vector "heading" if header = true
vector<ld> readFileToVec(string name, vector<string> &heading, bool header = true){
    ifstream file(name);
    string str;
    int count = 5;
    if(!header)
        count = 0;
    vector<ld> v;
	heading.clear();
    while (std::getline(file, str)){
        count--;
        if(count < 0)
            v.pb(stoi(str));
		else
			heading.pb(str);
    }
    return v;
}

// Store the data in vector v in the file "name"
// SAMPLES attribute in "heading" is updated
void writeVecToFile(vector<ld> v, string name, vector<string> heading){
	if(heading.size()){
		while(isdigit(heading[0].back()))
			heading[0].pop_back();
		heading[0] += to_string((ll)v.size());
	}
    ofstream myfile;
    myfile.open (name);
	for(int i=0;i<heading.size();i++)
		myfile<<heading[i]<<endl;
    for(int i=0;i<v.size();i++)
        myfile<<(ll)v[i]<<endl;
    myfile.close();
}

// Correct the DC Shift in vector v by subtracting average
void correctDCShift(vector<ld> &v){
    ld avg,sum=0;
    for(int i=0;i<v.size();i++)
        sum+=v[i];
    avg=sum/(ld)v.size();
    for(int i=0;i<v.size();i++)
        v[i]-=avg;
}

// Normalize the data in vector v to contain entries
// from -range to +range
void normalize(vector<ld> &v, ld range = 5000.0){
    ld maxm=v[0];
    for(int i=0;i<v.size();i++)
        maxm=max(maxm,abs(v[i]));
    for(int i=0;i<v.size();i++)
        v[i]*=range/maxm;
}

// Put the markers in the beginning and end of useful entries
// and delete the useless entries before first marker
// and after the next marker
void putMarkers(vector<ld> &v){
    ld maxm=v[0];
    for(int i=0;i<v.size();i++)
        maxm=max(maxm,abs(v[i]));
    maxm/=10.0;
    while(v[0] < maxm)
        v.erase(v.begin());
    while(v.back() < maxm)
        v.pop_back();
}

// Return vector of framewise energy for signal in vector v
// by dividing v into frames of size fsize and dropping
// the last frame
vector<ld> calculateFramewiseE(vector<ld> v, int fsize = 100){
    vector<ld> e;
    ld temp;
    for(int i=0;i<v.size()-fsize;i+=fsize){
        temp = 0.0;
        for(int j=0;j<fsize;j++)
            temp+=(v[i+j]*v[i+j]);
        temp/=(ld)fsize;
        e.pb(temp);
    }
    return e;
}

// Return vector of framewise ZCR for signal in vector v
// by dividing v into frames of size fsize and dropping
// the last frame
vector<int> calculateFramewiseZCR(vector<ld> v, int fsize = 100){
    vector<int> zcr;
    int temp;
    for(int i=0;i<v.size()-fsize;i+=fsize){
        temp = 0;
        for(int j=1;j<fsize;j++)
            if(v[i+j]*v[i+j-1] < 0)
                temp++;
        zcr.pb(temp);
    }
    return zcr;
}

int _tmain(int argc, _TCHAR* argv[])
{
	// These parameters have to be set correctly before code execution
	string keywords[2] = {"yes", "no"};
    int nSamples = 10;
	string dir = "C:/Users/Sanjai Kumar Agrawal/Desktop/CS566/Data/";
    string name;
    vector<ld> v;
	vector<string> heading;
    
    for(int i=0;i<2;i++){
        for(int j=1;j<=nSamples;j++){
            name = keywords[i] + "_" + to_string((ll)j);
            cout<<name<<endl;
            v = readFileToVec(dir + name + ".txt", heading);
            correctDCShift(v);
            normalize(v);
            putMarkers(v);
			writeVecToFile(v, dir + name + "_seg.txt", heading);
			v = readFileToVec(dir + name + "_seg.txt", heading);
            vector<ld> e;
            vector<int> zcr;
            e = calculateFramewiseE(v);
            zcr = calculateFramewiseZCR(v);
            //for(int i=e.size()-8;i<e.size();i++)
            //    cout<<e[i]<<" ";
            //cout<<endl;
			double zcrAvg = 0;
            for(int i=zcr.size()-6;i<zcr.size();i++)
                zcrAvg += (double)zcr[i];
			zcrAvg /= 6.0;
			//cout<<zcrAvg;
			//cout<<endl;
			cout<<"Predicted: ";
			if(zcrAvg < 7.0)
				cout<<"No\n";
			else
				cout<<"Yes\n";
        }
    }
	getchar();
	return 0;
}

