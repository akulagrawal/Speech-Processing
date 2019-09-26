// LBG.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stdio.h"
#include "iostream"
#include "fstream"
#include "vector"
#include "cmath"
#include "string"
#include "iomanip"

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
#define PI 3.14159265358979323846  /* pi */

#define K 16

/* Read CoolEdit data from file "name" and store it in a vector which is
   returned by the function. Heading comprising of top 5 lines is stored
   in vector "heading" if header = true
   if nSamples < 0, whole file is stored, else only first nSamples lines are read and stored
*/
vector<ld> readFileToVec(string name, bool header = true, int nSamples = -1, vector<string> &heading = vector<string>()){
    ifstream file(name);
    string str;
    int count = 5;
    if(!header)
        count = 0;
    vector<ld> v;
	heading.clear();
    while (std::getline(file, str)){
        count--;
        if(count < 0){
			v.pb(stold(str));
			nSamples--;
			if(nSamples == 0)
				break;
		}
		else
			heading.pb(str);
    }
    return v;
}

/* Store the data in vector v in the file "name"
   SAMPLES attribute in "heading" is updated
*/
void writeVecToFile(vector<ld> v, string name, vector<string> heading = vector<string>()){
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
        myfile<<v[i]<<endl;
    myfile.close();
}

/* Correct the DC Shift in vector v by subtracting average
*/
void correctDCShift(vector<ld> &v){
    ld avg,sum=0;
    for(int i=0;i<v.size();i++)
        sum+=v[i];
    avg=sum/(ld)v.size();
    for(int i=0;i<v.size();i++)
        v[i]-=avg;
}

/* Normalize the data in vector v to contain entries
   from -range to +range
*/
void normalize(vector<ld> &v, ld range = 5000.0){
    ld maxm=v[0];
    for(int i=0;i<v.size();i++)
        maxm=max(maxm,abs(v[i]));
    for(int i=0;i<v.size();i++)
        v[i]*=range/maxm;
}

/* Put the markers in the beginning and end of useful entries
   and delete the useless entries before first marker
   and after the next marker
*/
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

void selectFrames(vector<vector<ld> > &frame, vector<ld> v, int N, int n){
	if(N*n > v.size())
		return;
	vector<ld> temp;
	frame.clear();
	for(int i=0;i<n;i++){
		temp.clear();
		for(int j=0;j<N;j++)
			temp.pb(v[(N*i)+j]);
		frame.pb(temp);
	}
}

/* Function which applies hamming window to first frame of vector v
   consisting of N elements and returns the resulting vector.
*/
void applyHammingWindow(vector<ld> &v, int N){
	ld w;
	//cout<<"Hamming Window: \n";
	for(int i=0;i<N;i++){
		w = 0.54 - 0.46*cos(((ld)2.0*(ld)PI*(ld)i)/(ld)(N-1));
		v[i]*=w;
		//cout<<w<<endl;
	}
}

/* Function to calculate and return R(0), R(1), ..., R(p)
   s: signal vector (after windowing if windowing is used)
   N: number of samples in the frame
   p: number of past samples used in LPC
*/
vector<ld> calculateR(vector<ld> s, int N, int p){
	vector<ld> r;
	for(int k=0;k<=p;k++){
		ld temp = 0;
		for(int m=0;m<N-k;m++)
			temp += s[m]*s[m+k];
		r.pb(temp);
	}
	return r;
}

/* Function that implements Durbin's Algorithm taking vector R(0), R(1), ..., R(p)
   as input and solve the equations to find and return the values of a(1), a(2), ..., a(p)
*/
vector<ld> durbin(vector<ld> r, int p){
	vector<vector<ld> > alpha(p+1, vector<ld>(p+1, 0.0));
	vector<ld> k;
	vector<ld> e;
	vector<ld> a;
	e.pb(r[0]);
	k.pb(0);
	for(int i=1;i<=p;i++){
		ld temp = 0.0;
		for(int j=1;j<i;j++)
			temp += alpha[j][i-1]*r[i-j];
		k.pb((r[i]-temp)/e[i-1]);
		alpha[i][i]=k[i];
		for(int j=1;j<i;j++)
			alpha[j][i] = alpha[j][i-1] - k[i]*alpha[i-j][i-1];
		e.pb((1-k[i]*k[i])*e[i-1]);
	}
	a.pb(0);
	for(int i=1;i<=p;i++)
		a.pb(alpha[i][p]);
	return a;
}

/* Function to calculate Cepstral Coefficients from
   a(1), a(2), ..., a(p) obtained as a solution to optimality equations.
   R(0) is required to approximate the gain term
*/
vector<ld> calculateC(vector<ld> a, vector<ld> r, int p){
	vector<ld> c;
	c.pb(log(r[0]*r[0])/log((ld)2));
	for(int m=1;m<=p;m++){
		c.pb(a[m]);
		for(int k=1;k<m;k++)
			c[m] += (ld)k*c[k]*a[m-k]/(ld)m;
	}
	return c;
}

/* Function to print the final output
*/
void print_output(vector<ld> r, vector<ld> a, vector<ld> c, int p){
	cout<<std::fixed;
    cout<<std::setprecision(6);
	cout<<"Ri: ";
	for(int i=0;i<=p;i++)
		cout<<r[i]<<"  ";
	cout<<endl;
	cout<<"Ai: ";
	for(int i=1;i<=p;i++)
		cout<<a[i]<<"  ";
	cout<<endl;
	cout<<"Ci: ";
	for(int i=1;i<=p;i++)
		cout<<c[i]<<"  ";
	cout<<endl;
}

/* Function which applies raised sine window to first frame of vector v
   consisting of N elements and returns the resulting vector.
*/
void applyRaisedSineWindow(vector<ld> &v, int N){
	ld w;
	//cout<<"Raised Sine Window: \n";
	for(int m=1;m<=N;m++){
		w = 1 + ((ld)N/(ld)2.0)*sin(((ld)PI*(ld)m)/(ld)N);
		v[m-1]*=w;
		//cout<<w<<endl;
	}
}

vector<vector<ld> > calculateCFromScratch(string fileName, int N, int p, int nFrames, bool isHeader){
	vector<ld> v;
	vector<vector<ld> > frame;
	vector<ld> r;
	vector<ld> a;
	vector<vector<ld> > c(nFrames+1, vector<ld>(p+1, 0.0));

	v = readFileToVec(fileName, isHeader);
	correctDCShift(v);
	normalize(v);
	putMarkers(v);
	selectFrames(frame, v, N, nFrames);
	for(int k=0;k<nFrames;k++){
		r = calculateR(frame[k], N, p);
		a = durbin(r, p);
		c[k] = calculateC(a, r, p);
		applyRaisedSineWindow(c[k], p);
	}
	return c;
}

ld calculateTokhuraDistance(vector<ld> v1, vector<ld> v2){
	ld wt[13] = {0.0, 1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
	ld sum = 0;
	for(int i=0;i<13;i++)
		sum += wt[i]*(v1[i]-v2[i])*(v1[i]-v2[i]);
	return sum;
}

void runTests(int N, int p, bool isHeader = true){
	string vowel[5] = {"a", "e", "i", "o", "u"};
	int nFrames = 5;
	int nTests = 10;
	string testName;
	string fileName;

	vector<vector<ld> > c;
	vector<ld> v;
	int correct = 0;
    
    for(int i=0;i<5;i++){
		for(int j=0;j<nTests;j++){
			testName = "./Tests/160101085_" + vowel[i] + "_1" + to_string((ll)j) + ".txt";
			c = calculateCFromScratch(testName, N, p, nFrames, isHeader);
			ld dist = numeric_limits<ld>::max();
			string result;
			for(int k=0;k<5;k++){
				fileName = "./Vowels/c_" + vowel[k] + ".txt";
				v = readFileToVec(fileName, false);
				ld sum = 0;
				for(int n=0;n<nFrames;n++){
					vector<ld> v2 = vector<ld>(v.begin () + n*(p+1), v.begin () + (n+1)*(p+1));
					sum += calculateTokhuraDistance(c[n+1],v);
				}
				if(sum < dist){
					dist = sum;
					result = vowel[k];
				}
			}
			cout<<testName<<endl<<"Predicted: "<<result<<endl;
			if(result == vowel[i])
				correct++;
		}
	}
	cout<<"Accuracy: "<<(ld)correct* 100.0/50.0<<"%"<<endl;
}

void runRealTime(int N, int p, bool isHeader = false){
	const char* exeCommand = "C:/Recording_Module/Recording_Module.exe 3 rec.wav rec.txt";
	system(exeCommand);

	string vowel[5] = {"a", "e", "i", "o", "u"};
	int nFrames = 5;
	string fileName;

	vector<vector<ld> > c;
	vector<ld> v;
    
	c = calculateCFromScratch("rec.txt", N, p, nFrames, isHeader);
	ld dist = numeric_limits<ld>::max();
	string result;
	for(int k=0;k<5;k++){
		fileName = "./Vowels/c_" + vowel[k] + ".txt";
		v = readFileToVec(fileName, false);
		ld sum = 0;
		for(int n=0;n<nFrames;n++){
			vector<ld> v2 = vector<ld>(v.begin () + n*(p+1), v.begin () + (n+1)*(p+1));
			sum += calculateTokhuraDistance(c[n+1],v);
		}
		if(sum < dist){
			dist = sum;
			result = vowel[k];
		}
	}
	cout<<"Predicted: "<<result<<endl;
}

string chooseMode(){
	string mode;
	cout<<"\t\t\tMENU\n";
	cout<<"1. Run pre-recorded Tests\n";
	cout<<"2. Real-Time recording\n";
	while(1){
		cout<<"Enter 1 or 2: ";
		cin>>mode;
		if((mode!="1") && (mode!="2"))
			cout<<"Invalid input, please try again\n";
		else
			break;
	}
	return mode;
}


vector<vector<ld> > centroid(vector<int> idx, vector<vector<vector<ld> > > x) {
	vector<vector<ld> > avg(x[0].size(), vector<ld>(x[0][0].size(), 0.0));
	if(!idx.size())
		return avg;
	for(int i=0;i<idx.size();i++) {
		for(int j=0;j<x[0].size();j++) {
			for(int l=0;l<x[0][0].size();l++) {
				avg[j][l] += x[idx[i]][j][l];
			}
		}
	}
	for(int j=0;j<x[0].size();j++) {
		for(int l=0;l<x[0][0].size();l++) {
			avg[j][l] /= idx.size();
		}
	}
	return avg;
}

vector<vector<vector<ld> > > KMeans(vector<vector<vector<ld> > > x)
{
	vector<vector<int> > bucket(K+1, vector<int>(1, 0));
	bool mark[K+1];
	ld d, dmin, idx;
	vector<vector<vector<ld> > > y(K+1, vector<vector<ld> >(x[0].size(), vector<ld>(x[0][0].size(), 0.0)));
	vector<vector<vector<ld> > > c(K+1, vector<vector<ld> >(x[0].size(), vector<ld>(x[0][0].size(), 0.0)));
	memset(mark, 0, sizeof(mark));
	for(int i=0;i<K;i++) {
		while(1) {
			int idx = rand()%x.size();
			if((!mark[idx]) && (idx > 0)) {
				y[i+1] = x[idx];
				mark[idx] = 1;
				break;
			}
		}
	}

	ld dist = numeric_limits<ld>::max();
	ld new_dist = 0;
	while(1) {
		for(int i=0;i<K+1;i++) {
			bucket[i].clear();
		}
		for(int i=1;i<x.size();i++) {
			ld d, dmin = numeric_limits<ld>::max();
			int idx = 1;
			for(int j=0;j<K;j++) {
				d = 0;
				for(int l=0;l<x[0].size();l++) {
					d += calculateTokhuraDistance(x[i][l], y[j+1][l]);
					if(d < dmin) {
						dmin = d;
						idx = j+1;
					}
				}
			}
			bucket[idx].pb(i);
		}
		
		for(int i=0;i<K;i++) {
			if(bucket[i+1].size())
				c[i+1] = centroid(bucket[i+1], x);
		}

		new_dist = 0;
		for(int i=0;i<K;i++) {
			for(int j=0;j<bucket[i+1].size();j++) {
					int idx = bucket[i+1][j];
					for(int l=0;l<x[0].size();l++) {
						new_dist += calculateTokhuraDistance(x[idx][l], c[i+1][l]);
					}
			}
		}
		cout<<dist<<" "<<new_dist<<endl;
		if(dist - new_dist < 1.0) {
			break;
		}
		else {
			y = c;
			dist = new_dist;
		}
	}
	return y;
}


int _tmain(int argc, _TCHAR* argv[])
{
	//These values need to be set correctly before running the program
	ld epsilon = 10.0;
	bool applyWindow = false;
	bool isHeader = true;
	int N = 320;
	int p = 12;

	string vowel[5] = {"a", "e", "i", "o", "u"};
    int nSamples = 5;
	int nFrames = 5;
	int nTests = 10;
	string mode;

	vector<vector<vector<ld> > > x(5*nSamples+1, vector<vector<ld> >(nFrames+1, vector<ld>(p+1, 0.0)));
	vector<int> cluster[K+1];
	for(int i=0;i<K+1;i++)
		cluster[i].clear();
    
    for(int i=0;i<5;i++){
		for(int j=1;j<=nSamples;j++){
			string fileName = "./Vowels/160101085_" + vowel[i] + "_" + to_string((ll)j) + ".txt";
			x[j+(i*nSamples)] = calculateCFromScratch(fileName, N, p, nFrames, isHeader);
		}
	}

	//START LBG CLUSTERING
	vector<vector<int> > bucket(K+1, vector<int>(1, 0));
	ld d, dmin, idx;
	vector<vector<vector<ld> > > y(K+1, vector<vector<ld> >(x[0].size(), vector<ld>(x[0][0].size(), 0.0)));
	vector<vector<vector<ld> > > c(K+1, vector<vector<ld> >(x[0].size(), vector<ld>(x[0][0].size(), 0.0)));
	while(1) {
		int idx = rand()%x.size();
		if(idx > 0) {
			y[1] = x[idx];
			break;
		}
	}

	int k=1;

	while(k <= K/2) {

		ld dist = numeric_limits<ld>::max();
		ld new_dist = 0;
		while(1) {
			for(int i=0;i<k+1;i++) {
				bucket[i].clear();
				cluster[i].clear();
			}
			for(int i=1;i<x.size();i++) {
				ld d, dmin = numeric_limits<ld>::max();
				int idx = 1;
				for(int j=0;j<k;j++) {
					d = 0;
					for(int l=0;l<x[0].size();l++) {
						d += calculateTokhuraDistance(x[i][l], y[j+1][l]);
						if(d < dmin) {
							dmin = d;
							idx = j+1;
						}
					}
				}
				bucket[idx].pb(i);
			}
			
			for(int i=0;i<k;i++) {
				if(bucket[i+1].size())
					c[i+1] = centroid(bucket[i+1], x);
			}

			new_dist = 0;
			for(int i=0;i<k;i++) {
				for(int j=0;j<bucket[i+1].size();j++) {
					int idx = bucket[i+1][j];
					for(int l=0;l<x[0].size();l++) {
						new_dist += calculateTokhuraDistance(x[idx][l], c[i+1][l]);
					}
				}
			}
			cout<<dist<<" "<<new_dist<<endl;
			if(dist - new_dist <= 0.0) {
				break;
			}
			else {
				y = c;
				dist = new_dist;
			}
		}

		for(int i=1;i<x.size();i++) {
			ld d, dmin = numeric_limits<ld>::max();
			int idx = 1;
			for(int j=0;j<k;j++) {
				d = 0;
				for(int l=0;l<x[0].size();l++) {
					d += calculateTokhuraDistance(x[i][l], y[j+1][l]);
					if(d < dmin) {
						dmin = d;
						idx = j+1;
					}
				}
			}
			cluster[idx].pb((i-1)/5);
		}

		for(int i=1;i<k+1;i++) {
			cout<<i<<": ";
			for(int j=0;j<cluster[i].size();j++) {
				cout<<vowel[cluster[i][j]]<<" ";
			}
			cout<<endl;
		}
		for(int i=1;i<=k;i++) {
			y[i+k] = y[i];
			for(int j=0;j<nFrames;j++) {
				for(int l=0;l<p;l++) {
					y[i][j][l] -= epsilon;
					y[i][j][l] += epsilon;
				}
			}
		}
		k*=2;
	}


	getchar();
	return 0;
}