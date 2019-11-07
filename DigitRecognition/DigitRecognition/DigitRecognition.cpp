// DigitRecognition.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<iomanip>
#include<algorithm>
#include<math.h>

#define ll long long
#define ld long double
#define pb push_back
#define THRESHOLD pow((ld)10,-30)
#define PI 3.14159265358979323846  /* pi */

using namespace std;

int N=5, M=32, T=500;
ld p, pBar, pStar;

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

void selectFrames(vector<vector<ld> > &frame, vector<ld> v, int lenFrame, int window){
	vector<ld> temp;
	frame.clear();
	int i=0;
	while(1){
		temp.clear();
		if((i+lenFrame-1)>=v.size())
			break;
		for(int j=0;j<lenFrame;j++)
			temp.pb(v[i+j]);
		frame.pb(temp);
		i+=window;
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

vector<vector<ld> > calculateCFromScratch(string fileName, int lenFrame, int p, int window, bool isHeader){
	vector<ld> v;
	vector<vector<ld> > frame;
	vector<ld> r;
	vector<ld> a;

	v = readFileToVec(fileName, isHeader);
	correctDCShift(v);
	normalize(v);
	putMarkers(v);
	selectFrames(frame, v, lenFrame, window);
	int nFrames = frame.size();
	vector<vector<ld> > c(nFrames+1, vector<ld>(p+1, 0.0));
	for(int k=0;k<nFrames;k++){
		r = calculateR(frame[k], lenFrame, p);
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


vector<ld> centroid(vector<int> idx, vector<vector<ld> > x) {
	vector<ld> avg(x[0].size(), 0.0);
	if(!idx.size())
		return avg;
	for(int i=0;i<idx.size();i++) {
		for(int j=0;j<x[0].size();j++) {
			avg[j] += x[idx[i]][j];
		}
	}
	for(int j=0;j<x[0].size();j++) {
		avg[j] /= idx.size();
	}
	return avg;
}

void fillEmptyBuckets(vector<vector<int> > &bucket, int K) {

	for(int i=0;i<K;i++) {
		if(bucket[i].size() == 0) {
			int maxm=0;
			int idx=0;
			for(int j=0;j<bucket.size();j++) {
				if(bucket[j].size() > maxm) {
					maxm = bucket[j].size();
					idx = j;
				}
			}
			maxm /= 2;
			while(maxm--) {
				bucket[i].pb(bucket[idx][0]);
				bucket[idx].erase(bucket[idx].begin());
			}
		}
	}
}

vector<vector<ld> > LBG(vector<vector<ld> > x, int K, int p, int epsilon, ld minSeparation) {
	int nSamples = x.size();
	ld d, dmin;
	int k = 1;
	vector<vector<int> > bucket(K, vector<int>(1, 0));
	vector<vector<ld> > y(K, vector<ld>(p, 0.0));
	vector<vector<ld> > c(K, vector<ld>(p, 0.0));
	int idx = rand()%x.size();
	y[0] = x[idx];

	ld dist = numeric_limits<ld>::max();
	ld new_dist = 0;
	while(k <= K) {
		dist = numeric_limits<ld>::max();
		new_dist = 0;
		while(1) {
			for(int i=0;i<k;i++) {
				bucket[i].clear();
			}
			for(int i=0;i<nSamples;i++) {
				ld d, dmin = numeric_limits<ld>::max();
				int idx = 0;
				for(int j=0;j<k;j++) {
					d = calculateTokhuraDistance(x[i], y[j]);
					if(d < dmin) {
						dmin = d;
						idx = j;
					}
				}
				bucket[idx].pb(i);
			}
			fillEmptyBuckets(bucket, k);
		
			for(int i=0;i<k;i++) {
				c[i] = centroid(bucket[i], x);
			}

			new_dist = 0;
			for(int i=0;i<k;i++) {
				for(int j=0;j<bucket[i].size();j++) {
						int idx = bucket[i][j];
						new_dist += calculateTokhuraDistance(x[idx], c[i]);
				}
			}
			cout<<dist<<" "<<new_dist<<endl;
			if(dist - new_dist <= minSeparation) {
				break;
			}
			else {
				y = c;
				dist = new_dist;
			}
		}
		if(k == K)
			break;
		for(int i=0;i<k;i++) {
			y[i+k] = y[i];
			for(int j=0;j<p;j++) {
				y[i][j] = y[i][j] + epsilon;
				y[i+k][j] = y[i+k][j] - epsilon;
			}
		}
		k = k*2;
	}
	return y;
}


vector<ld> readDim1(string filename){
    ifstream file(filename);
    string str;
    vector<ld> v;
	v.pb(0.0);
    while (file >> str)
      v.pb(stold(str));
    return v;
}

vector<vector<ld> > readDim2(string filename, int M){
    ifstream file(filename);
    string str;
    vector<vector<ld> > v;
	vector<ld> temp;
	v.pb(temp);
	temp.pb(0.0);
	int i=0;
    while (file >> str){
      temp.pb(stold(str));
	  i=(i+1)%M;
	  if(i==0){
		  v.pb(temp);
		  temp.clear();
		  temp.pb(0.0);
	  }
	}
    return v;
}

vector<int> readO(string filename){
    ifstream file(filename);
    string str;
    vector<int> v;
	v.pb(0);
    while (file >> str)
      v.pb(stoi(str));
    return v;
}

void print(vector<ld> pi, vector<vector<ld> > a, vector<vector<ld> > b, vector<int> o){
	cout<<"\nPI:\n";
	for(int i=1;i<=N;i++)
		cout<<'|'<<setw(10)<<pi[i];
	cout<<'|'<<endl;
	cout<<"\nA:\n";
	for(int i=1;i<=N;i++){
		for(int j=1;j<=N;j++)
			cout<<'|'<<setw(10)<<a[i][j];
		cout<<'|'<<endl;
	}
	cout<<"\nB:\n";
	for(int i=1;i<=N;i++){
		for(int j=1;j<=M;j++)
			cout<<'|'<<setw(10)<<b[i][j];
		cout<<'|'<<endl;
	}
	cout<<"\nO:\n";
	for(int i=1;i<=T;i++)
		cout<<'|'<<setw(7)<<o[i];
	cout<<'|'<<endl;
	cout<<endl;
}

// Store the data in vector v in the file "name"
// SAMPLES attribute in "heading" is updated
void writeOutputToFile(vector<ld> pi, vector<vector<ld> > a, vector<vector<ld> > b, int digit){
    ofstream myfile;
    myfile.open ("model_" + to_string((ll)digit) + ".txt");
	myfile<<"\nPI:\n";
	for(int i=1;i<=N;i++)
		myfile<<'|'<<setw(10)<<pi[i];
	myfile<<'|'<<endl;
	myfile<<"\nA:\n";
	for(int i=1;i<=N;i++){
		for(int j=1;j<=N;j++)
			myfile<<'|'<<setw(10)<<a[i][j];
		myfile<<'|'<<endl;
	}
	myfile<<"\nB:\n";
	for(int i=1;i<=N;i++){
		for(int j=1;j<=M;j++)
			myfile<<'|'<<setw(10)<<b[i][j];
		myfile<<'|'<<endl;
	}
	myfile<<endl;
    myfile.close();
}

void solution1(vector<ld> pi, vector<vector<ld> > a, vector<vector<ld> > b, vector<int> o,
			   vector<vector<ld> > &alpha, vector<vector<ld> > &beta){
	//FORWARD PROCEDURE
	for(int i=1;i<=N;i++)
			alpha[1][i] = pi[i]*b[i][o[1]];
	for(int t=1;t<=T-1;t++){
		for(int j=1;j<=N;j++){
			ld sum=0.0;
			for(int i=1;i<=N;i++)
				sum+=alpha[t][i]*a[i][j];
			alpha[t+1][j]=sum*b[j][o[t+1]];
		}
	}
	pBar = 0.0;
	for(int i=1;i<=N;i++)
		pBar+=alpha[T][i];
	//BACKWARD PROCEDURE
	for(int i=1;i<=N;i++)
		beta[T][i] = 1.0;
	for(int t=T-1;t>=1;t--){
		for(int i=1;i<=N;i++){
			ld sum=0.0;
			for(int j=1;j<=N;j++)
				sum+=a[i][j]*b[j][o[t+1]]*beta[t+1][j];
			beta[t][i]=sum;
		}
	}
}

void solution2(vector<ld> pi, vector<vector<ld> > a, vector<vector<ld> > b, vector<int> o,
			   vector<int> &qStar){

	vector<vector<int> > psi(T+1, vector<int>(N+1, 0));
	vector<vector<ld> > delta(T+1, vector<ld>(N+1, 0.0));

	for(int i=1;i<=N;i++){
		delta[1][i] = pi[i]*b[i][o[1]];
		psi[1][i] = 0;
	}
	for(int t=2;t<=T;t++){
		for(int j=1;j<=N;j++){
			ld maxm=0.0;
			int idx=0;
			for(int i=1;i<=N;i++){
				if(maxm<(delta[t-1][i]*a[i][j])){
					maxm=delta[t-1][i]*a[i][j];
					idx=i;
				}
			}
			delta[t][j]=maxm*b[j][o[t]];
			psi[t][j]=idx;
		}
	}
	pStar=*max_element(delta[T].begin()+1,delta[T].end());

	qStar[T]=max_element(delta[T].begin()+1,delta[T].end()) - delta[T].begin();
	for(int t=T-1;t>=1;t--){
		qStar[t]=psi[t+1][qStar[t+1]];
	}
}

void solution3(vector<ld> pi, vector<vector<ld> > a, vector<vector<ld> > b, vector<int> o,
			   vector<ld> &piBar, vector<vector<ld> > &aBar, vector<vector<ld> > &bBar,
			   vector<vector<ld> > alpha, vector<vector<ld> > beta){

	vector<vector<vector<ld> > > xi(T+1, vector<vector<ld> >(N+1, vector<ld>(N+1, 0.0)));
	vector<vector<ld> > gamma(T+1, vector<ld>(N+1, 0.0));

	for(int t=1;t<=T-1;t++){
		ld sum=0.0;
		for(int i=1;i<=N;i++)
			for(int j=1;j<=N;j++)
				sum+=alpha[t][i]*a[i][j]*b[j][o[t+1]]*beta[t+1][j];
		for(int i=1;i<=N;i++)
			for(int j=1;j<=N;j++)
				xi[t][i][j]=(alpha[t][i]*a[i][j]*b[j][o[t+1]]*beta[t+1][j])/sum;
	}
	for(int t=1;t<=T;t++){
		ld sum=0.0;
		for(int i=1;i<=N;i++)
			sum+=alpha[t][i]*beta[t][i];
		for(int i=1;i<=N;i++)
			gamma[t][i]=(alpha[t][i]*beta[t][i])/sum;
	}
	for(int i=1;i<=N;i++)
		piBar[i]=gamma[1][i];
	for(int i=1;i<=N;i++){
		for(int j=1;j<=N;j++){
			ld Num=0.0;
			ld Den=0.0;
			for(int t=1;t<=T-1;t++){
				Num+=xi[t][i][j];
				Den+=gamma[t][i];
			}
			aBar[i][j]=Num/Den;
		}
	}
	for(int j=1;j<=N;j++){
		for(int k=1;k<=M;k++){
			ld Num=0.0;
			ld Den=0.0;
			for(int t=1;t<=T;t++){
				if(o[t]==k)
					Num+=gamma[t][j];
				Den+=gamma[t][j];
			}
			bBar[j][k]=Num/Den;
		}
	}
}

void manipulateB(vector<vector<ld> > &b){
	for(int i=1;i<=N;i++){
		for(int j=1;j<=M;j++){
			if(b[i][j] < THRESHOLD){
				ld dif = THRESHOLD - b[i][j];
				b[i][j] = THRESHOLD;
				ld maxm=0.0;
				int idx=1;
				int k;
				for(k=1;k<=M;k++){
					if(b[i][k] > maxm){
						maxm=b[i][j];
						idx=k;
					}
				}
				b[i][idx] -= dif;
			}
		}
	}
}

vector<int> makeObservationSeq(vector<vector<ld> > c, vector<vector<ld> > codebook, int cbSize){
	int nFrames=c.size();
	vector<int> obs;
	obs.pb(0);
	for(int i=0;i<nFrames;i++){
		int idx=1;
		ld minm=numeric_limits<ld>::max();
		for(int j=1;j<=cbSize;j++){
			ld dist=calculateTokhuraDistance(c[i],codebook[j]);
			if(dist<minm){
				idx=j;
				minm=dist;
			}
		}
		obs.pb(idx);
	}
	return obs;
}

void findAvg(vector<vector<ld> > pi, vector<vector<vector<ld> > > a, vector<vector<vector<ld> > > b, int nTrain,
			 vector<ld> &piAvg, vector<vector<ld> > &aAvg, vector<vector<ld> > &bAvg){
	
	for(int i=1;i<=N;i++){
		ld sum=0.0;
		for(int j=1;j<=nTrain;j++)
			sum+=pi[j][i];
		piAvg[i]=sum/nTrain;
	}
	for(int i=1;i<=N;i++){
		for(int j=1;j<=N;j++){
			ld sum=0.0;
			for(int k=1;k<=nTrain;k++)
				sum+=a[k][i][j];
			aAvg[i][j]=sum/nTrain;
		}
	}
	for(int i=1;i<=N;i++){
		for(int j=1;j<=M;j++){
			ld sum=0.0;
			for(int k=1;k<=nTrain;k++)
				sum+=b[k][i][j];
			bAvg[i][j]=sum/nTrain;
		}
	}
}

string chooseMode(){
	string mode;
	cout<<"\t\t\tMENU\n";
	cout<<"1. Run pre-recorded Tests\n";
	cout<<"2. Real-Time recording\n";
	cout<<"Enter 1 or 2: ";
	cin>>mode;
	return mode;
}

void runTests(vector<vector<ld> > piAvg, vector<vector<vector<ld> > > aAvg, vector<vector<vector<ld> > > bAvg){
	T = 500;
	vector<int> o;
	vector<vector<ld> > alpha(T+1, vector<ld>(N+1, 0.0));
	vector<vector<ld> > beta(T+1, vector<ld>(N+1, 0.0));
	int nDigits = 10;
	int nTest = 10;
	int nTrain = 20;
	int C_N=320;
	int C_p=12;
	int C_window=80;
	int C_cbSize=32;
	bool isHeader = true;

	int correct=0;
	for(int digit=0;digit<nDigits;digit++){
		for(int i=1;i<=nTest;i++){
			string testName = "./Digits_2/160101085_" + to_string((ll)digit) + "_" + to_string((ll)(nTrain+i))+".txt";
			vector<vector<ld> > c = calculateCFromScratch(testName, C_N, C_p, C_window, isHeader);
			vector<vector<ld> > codebook = readDim2("codebook.txt", C_p);
			o = makeObservationSeq(c,codebook,C_cbSize);
			T = o.size()-1;
			int num=0;
			ld prob=numeric_limits<ld>::min();
			for(int dig=0;dig<nDigits;dig++){
				solution1(piAvg[dig], aAvg[dig], bAvg[dig], o, alpha, beta);
				if(pBar >= prob){
					prob=pBar;
					num=dig;
				}
			}
			cout<<testName<<": "<<num<<" "<<prob<<endl;
			if(num==digit)
				correct++;
		}
	}
	cout<<"Accuracy: "<<correct<<"%"<<endl;
}

void runRealTime(vector<vector<ld> > piAvg, vector<vector<vector<ld> > > aAvg, vector<vector<vector<ld> > > bAvg){
	T = 500;
	vector<int> o;
	vector<vector<ld> > alpha(T+1, vector<ld>(N+1, 0.0));
	vector<vector<ld> > beta(T+1, vector<ld>(N+1, 0.0));
	int nDigits = 10;
	int C_N=320;
	int C_p=12;
	int C_window=80;
	int C_cbSize=32;
	bool isHeader = true;

	const char* exeCommand = "C:/Recording_Module/Recording_Module.exe 3 rec.wav rec.txt";
	system(exeCommand);

	vector<vector<ld> > c = calculateCFromScratch("rec.txt", C_N, C_p, C_window, isHeader);
	vector<vector<ld> > codebook = readDim2("codebook.txt", C_p);
	o = makeObservationSeq(c,codebook,C_cbSize);
	T = o.size()-1;
	int num=0;
	ld prob=numeric_limits<ld>::min();
	for(int dig=0;dig<nDigits;dig++){
		solution1(piAvg[dig], aAvg[dig], bAvg[dig], o, alpha, beta);
		if(pBar >= prob){
			prob=pBar;
			num=dig;
		}
	}
	cout<<"Predicted: "<<num<<" "<<prob<<endl;
}

int _tmain(int argc, _TCHAR* argv[])
{
	int nTrain=20;
	int nTest=10;
	int nDigits=10;
	vector<int> o;
	vector<ld> pi;
	vector<vector<ld> > a;
	vector<vector<ld> > b;
	vector<vector<ld> > piBar(nTrain+1, vector<ld>(N+1, 0.0));
	vector<vector<vector<ld> > > aBar(nTrain+1, vector<vector<ld> >(N+1, vector<ld>(N+1, 0.0)));
	vector<vector<vector<ld> > > bBar(nTrain+1, vector<vector<ld> >(N+1, vector<ld>(M+1, 0.0)));
	vector<vector<ld> > piAvg(nDigits+1, vector<ld>(N+1, 0.0));
	vector<vector<vector<ld> > > aAvg(nDigits+1, vector<vector<ld> >(N+1, vector<ld>(N+1, 0.0)));
	vector<vector<vector<ld> > > bAvg(nDigits+1, vector<vector<ld> >(N+1, vector<ld>(M+1, 0.0)));
	vector<int> qStar(T+1, 0);
	vector<vector<ld> > alpha(T+1, vector<ld>(N+1, 0.0));
	vector<vector<ld> > beta(T+1, vector<ld>(N+1, 0.0));


	bool isHeader=true;
	int C_N=320;
	int C_p=12;
	int C_window=80;
	int C_cbSize=32;
	int n_iter=20;
	int n_cycles=3;

	/***************************TRAINING MODULE STARTS************************************************/
	for(int digit=0;digit<nDigits;digit++){
		int cycle=0;
		while(cycle<n_cycles){
			for(int i=1;i<=nTrain;i++){
				string testName = "./Digits_2/160101085_" + to_string((ll)digit) + "_" + to_string((ll)i)+".txt";
				cout<<testName<<endl;
				vector<vector<ld> > c = calculateCFromScratch(testName, C_N, C_p, C_window, isHeader);
				vector<vector<ld> > codebook = readDim2("codebook.txt", C_p);
				o = makeObservationSeq(c,codebook,C_cbSize);
				T = o.size()-1;
				if(cycle==0){
					pi = readDim1("pi_ini.txt");
					a = readDim2("a_ini.txt", N);
					b = readDim2("b_ini.txt", M);
				}
				else{
					pi = piAvg[digit];
					a = aAvg[digit];
					b = bAvg[digit];
				}
				manipulateB(b);

				bool flag = true;
				int itr=1;
				p=0.0;
				pBar=0.0;
				while((pBar > p) || (flag)){
					flag = false;
					p=pBar;
					solution1(pi, a, b, o, alpha, beta);
					solution2(pi, a, b, o, qStar);
					solution3(pi, a, b, o, piBar[i], aBar[i], bBar[i], alpha, beta);
					a = aBar[i];
					b = bBar[i];
					pi = piBar[i];
					manipulateB(b);
					cout<<"Cycle: "<<cycle+1<<" | Digit: "<<digit<<" | Sample: "<<i<<" | Iteration #"<<itr<<": "<<p<<" "<<pBar<<endl;
					if(itr>=n_iter)
						break;
					itr++;
				}
			}
			findAvg(piBar, aBar, bBar, nTrain, piAvg[digit], aAvg[digit], bAvg[digit]);
			cycle++;
		}
		writeOutputToFile(piAvg[digit], aAvg[digit], bAvg[digit], digit);
	}
	/***************************TRAINING MODULE ENDS************************************************/
	//Individual Model parameters (Pi,A,B) are stored in files model_<digit>.txt
	
	/***************************TESTING MODULE STARTS************************************************/
	while(1){
		string mode = chooseMode();
		if(mode == "1")
			runTests(piAvg, aAvg, bAvg);
		else if(mode == "2")
			runRealTime(piAvg, aAvg, bAvg);
		else
			break;
	}
	/***************************TESTING MODULE STARTS************************************************/
	
	getchar();
	return 0;
}