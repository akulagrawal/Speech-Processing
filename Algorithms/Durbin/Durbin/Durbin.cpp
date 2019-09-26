// Durbin.cpp : Defines the entry point for the console application.
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

/* Read CoolEdit data from file "name" and store it in a vector which is
   returned by the function. Heading comprising of top 5 lines is stored
   in vector "heading" if header = true
   if nSamples < 0, whole file is stored, else only first nSamples lines are read and stored
*/
vector<ld> readFileToVec(string name, vector<string> &heading, bool header = true, int nSamples = -1){
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

int _tmain(int argc, _TCHAR* argv[])
{
	//These values need to be set correctly before running the program
	bool applyWindow = false;
	bool isHeader = false;
	int N = 320;
	int p = 12;
	string dir = "./";
	string fileName = "test.txt";

	vector<ld> v;
	vector<ld> r;
	vector<ld> a;
	vector<ld> c;
	vector<string> heading;
	
	v = readFileToVec(dir + fileName, heading, isHeader, N);
	if(applyWindow)
		applyHammingWindow(v, N);
	r = calculateR(v, N, p);
	a = durbin(r, p);
	c = calculateC(a, r, p);
	print_output(r, a, c, p);
	getchar();
	return 0;
}

