// LBGAlgo.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stdio.h"
#include "iostream"
#include "fstream"
#include "vector"
#include "cmath"
#include "string"
#include "iomanip"
#include "sstream"

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

vector<vector<ld> > readFileToVec() {
	vector<vector<ld> > x;
	vector<ld> temp;
	// File pointer 
    ifstream fin; 
  
    // Open an existing file 
    fin.open("Universe.csv");
	string unused;
    while ( getline(fin, unused) ){
		temp.clear();
		stringstream iss(unused);
		while(getline (iss, unused, ','))
			temp.pb(stold(unused));
		x.pb(temp);
    }
	return x;
}

ld calculateTokhuraDistance(vector<ld> v1, vector<ld> v2){
	ld wt[12] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
	ld sum = 0;
	for(int i=0;i<12;i++)
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

int _tmain(int argc, _TCHAR* argv[])
{
	int p = 12;
	int K = 8;
	int epsilon = 0.03;
	vector<vector<ld> > x = readFileToVec();
	vector<vector<ld> > centroids(K, vector<ld>(p, 0.0));
	centroids = LBG(x, K, p, epsilon, 10.0);
	for(int i=0;i<K;i++) {
		cout<<i+1<<": ";
		for(int j=0;j<p;j++) {
			cout<<centroids[i][j]<<" ";
		}
		cout<<endl;
	}
	getchar();
	return 0;
}

