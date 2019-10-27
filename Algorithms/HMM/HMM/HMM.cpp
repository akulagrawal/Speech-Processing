// HMM.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<iomanip>
#include<algorithm>
#include<math.h>

#define ld long double
#define pb push_back
#define THRESHOLD pow((ld)10,-30)

using namespace std;

int N=5, M=32, T=85;
ld p, pBar, pStar;
vector<int> o;
vector<ld> pi;
vector<vector<ld> > a;
vector<vector<ld> > b;
vector<ld> piBar(N+1, 0.0);
vector<vector<ld> > aBar(N+1, vector<ld>(N+1, 0.0));
vector<vector<ld> > bBar(N+1, vector<ld>(M+1, 0.0));
vector<int> qStar(T+1, 0);
vector<vector<ld> > alpha(T+1, vector<ld>(N+1, 0.0));
vector<vector<ld> > beta(T+1, vector<ld>(N+1, 0.0));
vector<vector<ld> > delta(T+1, vector<ld>(N+1, 0.0));
vector<vector<int> > psi(T+1, vector<int>(N+1, 0));
vector<vector<vector<ld> > > xi(T+1, vector<vector<ld> >(N+1, vector<ld>(N+1, 0.0)));
vector<vector<ld> > gamma(T+1, vector<ld>(N+1, 0.0));

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

void print(){
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

void solution1(){
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

void solution2(){
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

void solution3(){
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

void manipulateB(){
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

int _tmain(int argc, _TCHAR* argv[])
{
	pi = readDim1("pi_ini.txt");
	a = readDim2("a_ini.txt", N);
	b = readDim2("b_ini.txt", M);
	o = readO("HMM_OBSERVATION_SEQUENCE_1.txt");
	manipulateB();
	print();
	bool flag = true;

	//while((pBar > p) || (flag)){    //uncomment this, and comment next line for full training
	for(int i=0;i<20;i++){
		flag = false;

		p=pBar;

		// SOLUTION 1
		solution1();

		// SOLUTION 2
		solution2();

		// SOLUTION 3
		solution3();
		

		a = aBar;
		b = bBar;
		pi = piBar;
		manipulateB();
		cout<<p<<" "<<pBar<<endl;
	}
	
	cout<<"\nQ:\n";
	for(int i=1;i<=T;i++){
		cout<<qStar[i]<<" ";
	}
	cout<<endl;
	print();


	getchar();
	return 0;
}