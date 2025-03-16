#include "PolynomialMap.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#define EPSION 1.0e-10
using namespace std;

PolynomialMap::PolynomialMap(const PolynomialMap& other) {
    // TODO
	m_Polynomial = other.m_Polynomial;

}

PolynomialMap::PolynomialMap(const string& file) {
    ReadFromFile(file);
}

PolynomialMap::PolynomialMap(const double* cof, const int* deg, int n) {
	// TODO
	for(int i=0;i<n;i++){
		coff(deg[i]) = cof[i];//double & : this can be changed
	}
}

PolynomialMap::PolynomialMap(const vector<int>& deg, const vector<double>& cof) {
	assert(deg.size() == cof.size());
	// TODO
	for (int i = 0; i < cof.size(); i++)
	{
		coff(deg[i])=cof[i];
	}
	
}

double PolynomialMap::coff(int i) const {
	// TODO
	auto n = m_Polynomial.find(i);
	if(n == m_Polynomial.end()){
		return 0.;
	}
	else{
		return n->second;
	}
	return 0.f; // you should return a correct value
}

double& PolynomialMap::coff(int i) {
	// TODO
	//static double ERROR; // you should delete this line
	return m_Polynomial[i]; // you should return a correct value
}

void PolynomialMap::compress() {
	// TODO
	map<int,double> temp = m_Polynomial;
	m_Polynomial.clear();
	for(const auto& t:temp){
		if (fabs(t.second)>EPSION)
		{
			coff(t.first)=t.second;
			/* code */
		}
		
	}
}

PolynomialMap PolynomialMap::operator+(const PolynomialMap& right) const {
	// TODO
	PolynomialMap pol(right);
	for(const auto& temp:m_Polynomial){
		pol.coff(temp.first)+=temp.second;
	}
	pol.compress();
	return pol; // you should return a correct value
}

PolynomialMap PolynomialMap::operator-(const PolynomialMap& right) const {
	// TODO
	PolynomialMap pol(*(this));
	for(const auto& temp:right.m_Polynomial){
		pol.coff(temp.first)+=(-temp.second);
	}
	pol.compress();	
	return pol; // you should return a correct value
}

PolynomialMap PolynomialMap::operator*(const PolynomialMap& right) const {
	// TODO
	PolynomialMap pol;
	for(const auto& temp1:m_Polynomial){
		for(const auto& temp2:right.m_Polynomial){
			int deg = temp1.first+temp2.first;
			double cof = temp1.second*temp2.second;
			pol.coff(deg) += cof;//remember get the same deg so it is += not =
		}
	}
	return pol; // you should return a correct value
}

PolynomialMap& PolynomialMap::operator=(const PolynomialMap& right) {
	// TODO
	m_Polynomial = right.m_Polynomial;
	return *this;
}

void PolynomialMap::Print() const {
    auto itr = m_Polynomial.begin();
    if (itr == m_Polynomial.end()) {
        cout << "0" << endl;
        return;
    }

    for (; itr != m_Polynomial.end(); itr++) {
        if (itr != m_Polynomial.begin()) {
            cout << " ";
            if (itr->second > 0)
                cout << "+";
        }

        cout << itr->second;

        if (itr->first > 0)
            cout << "x^" << itr->first;
    }
    cout << endl;
}

bool PolynomialMap::ReadFromFile(const string& file)
{
    m_Polynomial.clear();

    ifstream inp;
    inp.open(file.c_str());
    if (!inp.is_open()) {
        cout << "ERROR::PolynomialList::ReadFromFile:" << endl
            << "\t" << "file [" << file << "] opens failed" << endl;
        return false;
    }

    char ch;
    int n;
    inp >> ch;
    inp >> n;
    for (int i = 0; i < n; i++) {
        int deg;
        double cof;
        inp >> deg;
        inp >> cof;
        coff(deg) = cof;
    }

    inp.close();

    return true;
}
