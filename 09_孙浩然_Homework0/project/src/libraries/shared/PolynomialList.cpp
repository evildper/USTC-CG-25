#include "PolynomialList.h"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>

#include <assert.h>

#define EPSILON 1.0e-10	// zero double

using namespace std;



PolynomialList::PolynomialList(const PolynomialList& other) {
    // TODO
    m_Polynomial = other.m_Polynomial;
}

PolynomialList::PolynomialList(const string& file) {
    // TODO
    ReadFromFile(file);
}

PolynomialList::PolynomialList(const double* cof, const int* deg, int n) {
    // 
    for(int i = 0;i<n;i++)
        AddOneTerm(Term(deg[i],cof[i]));
}

PolynomialList::PolynomialList(const vector<int>& deg, const vector<double>& cof) {
    // TODO
    assert(deg.size()==cof.size());
    for(int i = 0;i<cof.size();i++)
        AddOneTerm(Term(deg[i],cof[i]));
}

double PolynomialList::coff(int i) const {
    // TODO
    //bravo
    for(const Term& term : m_Polynomial){
        if(term.deg<i){
            break;
        }
        if(term.deg ==i){
            return term.cof;
        }
    }
    return 0.; // you should return a correct value
}

double& PolynomialList::coff(int i) {
    // TODO
    static double ERROR; // you should delete this line
        for(Term& term : m_Polynomial){
        if(term.deg<i){
            break;
        }
        if(term.deg ==i){
            return term.cof;
        }
    }
    return ERROR; // you should return a correct value
}
//clear the very little item to make the polynomial more looking 
void PolynomialList::compress() {
    // TODO
    auto lr = m_Polynomial.begin();//lr is a address
    while(lr != m_Polynomial.end()){
        if(fabs((*lr).cof) < EPSILON)//more ()
            lr = m_Polynomial.erase(lr);//the method of the std::vector class
        
        else
            lr++;
        
    }
}



PolynomialList PolynomialList::operator+(const PolynomialList& right) const {
//     // TODO
    PolynomialList pol(*(this));
    for(const Term& term:right.m_Polynomial){

    
        pol.AddOneTerm(term);
    }    
    pol.compress();
    return pol; // you should return a correct value
}

PolynomialList PolynomialList::operator-(const PolynomialList& right) const {
    PolynomialList pol(*(this));//this point at oneself's address
    for(const Term& term:right.m_Polynomial)
        pol.AddOneTerm(Term(term.deg,-term.cof));
    pol.compress();
    return pol; // you should return a correct value
}

PolynomialList PolynomialList::operator*(const PolynomialList& right) const {
    // TODO
    PolynomialList pol;
    int deg = 0;
    double cof = 0.;
    for(const Term& term_1:m_Polynomial){
        for(const auto& term_2:right.m_Polynomial){//auto is more flexible
            deg = term_1.deg + term_2.deg;
            cof = term_1.cof*term_2.cof;
            pol.AddOneTerm(Term(deg,cof));
        }
    }
    return pol; // you should return a correct value
}

PolynomialList& PolynomialList::operator=(const PolynomialList& right) {
    // TODO
    m_Polynomial = right.m_Polynomial;
    return *this;
}

void PolynomialList::Print() const {
    auto itr = m_Polynomial.begin();
    if (itr == m_Polynomial.end()) {
        cout << "0" << endl;
        return;
    }

    for (; itr != m_Polynomial.end(); itr++) {
        if (itr != m_Polynomial.begin()) {
            cout << " ";
            if (itr->cof > 0)
                cout << "+";
        }

        cout << itr->cof;

        if (itr->deg > 0)
            cout << "x^" << itr->deg;
    }
    cout << endl;
}

bool PolynomialList::ReadFromFile(const string& file) {
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
        Term nd;
        inp >> nd.deg;
        inp >> nd.cof;

        AddOneTerm(nd);
    }

    inp.close();

    return true;
}

PolynomialList::Term& PolynomialList::AddOneTerm(const Term& term) {
    auto itr = m_Polynomial.begin();
    for (; itr != m_Polynomial.end(); itr++) {
        if (itr->deg == term.deg) {
            itr->cof += term.cof;
            return *itr;
        }

        if (itr->deg > term.deg)
            break;
    }
    return *m_Polynomial.insert(itr, term);
}
