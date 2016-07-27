#pragma warning(disable:4996)
#include <iostream>
#include <set>
#include <vector>
#include <fstream>
#include <string>
using namespace std;
class Berlekamp_Massey_Struct {
public:
	unsigned int linear_span = 0;
	set<unsigned int> polynomial;
	bool primitive_polynomial = false;
	unsigned int constant = 0;
	void clear();
};
class Berlekamp_Massey : public Berlekamp_Massey_Struct {
public:
	Berlekamp_Massey_Struct Berlekamp_Massey_Algorithm(unsigned int sequence[], unsigned int N);
	string print_poly(set<unsigned int> polynomial);
	set<unsigned int> Berlekamp_Massey::sysmetric_different(set<unsigned int> ob1, set<unsigned int> ob2);
};
Berlekamp_Massey_Struct Berlekamp_Massey::Berlekamp_Massey_Algorithm(unsigned int sequence[], unsigned int N) {
	Berlekamp_Massey_Struct result;
	set<unsigned int> f, g;
	unsigned int k = 0;
	for (k = 0; k < N; ++k) {
		if (*(sequence + k) == 1) {
			break;
		}
	}
	// use a set to denote polynomial
	f.insert(0);
	f.insert(k + 1);
	unsigned int l = k + 1;

	g.insert(0);
	unsigned int a = k;
	unsigned int b = 0;

	for (unsigned int n = k + 1, d = 0; n < N; ++n) {
		d = 0;
		for (unsigned int ele : f) {
			d ^= sequence[ele + n - l];
		}
		if (d == 0) {
			b += 1;
		}
		else {
			if ((2 * l) > n) {
				set<unsigned int> mt;
				for (unsigned int t : g) {
					mt.insert(a - b + t);
				}
				f = sysmetric_different(f, mt);
				b += 1;
			}
			else {
				set<unsigned int> temp = f;
				set<unsigned int> mt;
				for (unsigned int t : f) {
					mt.insert(b - a + t);
				}
				f.clear();
				for (unsigned int t : mt) {
					f.insert(t);
				}
				for (unsigned int t : g) {
					f.insert(t);
				}
				l = n + 1 - l;
				g = temp;
				a = b;
				b = n - l + 1;
			}
		}
	}
	// save result
	result.linear_span = l;
	result.polynomial = f;
	result.constant = f.size();
	result.primitive_polynomial = ((result.constant & 1) ? true : false)&f.find(0) != f.end();
	return result;
}
string Berlekamp_Massey::print_poly(set<unsigned int> polynomial) {
	string result = "";
	unsigned int i = 0;
	for (set<unsigned int>::reverse_iterator rit = polynomial.rbegin(); rit != polynomial.rend(); ++rit) {
		if (*rit == 0) {
			result += '1';
		}
		else {
			result.append("x^");
			result.append(to_string(*rit));
			result.append("+");
		}
	}
	return result;
}
// 대칭차집합
set<unsigned int> Berlekamp_Massey::sysmetric_different(set<unsigned int> ob1, set<unsigned int> ob2)
{
	set<unsigned int> newset = ob1;
	for (unsigned int i : ob2) {
		newset.insert(i);
	}
	set<unsigned int> newset2;
	for (unsigned int i : ob1)
	{
		if (ob2.find(i) != ob2.end()) {
			newset2.insert(i);
		}
	}
	for (unsigned int i : newset2) {
		newset.erase(i);
	}
	return newset;
}
void Berlekamp_Massey_Struct::clear() {
	linear_span = 0;
	primitive_polynomial = false;
	constant = 0;
}
void main() {
	unsigned int sequence[] = { 1,0,0,0,0,0,0,1,1,0,0,1,1 };	// initalize array value
	const unsigned int N = sizeof(sequence) / sizeof(sequence[0]);
	Berlekamp_Massey_Struct mBMs;
	Berlekamp_Massey BM;
	
	mBMs = BM.Berlekamp_Massey_Algorithm(sequence, N);
	cout << "다항식 길이 : " << mBMs.linear_span << "\t 다항식 차수 갯수 : " << mBMs.constant << endl;
	cout << BM.print_poly(mBMs.polynomial) << endl;
	cout << "=============== Search Finish!!! ===============" << endl;
}