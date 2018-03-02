#ifndef HEADERFILE_H
#define HEADERFILE_H

#include <vector>
#include <map>
using namespace std;


class Rating{
public:
	int userid;
	int itemid;
	int rating;
	int time;
};


extern string dataset_name;
extern string experiment_name;
//load_data.cpp
void loading(const string & filename, vector<Rating> & ratings, map<int, int> & user_map, map<int, int> & item_map, int & T0);
void sort_ratings_by_time(vector<Rating> & ratings);
void sort_ratings_random(vector<Rating> & ratings);
void divide_80_20(const vector<Rating> & ratings, vector<Rating> & rtraining, vector<Rating> & rtest);
void timesrand();
double random(double a, double b);
int random(int n);
//expm_logm.cpp
double * matrix_expm(double * X, int K);
double * matrix_logm(double * X, int K);
//CTHMM.cpp
void CTHMM_mainfunction(vector<Rating> & ratings, int usern, int itemn);
void CTHMM_mainfunction_batch(vector<Rating> & rtraining, vector<Rating> & rtest, int usern, int itemn, bool classical);
//MF.cpp
void MF_mainfunction(vector<Rating> & ratings, int usern, int itemn);
void MF_mainfunction_batch(vector<Rating> & rtraining, vector<Rating> & rtest, int usern, int itemn);
//CKF.cpp
void CKF_mainfunction(vector<Rating> & ratings, int usern, int itemn);
void CKF_mainfunction_batch(vector<Rating> & rtraining, vector<Rating> & rtest, int usern, int itemn);
//TWCF.cpp
void TWCF_mainfunction(vector<Rating> & ratings, int usern, int itemn);
void TWCF_mainfunction_batch(vector<Rating> & rtraining, vector<Rating> & rtest, int usern, int itemn);







#endif
