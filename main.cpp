#include <vector>
#include <iostream>
#include <ctime>
#include <cmath>
#include <set>
#include <map>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
using namespace std;

#include "headerfile.h"


string dataset_name;
string experiment_name;
string algorithm_name;

int main() {
	{
		int experiment_num, dataset_num, algorithm_num;
		cout << "Please choose experiment type: 1.classical 2.time_order 3.time_order_online" << endl<<"> ";
		cin >> experiment_num;
		cout << "Please choose dataset: 1.MovieLens100k 2.Epinions" << endl<<"> ";
		cin >> dataset_num;
		cout << "Please choose algorithm: 1.CTHMM 2.MF 3.CKF 4.TWCF" << endl<<"> ";
		cin >> algorithm_num;
		experiment_name=(experiment_num == 1?"classical":(experiment_num == 2?"time_order":"time_order_online"));
		dataset_name = (dataset_num == 1 ? "ML100k" : "Epinions");
		algorithm_name = (algorithm_num == 1 ? "CTHMM" :(algorithm_num == 2?"MF":(algorithm_num == 3?"CKF":"TWCF")) );
		cout << "Experitment:" << experiment_name << " " << "Dataset:" << dataset_name << " " << "Algorithm:" << algorithm_name << endl;
	}
	int T0;
	map<int, int> user_map, item_map;
	vector<Rating> ratings,rtraining,rtest;
	{
		if (dataset_name == "ML100k") {
			{
				fstream f1;
				f1.open("Datasets\\MovieLens100k\\u.data", ios::in);
				if (f1.fail()) {
					cout << "Please download MovieLens100k dataset. See Datasets\\README.txt" << endl;
					system("pause");
					return 0;
				}
				f1.close();
			}
			loading("Datasets\\MovieLens100k\\u.data", ratings, user_map, item_map, T0);
		}
		if (dataset_name == "Epinions")
			loading("Datasets\\Epinions\\epinions.txt",ratings, user_map, item_map, T0);
	}
	if (experiment_name=="classical") {
		cout << "Classical experiment."<<endl;
		sort_ratings_random(ratings);
		divide_80_20(ratings, rtraining, rtest);
	}
	else if (experiment_name=="time_order") {
		cout << "Time order experiment."<<endl;
		sort_ratings_by_time(ratings);
		divide_80_20(ratings, rtraining, rtest);
	}
	else if(experiment_name=="time_order_online"){
		cout << "Time order online experiment."<<endl;
		sort_ratings_by_time(ratings);
	}
	timesrand();
	if (algorithm_name == "CTHMM") {
		if (experiment_name == "classical" || experiment_name == "time_order")
			CTHMM_mainfunction_batch(rtraining, rtest, user_map.size(), item_map.size(), true);
		if (experiment_name == "time_order_online")
			CTHMM_mainfunction(ratings, user_map.size(), item_map.size());
	}else if(algorithm_name == "MF") {
		if (experiment_name == "classical" || experiment_name == "time_order")
			MF_mainfunction_batch(rtraining, rtest, user_map.size(), item_map.size());
		if (experiment_name == "time_order_online")
			MF_mainfunction(ratings, user_map.size(), item_map.size());
	}
	else if(algorithm_name == "CKF") {
		if (experiment_name == "classical" || experiment_name == "time_order")
			CKF_mainfunction_batch(rtraining, rtest, user_map.size(), item_map.size());
		if (experiment_name == "time_order_online")
			CKF_mainfunction(ratings, user_map.size(), item_map.size());
	}
	else if(algorithm_name == "TWCF") {
		if (experiment_name == "classical" || experiment_name == "time_order")
			TWCF_mainfunction_batch(rtraining, rtest, user_map.size(), item_map.size());
		if (experiment_name == "time_order_online")
			TWCF_mainfunction(ratings, user_map.size(), item_map.size());
	}
	system("pause");
	return 0;
}
