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

void loading(const string & filename, vector<Rating> & ratings, map<int, int> & user_map, map<int, int> & item_map, int & T0) {
	user_map.clear();
	item_map.clear();
	ratings.clear();
	fstream f1;
	f1.open(filename.c_str(), ios::in);
	while (!f1.eof()) {
		string line;
		getline(f1, line);
		if (line.empty())
			continue;
		stringstream sstr(line);
		int userid, itemid, rating, time;
		sstr >> userid;
		if (sstr.fail())continue;
		sstr >> itemid;
		if (sstr.fail())continue;
		sstr >> rating;
		if (sstr.fail())continue;
		sstr >> time;
		if (sstr.fail())continue;
		if (!user_map.count(userid)) {
			int num = user_map.size();
			user_map[userid] = num;
		}
		if (!item_map.count(itemid)) {
			int num = item_map.size();
			item_map[itemid] = num;
		}
		Rating r;
		r.userid = user_map[userid];
		r.itemid = item_map[itemid];
		r.rating = rating;
		r.time = time;
		ratings.push_back(r);
	}
	if(!ratings.empty()){
		T0 = ratings[0].time;
		for (int n = 0; n != ratings.size(); n++) {
			if (ratings[n].time < T0)
				T0 = ratings[n].time;
		}
		for (int n = 0; n != ratings.size(); n++)
			ratings[n].time -= T0;
		map<int, int> user_map_swap, item_map_swap;
		for (map<int, int>::iterator it = user_map.begin(); it != user_map.end(); ++it)
			user_map_swap[it->second] = it->first;
		for (map<int, int>::iterator it = item_map.begin(); it != item_map.end(); ++it)
			item_map_swap[it->second] = it->first;
		user_map = user_map_swap;
		item_map = item_map_swap;
	}
	cout	<< "Load a dataset with " 
			<< user_map.size() << " users, " 
			<< item_map.size() << " items, " 
			<< ratings.size() << " ratings." << endl;
}


bool rating_earlier_than(const Rating & r1, const Rating & r2){
	return r1.time<r2.time;
}

void sort_ratings_by_time(vector<Rating> & ratings) {
	sort(ratings.begin(), ratings.end(), rating_earlier_than);
}

void sort_ratings_random(vector<Rating> & ratings) {
	srand(1);
	vector<pair<int, int> > ratings_num(ratings.size());
	for (int n = 0; n != ratings_num.size(); n++) {
		ratings_num[n].first = rand();
		ratings_num[n].second = n;
	}
	sort(ratings_num.begin(), ratings_num.end());
	vector<Rating> ratings_after_sort(ratings.size());
	for (int n = 0; n != ratings_num.size(); n++) {
		ratings_after_sort[n] = ratings[ratings_num[n].second];
	}
	ratings = ratings_after_sort;
}

void divide_80_20(const vector<Rating> & ratings, vector<Rating> & rtraining, vector<Rating> & rtest) {
	rtraining.clear();
	rtest.clear();
	int trainingN = ((int)(ratings.size()*0.8));
	for (int n = 0; n != ratings.size(); n++) {
		if (n < trainingN)
			rtraining.push_back(ratings[n]);
		else
			rtest.push_back(ratings[n]);
	}
}

void timesrand() {
	srand((int)time(0));
	rand(); rand();
}

double random(double a, double b){
	return a+(b-a)*(rand()*(RAND_MAX+1.0)+rand())/(RAND_MAX+1.0)/(RAND_MAX+1.0);
}

int random(int n) {
	return (rand()*(RAND_MAX + 1) + rand()) % n;
}