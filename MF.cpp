

#include <vector>
#include <iostream>
#include <ctime>
#include <cmath>
#include <set>
#include <map>
#include <iomanip>
#include <fstream>
#include <string>
using namespace std;

#include "headerfile.h"


class MF{
public:
	MF(){;}
	MF(int _usern, int _itemn, int totalsize):usern(_usern),itemn(_itemn),Totalsize(totalsize){
		K = 10;
		alpha=0.01;
		lambda=0.1;
		TrainingNum=100;
		TrainingLoops=100;
		//select hyper-parameters in {0.1,0.01,0.001,...} to get best performance
		if(dataset_name=="ML100k" && experiment_name=="classical"){
			alpha=0.01;
			lambda=0.1;
		}
		if(dataset_name=="ML100k" && experiment_name=="time_order"){
			alpha=0.01;
			lambda=0.1;
		}
		if(dataset_name=="ML100k" && experiment_name=="time_order_online"){
			alpha=0.01;
			lambda=0.1;
		}
		if(dataset_name=="Epinions" && experiment_name=="classical"){
			alpha=0.001;
			lambda=0.001;
		}
		if(dataset_name=="Epinions" && experiment_name=="time_order"){
			alpha=0.01;
			lambda=0.1;
		}
		if(dataset_name=="Epinions" && experiment_name=="time_order_online"){
			alpha=0.01;
			lambda=0.1;
		}

		//
		user_count=vector<int>(usern,0);
		item_count=vector<int>(itemn,0);
		pu=vector<double*>(usern);
		for(int m=0;m!=usern;m++)
			pu[m]=(double*)malloc(sizeof(double)*K);
		pi=vector<double*>(itemn);
		for(int m=0;m!=itemn;m++)
			pi[m]=(double*)malloc(sizeof(double)*K);
		global_average = 0;
		user_average = vector<double>(usern, 0);
		item_average = vector<double>(itemn, 0);
	}
	void init(double * p){
		for(int k=0;k!=K;k++)
			p[k]=random(-0.01,0.01);
	}
	void update(int userid, int itemid, int rating){
		Rating r;
		r.userid=userid;
		r.itemid=itemid;
		r.rating=rating;
		r.time=-1;
		trainingset.push_back(r);
		if(trainingset.size()%(Totalsize/TrainingNum)==0)
			retrain();
	}
	void push_back(int userid, int itemid, int rating){
		Rating r;
		r.userid=userid;
		r.itemid=itemid;
		r.rating=rating;
		r.time=-1;
		trainingset.push_back(r);
	}
	void retrain(bool show=false){
		global_average = 0;
		user_average = vector<double>(usern, 0);
		item_average = vector<double>(itemn, 0);
		user_count=vector<int>(usern,0);
		item_count=vector<int>(itemn,0);
		for (int m = 0; m != usern; m++) {
			init(pu[m]);
		}
		for (int m = 0; m != itemn; m++) {
			init(pi[m]);
		}
		for(int m=0;m!=trainingset.size();m++){
			int user=trainingset[m].userid;
			int item=trainingset[m].itemid;
			user_count[user]+=1;
			item_count[item]+=1;
			global_average += trainingset[m].rating;
			user_average[user]+=trainingset[m].rating;
			item_average[item]+=trainingset[m].rating;
		}
		if (trainingset.size() != 0)
			global_average /= trainingset.size();
		for (int m = 0; m != usern; m++) {
			if (user_count[m] != 0) {
				user_average[m] /= user_count[m];
			}
		}
		for (int m = 0; m != itemn; m++) {
			if (item_count[m] != 0) {
				item_average[m] /= item_count[m];
			}
		}
		//
		for(int n=0;n!=TrainingLoops;n++){
			if (show) {
				cout << ".";
				if ((n + 1) % 10 == 0) {
					cout << " " << (n + 1) << "%" << endl;
				}
			}
			for(int m=0;m!=trainingset.size();m++){
				int sn=random((int)trainingset.size());
				SGD(trainingset[sn].userid,trainingset[sn].itemid,trainingset[sn].rating);
			}
		}
	}
	double calculate(int userid, int itemid){
		if (user_count[userid] == 0 || item_count[itemid] == 0) {
			if (user_count[userid] == 0 && item_count[itemid] == 0)
				return global_average;
			else if (user_count[userid] == 0)
				return item_average[itemid];
			else 
				return user_average[userid];
		}
		double r=0;
		for(int k=0;k!=K;k++)
			r+=pu[userid][k]*pi[itemid][k];
		return r;
	}
	void SGD(int userid, int itemid, int rating){
		double r=0;
		for(int k=0;k!=K;k++)
			r+=pu[userid][k]*pi[itemid][k];
		double delta=rating-r;
		double du_temp,di_temp;
		for(int k=0;k!=K;k++){
			du_temp=(delta*pi[itemid][k]-lambda*pu[userid][k])*alpha;
			di_temp=(delta*pu[userid][k]-lambda*pi[itemid][k])*alpha;
			pu[userid][k]+=du_temp;
			pi[itemid][k]+=di_temp;
		}
		
	}
	int usern,itemn,K,Totalsize,TrainingNum,TrainingLoops;
	double alpha, lambda;
	vector<Rating> trainingset;
	vector<int> user_count,item_count;
	vector<double*> pu,pi;
	//
	double global_average;
	vector<double> user_average, item_average;
};


void MF_mainfunction(vector<Rating> & ratings, int usern, int itemn){

	MF model1(usern,itemn,ratings.size());
	double rmse=0;
	double rmse1=0,num1=0,num_all=0;
	set<int> appeared_users, appeared_items;
	for(int m=0;m!=ratings.size();++m){
		if((ratings.size()/100)>0 && (m+1)%(ratings.size()/100)==0){
			cout<<".";
		}
		if((ratings.size()/10)>0 && (m+1)%(ratings.size()/10)==0){
			cout<<" "<<(m+1)/(ratings.size()/10)*10<<"%";
			cout<<" Current piece RMSE: "<<sqrt(rmse1/num1);
			cout<<endl;
			rmse1=num1=0;
		}
		double r;
		if(appeared_users.count(ratings[m].userid) && appeared_items.count(ratings[m].itemid)){
			r=model1.calculate(ratings[m].userid,ratings[m].itemid);
			rmse+=pow(ratings[m].rating-r,2);
			rmse1+=pow(ratings[m].rating-r,2);
			num1+=1;
			num_all+=1;
		}
		appeared_users.insert(ratings[m].userid);
		appeared_items.insert(ratings[m].itemid);
		model1.update(ratings[m].userid,ratings[m].itemid,ratings[m].rating);
	}
	cout << "Overall RMSE: " << endl;
	cout<<sqrt(rmse/num_all)<<endl;


}


void MF_mainfunction_batch(vector<Rating> & rtraining, vector<Rating> & rtest, int usern, int itemn){

	MF model1(usern,itemn,rtraining.size()+rtest.size());
	set<int> appeared_users, appeared_items;
	cout << "Training" << endl;
	for(int m=0;m!=rtraining.size();++m){
		appeared_users.insert(rtraining[m].userid);
		appeared_items.insert(rtraining[m].itemid);
		model1.push_back(rtraining[m].userid,rtraining[m].itemid,rtraining[m].rating);
	}
	model1.retrain(true);
	cout << "Test" << endl;
	double rmse=0,num_all=0;
	for(int m=0;m!=rtest.size();++m){
		if(!appeared_users.count(rtest[m].userid) || !appeared_items.count(rtest[m].itemid))
			continue;
		double r=model1.calculate(rtest[m].userid,rtest[m].itemid);
		rmse+=pow(rtest[m].rating-r,2);
		num_all+=1;
	}
	
	cout<<num_all<<" ";
	cout << "RMSE: " ;
	cout<<sqrt(rmse/num_all)<<endl;

}



