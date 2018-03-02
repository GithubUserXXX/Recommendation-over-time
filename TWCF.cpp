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

extern string dataset_name;
extern string experiment_name;

class TWCF{
public:
	TWCF(){;}
	TWCF(int _usern, int _itemn):usern(_usern),itemn(_itemn){
		T=1.0;
		//select hyper-parameters in {1.0, 10.0, ...} to get best performance
		if(dataset_name=="ML100k" && experiment_name=="classical"){
			T=1.0;
		}
		if(dataset_name=="ML100k" && experiment_name=="time_order"){
			T=1.0;
		}
		if(dataset_name=="ML100k" && experiment_name=="time_order_online"){
			T=1.0;
		}
		if(dataset_name=="Epinions" && experiment_name=="classical"){
			T=10.0;
		}
		if(dataset_name=="Epinions" && experiment_name=="time_order"){
			T = 10.0;
		}
		if(dataset_name=="Epinions" && experiment_name=="time_order_online"){
			T = 10.0;
		}
		Avg=0;
		T0 = 86400;
		user_item=vector< map<int,double> >(usern);
		item_user=vector< map<int,double> >(itemn);
		user_item_time=vector< map<int,double> >(usern);
		user_avg=vector<double>(usern,Avg);
	}
	double f(double t){
		if(t<0)
			t=-t;
		return exp(-t/T);
	}
	void update(int userid, int itemid, double rating, double time){
		time /= T0;
		double rating_old;
		if(user_item[userid].count(itemid)){
			rating_old=user_item[userid][itemid];
		}else{
			rating_old=Avg;
		}
		user_avg[userid]+=(rating-rating_old)/itemn;
		user_item[userid][itemid]=rating;
		item_user[itemid][userid]=rating;
		user_item_time[userid][itemid]=time;
	}
	double calculate(int userid, int itemid, double time){
		time /= T0;
		double sim_p21=0;
		for(int user1=0;user1!=usern;user1++){
			if(item_user[itemid].count(user1)){
				sim_p21+=(item_user[itemid][user1]-user_avg[user1])*(item_user[itemid][user1]-user_avg[user1]);
			}else{
				sim_p21+=(Avg-user_avg[user1])*(Avg-user_avg[user1]);
			}
		}
		double p1=0,p2=0;
		for(map<int,double>::iterator it1=user_item[userid].begin();it1!=user_item[userid].end();++it1){
			int item1=it1->first;
			double rating1=it1->second;
			double similar=0;
			if(item1==itemid)
				similar=1;
			else if(sim_p21==0)
				similar=0;
			else {
				double sim_p1=0;
				double sim_p22=0;
				for(int user1=0;user1!=usern;user1++){
					if(item_user[item1].count(user1)){
						if(item_user[itemid].count(user1)){
							sim_p1+=(item_user[itemid][user1]-user_avg[user1])*(item_user[item1][user1]-user_avg[user1]);
						}else{
							sim_p1+=(Avg-user_avg[user1])*(item_user[item1][user1]-user_avg[user1]);
						}
						sim_p22+=(item_user[item1][user1]-user_avg[user1])*(item_user[item1][user1]-user_avg[user1]);
					}else{
						if(item_user[itemid].count(user1)){
							sim_p1+=(item_user[itemid][user1]-user_avg[user1])*(Avg-user_avg[user1]);
						}else{
							sim_p1+=(Avg-user_avg[user1])*(Avg-user_avg[user1]);
						}
						sim_p22+=(Avg-user_avg[user1])*(Avg-user_avg[user1]);
					}
				}
				if(sim_p22==0)
					similar=0;
				else similar=sim_p1/sqrt(sim_p21*sim_p22);
				if(similar<0)
					similar=0;
			}
			double func=f(time-user_item_time[userid][itemid]);
			p1+=rating1*similar*func;
			p2+=similar*func;
		}
		if(p2==0)
			return 0;
		return p1/p2;
	}
	int usern,itemn;
	double Avg,T,T0;
	vector< map<int,double> > user_item;
	vector< map<int,double> > item_user;
	vector< map<int,double> > user_item_time;
	vector<double> user_avg;
	
};



void TWCF_mainfunction(vector<Rating> & ratings, int usern, int itemn){
	TWCF model1(usern,itemn);
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
			r=model1.calculate(ratings[m].userid,ratings[m].itemid,ratings[m].time);
			rmse+=pow(ratings[m].rating-r,2);
			rmse1+=pow(ratings[m].rating-r,2);
			num1+=1;
			num_all+=1;
		}
		appeared_users.insert(ratings[m].userid);
		appeared_items.insert(ratings[m].itemid);
		model1.update(ratings[m].userid,ratings[m].itemid,ratings[m].rating,ratings[m].time);
	}
	cout << "Overall RMSE: " << endl;
	cout<<sqrt(rmse/num_all)<<endl;


}


void TWCF_mainfunction_batch(vector<Rating> & rtraining, vector<Rating> & rtest, int usern, int itemn){

	TWCF model1(usern,itemn);
	set<int> appeared_users, appeared_items;
	cout << "Training" << endl;
	for(int m=0;m!=rtraining.size();++m){
		appeared_users.insert(rtraining[m].userid);
		appeared_items.insert(rtraining[m].itemid);
		model1.update(rtraining[m].userid,rtraining[m].itemid,rtraining[m].rating,rtraining[m].time);
	}
	cout << "Test" << endl;
	double rmse=0,num_all=0;
	for(int m=0;m!=rtest.size();++m){
		if((rtest.size()/100)>0 && (m+1)%(rtest.size()/100)==0){
			cout<<".";
		}
		if((rtest.size()/10)>0 && (m+1)%(rtest.size()/10)==0){
			cout<<" "<<(m+1)/(rtest.size()/10)*10<<"%"<<endl;
		}
		if(!appeared_users.count(rtest[m].userid) || !appeared_items.count(rtest[m].itemid))
			continue;
		double r=model1.calculate(rtest[m].userid,rtest[m].itemid,rtest[m].time);
		rmse+=pow(rtest[m].rating-r,2);
		num_all+=1;
	}
	
	cout<<num_all<<" ";
	cout << "RMSE: " ;
	cout<<sqrt(rmse/num_all)<<endl;

}

