

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


double random(double a, double b);



double * zeroMatrix(int M, int N) {
	double * Y=(double *)malloc(M*N*sizeof(double));
	memset(Y, 0, M*N * sizeof(double));
	return Y;
}



class CTHMM {
public:
	CTHMM() { ; }
	CTHMM(int _usern, int _itemn, bool classical)
		:usern(_usern), itemn(_itemn), is_classical_experiment(classical) {
		//The hyper-parameters
		K = 10;				//The number of user types (and item types, we set J=K for programming convenience)
		N = 5;				//the number of levels of ratings ( ratings are in {1,2,...,N} )
		M = 100;			//max number of recent ratings
		lambda_1 = 1;		//regularization for pi and omega
		lambda_2 = 0.01;	//regularization for A and B
		T0 = 86400;			

		pu = vector< map<double, double *> >(usern);
		pi = vector< map<double, double *> >(itemn);
		ratings_user = vector< map<double, vector< pair<int, double> > > >(usern);
		ratings_item = vector< map<double, vector< pair<int, double> > > >(itemn);
		ratingnum_user = vector<int>(usern, 0);
		ratingnum_item = vector<int>(itemn, 0);
		B0_user = vector<double*>(usern);
		B0_item = vector<double*>(itemn);

		Xi_user = vector< double* >(usern);
		Xi_item = vector< double* >(itemn);
		for (int m = 0; m != usern; ++m) {
			Xi_user[m] = zeroMatrix(K, K);
			B0_user[m] = (double*)malloc(sizeof(double)*K);
			for (int k = 0; k != K; k++)
				B0_user[m][k] = 1.0;
		}
		for (int m = 0; m != itemn; ++m) {
			Xi_item[m] = zeroMatrix(K, K);
			B0_item[m] = (double*)malloc(sizeof(double)*K);
			for (int k = 0; k != K; k++)
				B0_item[m][k] = 1.0;
		}
		A_user = zeroMatrix(K, K);
		A_item = zeroMatrix(K, K);
		Pxy = zeroMatrix(K, K);
		Pxy1 = zeroMatrix(K, K);
		Pxy2 = zeroMatrix(K, K);
		Xi_users = zeroMatrix(K, K);
		Xi_items = zeroMatrix(K, K);
		P_users = zeroMatrix(K, 1);
		P_items = zeroMatrix(K, 1);
		Pi_users = zeroMatrix(K, 1);
		Pi_items = zeroMatrix(K, 1);

		for (int k = 0; k != K; ++k) {
			P_users[k] = lambda_1;
			P_items[k] = lambda_1;
		}
		P_i(P_users, Pi_users);
		P_i(P_items, Pi_items);

		for (int k1 = 0; k1 != K; ++k1) {
			for (int k2 = 0; k2 != K; ++k2) {
				Pxy1[k1*K + k2] = random(0.0, 1.0);
				Pxy2[k1*K + k2] = 1.0;
			}
		}
		P_xy();

		initialize_Xi(Xi_users);
		initialize_Xi(Xi_items);
		Xi_N_user = 1;
		Xi_N_item = 1;
		Matrix_A(true);
		Matrix_A(false);

		initialize_vector = (double*)malloc(sizeof(double)*(usern + itemn)*K);
		for (int n = 0; n != (usern + itemn); n++)
			initialize(initialize_vector + n * K);

		;
	}
	void initialize(double * p) {
		double r = 0;
		while (r == 0) {
			r = 0;
			for (int k = 0; k != K; ++k) {
				p[k] = random(0.0, 1.0);
				r += p[k];
			}
		}
		for (int k = 0; k != K; ++k)
			p[k] /= r;
	}
	void initialize(double * p, bool isuser, int userid) {
		memcpy(p, (isuser ? initialize_vector + userid * K : initialize_vector + usern + userid * K), sizeof(double)*K);
	}
	void initialize_Xi(double * xi_users) {
		for (int k1 = 0; k1 != K; ++k1) {
			double r = 0;
			for (int k2 = 0; k2 != K; ++k2) {
				xi_users[k1*K + k2] = random(0.0, 1.0);
				r += xi_users[k1*K + k2];
			}
			if (r == 0) {
				k1 -= 1;
				continue;
			}
			for (int k2 = 0; k2 != K; ++k2)
				xi_users[k1*K + k2] = lambda_2 * (xi_users[k1*K + k2] / r) + (1 - lambda_2)*(k1 == k2 ? 1.0 : 0.0);
		}
		double * dQ = matrix_logm(xi_users, K);
		memcpy(xi_users, dQ, sizeof(double)*K*K);
		free(dQ);
	}
	//Parameter Estimation: pi and omega
	void P_i(double * p_user, double * pi_user) {
		double r = 0;
		for (int k = 0; k != K; ++k)
			r += p_user[k];
		for (int k = 0; k != K; ++k) {
			if (r == 0)
				pi_user[k] = 1.0 / K;
			else
				pi_user[k] = p_user[k] / r;
		}
	}
	//Parameter Estimation: p
	void P_xy() {
		for (int k1 = 0; k1 != K; ++k1) {
			for (int k2 = 0; k2 != K; ++k2) {
				if (Pxy2[k1*K + k2] == 0)
					Pxy[k1*K + k2] = 0.5;
				else Pxy[k1*K + k2] = Pxy1[k1*K + k2] / Pxy2[k1*K + k2];
			}
		}
	}
	//Parameter Estimation: A and B
	void Matrix_A(bool isuser) {
		double * xi_users = (isuser ? Xi_users : Xi_items);
		double * A = (isuser ? A_user : A_item);
		double xi_N_user = (isuser ? Xi_N_user : Xi_N_item);
		for (int k1 = 0; k1 != K; ++k1)for (int k2 = 0; k2 != K; ++k2) {
			A[k1*K + k2] = xi_users[k1*K + k2] / xi_N_user;
		}
	}
	//Algorithm: UpdateParameter
	void retrain(int userid, bool isuser, pair<int, double> * newrecord) {
		map<double, double* > & p = (isuser ? pu[userid] : pi[userid]);
		double* pi_users = (isuser ? Pi_users : Pi_items);
		vector< map<double, double* > > & q = (isuser ? pi : pu);
		const map<double, vector< pair<int, double> > > & ratings = (isuser ? ratings_user[userid] : ratings_item[userid]);
		double* A = (isuser ? A_user : A_item);
		value(NULL, 0, 0, (isuser ? 1 : 0));
		double* xi_user = (isuser ? Xi_user[userid] : Xi_item[userid]);
		double* xi_users = (isuser ? Xi_users : Xi_items);
		double* b0_user = (isuser ? B0_user[userid] : B0_item[userid]);

		map<double, double*> A_expm;
		for (map<double, vector< pair<int, double> > >::const_iterator it = ratings.begin(); it != ratings.end(); ++it) {
			map<double, vector< pair<int, double> > >::const_iterator it2 = it;
			it2++;
			if (it2 == ratings.end())
				continue;
			double dt = it2->first - it->first;
			if (!A_expm.count(dt)) {
				double * dQ = zeroMatrix(K, K);
				for (int k = 0; k != K * K; k++)
					dQ[k] = dt * A[k];
				A_expm[dt] = matrix_expm(dQ, K);
				free(dQ);
			}
		}

		map<double, double* > b;
		for (map<double, vector< pair<int, double> > >::const_iterator it = ratings.begin(); it != ratings.end(); ++it) {
			b[it->first] = (double*)malloc(sizeof(double)*K);
			for (int k = 0; k != K; ++k)
				b[it->first][k] = 1;
			const vector< pair<int, double> > & ratings1 = it->second;
			for (int n = 0; n != ratings1.size(); ++n) {
				for (int k1 = 0; k1 != K; ++k1) {
					double r = 0;
					for (int k2 = 0; k2 != K; ++k2)
						r += B((N - 1), (int)ratings1[n].second, value(Pxy, k1, k2))*q[ratings1[n].first][it->first][k2];
					b[it->first][k1] *= r;
				}
			}
			double r = 0;
			for (int k1 = 0; k1 != K; ++k1)
				r += b[it->first][k1];
			for (int k1 = 0; k1 != K; ++k1) {
				if (r == 0)
					b[it->first][k1] = 1.0 / K;
				else
					b[it->first][k1] /= r;
			}
		}

		map<double, double* > alpha;
		for (map<double, vector< pair<int, double> > >::const_iterator it = ratings.begin(); it != ratings.end(); ++it) {
			if (it == ratings.begin()) {
				alpha[it->first] = (double*)malloc(sizeof(double)*K);
				for (int k1 = 0; k1 != K; ++k1)
					alpha[it->first][k1] = pi_users[k1] * b[it->first][k1] * b0_user[k1];
			}
			map<double, vector< pair<int, double> > >::const_iterator it2 = it;
			it2++;
			if (it2 == ratings.end())
				continue;
			double* An = A_expm[it2->first - it->first];
			alpha[it2->first] = zeroMatrix(K, 1);
			for (int k2 = 0; k2 != K; ++k2) {
				for (int k1 = 0; k1 != K; ++k1) {
					alpha[it2->first][k2] += alpha[it->first][k1] * An[k1*K + k2];
				}
				alpha[it2->first][k2] *= b[it2->first][k2];
			}
		}

		map<double, double* > beta;
		for (map<double, vector< pair<int, double> > >::const_iterator it = ratings.begin(); it != ratings.end(); ++it) {
			beta[it->first] = zeroMatrix(K, 1);
		}
		for (map<double, double* >::reverse_iterator it = beta.rbegin(); it != beta.rend(); ++it) {
			if (it == beta.rbegin()) {
				for (int k1 = 0; k1 != K; ++k1)
					beta[it->first][k1] = 1;
			}
			map<double, double* >::reverse_iterator it1 = it;
			it1++;
			if (it1 == beta.rend())
				continue;
			double* An = A_expm[it->first - it1->first];
			for (int k1 = 0; k1 != K; ++k1) {
				for (int k2 = 0; k2 != K; ++k2) {
					beta[it1->first][k1] += beta[it->first][k2] * An[k1*K + k2] * b[it->first][k2];
				}
			}
		}

		map<double, double* > gamma;
		for (map<double, vector< pair<int, double> > >::const_iterator it = ratings.begin(); it != ratings.end(); ++it) {
			gamma[it->first] = zeroMatrix(K, 1);
			double r = 0;
			for (int k1 = 0; k1 != K; ++k1) {
				gamma[it->first][k1] = alpha[it->first][k1] * beta[it->first][k1];
				r += gamma[it->first][k1];
			}
			for (int k1 = 0; k1 != K; ++k1) {
				if (r == 0)
					gamma[it->first][k1] = 1.0 / K;
				else
					gamma[it->first][k1] /= r;
			}
		}

		map<double, double* > xi;
		for (map<double, vector< pair<int, double> > >::const_iterator it = ratings.begin(); it != ratings.end(); ++it) {
			map<double, vector< pair<int, double> > >::const_iterator it2 = it;
			it2++;
			if (it2 == ratings.end())
				continue;
			xi[it->first] = zeroMatrix(K, K);
			double* An = A_expm[it2->first - it->first];
			double r = 0;
			for (int k1 = 0; k1 != K; ++k1)for (int k2 = 0; k2 != K; ++k2) {
				xi[it->first][k1*K + k2] = alpha[it->first][k1] * An[k1*K + k2] * beta[it2->first][k2] * b[it2->first][k2];
				r += xi[it->first][k1*K + k2];
			}
			for (int k1 = 0; k1 != K; ++k1)for (int k2 = 0; k2 != K; ++k2) {
				if (r == 0)
					xi[it->first][k1*K + k2] = 1.0 / (K*K);
				else
					xi[it->first][k1*K + k2] /= r;
			}
		}

		double* p_users = (isuser ? P_users : P_items);
		for (int k = 0; k != K; ++k) {
			p_users[k] += gamma.begin()->second[k] - p[gamma.begin()->first][k];
		}
		P_i(p_users, pi_users);

		for (map<double, vector< pair<int, double> > >::const_iterator it = ratings.begin(); it != ratings.end(); ++it) {
			const vector< pair<int, double> > & ratings1 = it->second;
			for (int n = 0; n != ratings1.size(); ++n) {
				int item = ratings1[n].first;
				double rating = ratings1[n].second;
				if (newrecord != &ratings1[n]) {
					for (int k1 = 0; k1 != K; ++k1)for (int k2 = 0; k2 != K; ++k2) {
						value(Pxy1, k1, k2) -= p[it->first][k1] * q[item][it->first][k2] * (rating / (N - 1));
						value(Pxy2, k1, k2) -= p[it->first][k1] * q[item][it->first][k2];
					}
				}
				for (int k1 = 0; k1 != K; ++k1)for (int k2 = 0; k2 != K; ++k2) {
					value(Pxy1, k1, k2) += gamma[it->first][k1] * q[item][it->first][k2] * (rating / (N - 1));
					value(Pxy2, k1, k2) += gamma[it->first][k1] * q[item][it->first][k2];
				}
			}
		}
		P_xy();

		for (int k1 = 0; k1 != K; ++k1)for (int k2 = 0; k2 != K; ++k2) {
			xi_users[k1*K + k2] -= xi_user[k1*K + k2];
		}
		double * xi_P = zeroMatrix(K, K);
		double weights = 0;
		for (map<double, vector< pair<int, double> > >::const_iterator it = ratings.begin(); it != ratings.end(); ++it) {
			map<double, vector< pair<int, double> > >::const_iterator it2 = it;
			it2++;
			if (it2 == ratings.end())
				continue;
			double weight = it2->first - it->first;
			for (int k1 = 0; k1 != K; ++k1) {
				for (int k2 = 0; k2 != K; ++k2)
					xi_P[k1*K + k2] += xi[it->first][k1*K + k2] * weight;
			}
			weights += weight;

		}
		if (weights != 0) {
			for (int k1 = 0; k1 != K; ++k1)for (int k2 = 0; k2 != K; ++k2)
				xi_P[k1*K + k2] /= weights;
		}
		for (int k1 = 0; k1 != K; ++k1) {
			double r = 0;
			for (int k2 = 0; k2 != K; ++k2)
				r += xi_P[k1*K + k2];
			for (int k2 = 0; k2 != K; ++k2)
				xi_P[k1*K + k2] = lambda_2 * (r == 0 ? 1.0 / K : xi_P[k1*K + k2] / r) + (1 - lambda_2)*(k1 == k2 ? 1 : 0);
		}
		double * xi_Q = matrix_logm(xi_P, K);
		for (int k1 = 0; k1 != K; ++k1)for (int k2 = 0; k2 != K; ++k2)
			xi_user[k1*K + k2] = xi_Q[k1*K + k2];
		free(xi_Q);
		for (int k1 = 0; k1 != K; ++k1)for (int k2 = 0; k2 != K; ++k2) {
			xi_users[k1*K + k2] += xi_user[k1*K + k2];
		}
		Matrix_A(isuser);

		for (map<double, vector< pair<int, double> > >::const_iterator it = ratings.begin(); it != ratings.end(); ++it) {
			memcpy(p[it->first], gamma[it->first], sizeof(double)*K);
		}

		for (map<double, double* >::iterator it = A_expm.begin(); it != A_expm.end(); ++it)
			free(it->second);
		for (map<double, double* >::iterator it = b.begin(); it != b.end(); ++it)
			free(it->second);
		for (map<double, double* >::iterator it = alpha.begin(); it != alpha.end(); ++it)
			free(it->second);
		for (map<double, double* >::iterator it = beta.begin(); it != beta.end(); ++it)
			free(it->second);
		for (map<double, double* >::iterator it = gamma.begin(); it != gamma.end(); ++it)
			free(it->second);
		for (map<double, double* >::iterator it = xi.begin(); it != xi.end(); ++it)
			free(it->second);
	}
	
	//Algorithm: ConvertRating 
	void check_remove(bool isuser, int userid) {
		vector<int> & ratingnum = (isuser ? ratingnum_user : ratingnum_item);
		if (ratingnum[userid] <= M)
			return;

		map<double, vector< pair<int, double> > > & ratings = (isuser ? ratings_user[userid] : ratings_item[userid]);
		vector< map<double, vector< pair<int, double> > > > & ratings_q = (isuser ? ratings_item : ratings_user);
		map<double, double* > & p = (isuser ? pu[userid] : pi[userid]);
		double* p_users = (isuser ? P_users : P_items);
		double* pi_users = (isuser ? Pi_users : Pi_items);
		vector< map<double, double* > > & q = (isuser ? pi : pu);
		double* b0_user = (isuser ? B0_user[userid] : B0_item[userid]);
		double & xi_N_user = (isuser ? Xi_N_user : Xi_N_item);
		value(NULL, 0, 0, (isuser ? 1 : 0));

		double t0 = ratings.begin()->first;
		if (ratings[t0].size() == 1) {
			map<double, double* >::iterator it = p.find(t0);
			it++;
			for (int k = 0; k != K; ++k) {
				p_users[k] += it->second[k] - p[t0][k];
			}
			P_i(p_users, pi_users);
			xi_N_user -= 1;
		}
		int itemid = ratings[t0][0].first;
		double rating = ratings[t0][0].second;

		{
			for (int k1 = 0; k1 != K; k1++) {
				double r = 0;
				for (int k2 = 0; k2 != K; k2++)
					r += B((N - 1), (int)rating, value(Pxy, k1, k2))*q[itemid][t0][k2];
				b0_user[k1] *= r;
			}
			{
				double r = 0;
				for (int k = 0; k != K; k++)
					r += b0_user[k];
				for (int k = 0; k != K; k++)
					b0_user[k] = (r == 0 ? 1.0 : b0_user[k] / r);
			}
		}
		ratings[t0].erase(ratings[t0].begin());
		if (ratings[t0].size() == 0) {
			ratings.erase(ratings.begin());
		}
		ratingnum[userid]--;
	}
	//Algorithm: Update
	void update(int userid, int itemid, double rating, double t) {
		t /= T0;
		ratings_user[userid][t].push_back(make_pair(itemid, rating - 1));
		ratings_item[itemid][t].push_back(make_pair(userid, rating - 1));
		ratingnum_user[userid]++;
		ratingnum_item[itemid]++;

		if (!pu[userid].count(t)) {
			double dt = (pu[userid].empty() ? 0 : t - pu[userid].rbegin()->first);
			double * p = distribution(true, userid, t);
			pu[userid][t] = p;
			Xi_N_user += 1;
		}
		if (!pi[itemid].count(t)) {
			double dt = (pi[itemid].empty() ? 0 : t - pi[itemid].rbegin()->first);
			double * p = distribution(false, itemid, t);
			pi[itemid][t] = p;
			Xi_N_item += 1;
		}

		if (!is_classical_experiment) {
			check_remove(true, userid);
			check_remove(false, itemid);
		}

		retrain(userid, true, &(ratings_user[userid][t].back()));
		retrain(itemid, false, NULL);

	}
	//the posterior distribution of the random variable for userid at time t
	double * distribution(bool isuser, int userid, double t) {
		map<double, double* > & p_u = (isuser ? pu[userid] : pi[userid]);
		double* A = (isuser ? A_user : A_item);
		double* p_users = (isuser ? P_users : P_items);
		double * q = (double*)malloc(sizeof(double)*K);
		if (p_u.empty()) {
			initialize(q, isuser, userid);
			for (int k = 0; k != K; k++)
				p_users[k] += q[k];
			return q;
		}
		//
		double * p = NULL;
		double dt;
		double t0 = p_u.rbegin()->first;
		for (map<double, double* >::reverse_iterator it = p_u.rbegin(); it != p_u.rend(); it++) {
			t0 = it->first;
			p = it->second;
			if (t >= t0) {
				dt = t - t0;
				break;
			}
		}
		if (t <= t0) {
			memcpy(q, p, sizeof(double)*K);
			return q;
		}
		//
		double * dQ = zeroMatrix(K, K);
		for (int k = 0; k != K * K; k++)
			dQ[k] = dt * A[k];
		double * An = matrix_expm(dQ, K);
		free(dQ);
		memset(q, 0, sizeof(double)*K);
		double r = 0;
		for (int k2 = 0; k2 != K; k2++) {
			for (int k1 = 0; k1 != K; k1++) {
				q[k2] += p[k1] * An[k1*K + k2];
			}
			r += q[k2];
		}
		for (int k2 = 0; k2 != K; k2++) {
			if (r == 0)
				q[k2] = 1.0 / K;
			else
				q[k2] /= r;
		}
		free(An);
		return q;
	}
	//binomial distribution
	double B(int n, int k, double p) {
		double m = pow(p, k)*pow(1 - p, n - k);
		for (int l = k + 1; l <= n; l++)
			m *= l;
		for (int l = 1; l <= n - k; l++)
			m /= l;
		return m;
	}
	//Algorithm: Predicting (return the expectation of predicted rating)
	double calculate(int userid, int itemid, double t) {
		t /= T0;
		double * x = distribution(true, userid, t);
		double * y = distribution(false, itemid, t);
		double r = 1;
		for (int k1 = 0; k1 != K; k1++)for (int k2 = 0; k2 != K; k2++) {
			r += (N - 1)*Pxy[k1*K + k2] * x[k1] * y[k2];
		}
		free(x);
		free(y);
		return r;
	}
	//function in using Pxy
	double & value(double * p, int k1, int k2, int isuser = -1) {
		static bool is_user = true;
		static double n = 0;
		if (isuser != -1) {
			if (isuser == 0)
				is_user = false;
			else is_user = true;
			return n;
		}
		if (is_user)
			return p[k1*K + k2];
		else return p[k2*K + k1];
	}
	//function to find a rating
	double find_rating(bool isuser, int userid, int t0, int itemid) {
		map<double, vector< pair<int, double> > > & ratings = (isuser ? ratings_user[userid] : ratings_item[userid]);
		if (ratings.count(t0)) {
			for (int m = 0; m != ratings[t0].size(); m++) {
				if (ratings[t0][m].first == itemid) {
					return ratings[t0][m].second;
				}
			}
		}
		return -1;
	}

	int usern, itemn;							//Total number of users and items
	int K, N, M;								//hyper-parameters
	double lambda_1, lambda_2;
	double T0;
	double *Pxy;								//parameter p
	double *A_user, *A_item;					//parameter A and B
	double *Pi_users, *Pi_items;				//parameter pi and omega
	vector< map<double, double* > > pu, pi;		//parameter f and g
	vector< map<double, vector< pair<int, double> > > > ratings_user, ratings_item;	//recent ratings of each user and item
	vector<int> ratingnum_user, ratingnum_item;	//current number of recent ratings
	//program variables for Parameter Estimation
	double Xi_N_user, Xi_N_item;
	double *P_users, *P_items;
	double *Pxy1, *Pxy2, *Xi_users, *Xi_items;
	vector< double* > Xi_user, Xi_item, B0_user, B0_item;
	//
	bool is_classical_experiment;	// in classical experiment, the early and recent ratings setting are invalid
	double * initialize_vector;		// used to initialize

};



void CTHMM_mainfunction(vector<Rating> & ratings, int usern, int itemn){

	CTHMM model1(usern,itemn,false);
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
			r = model1.calculate(ratings[m].userid, ratings[m].itemid, ratings[m].time);
			rmse+=pow(ratings[m].rating-r,2);
			rmse1+=pow(ratings[m].rating-r,2);
			num1+=1;
			num_all+=1;
		}
		appeared_users.insert(ratings[m].userid);
		appeared_items.insert(ratings[m].itemid);
		model1.update(ratings[m].userid, ratings[m].itemid, ratings[m].rating, ratings[m].time);
	}
	cout << "Overall RMSE: " << endl;
	cout<<sqrt(rmse/num_all)<<endl;


}


void CTHMM_mainfunction_batch(vector<Rating> & rtraining, vector<Rating> & rtest, int usern, int itemn, bool classical){

	CTHMM model1(usern,itemn,classical);
	set<int> appeared_users, appeared_items;
	cout << "Training" << endl;
	for(int m=0;m!=rtraining.size();++m){
		if((rtraining.size()/100)>0 && (m+1)%(rtraining.size()/100)==0){
			cout<<".";
		}
		if((rtraining.size()/10)>0 && (m+1)%(rtraining.size()/10)==0){
			cout<<" "<<(m+1)/(rtraining.size()/10)*10<<"%"<<endl;
		}
		appeared_users.insert(rtraining[m].userid);
		appeared_items.insert(rtraining[m].itemid);
		model1.update(rtraining[m].userid, rtraining[m].itemid, rtraining[m].rating, rtraining[m].time);
	}
	cout << "Test" << endl;
	double rmse=0,num_all=0;
	for(int m=0;m!=rtest.size();++m){
		if(!appeared_users.count(rtest[m].userid) || !appeared_items.count(rtest[m].itemid))
			continue;
		double r = model1.calculate(rtest[m].userid, rtest[m].itemid, rtest[m].time);
		rmse+=pow(rtest[m].rating-r,2);
		num_all+=1;
	}
	
	cout<<num_all<<" ";
	cout << "RMSE: " ;
	cout<<sqrt(rmse/num_all)<<endl;

}
