#include <vector>
#include <random>
#include <algorithm>

using namespace std;

double l2_dist(vector<double> & a,vector<double> & b)
{
	double s = 0.0;
	for(int d=0;d<a.size();d++)
		s += (a[d]-b[d])*(a[d]-b[d]);
	return s;
}

vector<vector<double> > wsm_cpp(vector<vector<double> > & X, vector<double> & w, int K)
{
	
	int N = X.size();
	int D = X[0].size();

	//Random number generators
	random_device rd;
	mt19937 gen(rd());
    uniform_int_distribution<> randint;
	normal_distribution<>      randn  (0.0,1.0);
	
	//Pairwise distances between clusters and a permutation  
	vector<vector<double> > dist(K,vector<double>(K,0.0));
	vector<vector<int> >    perm(K,vector<int>   (K,  0));

	//Clusters
	vector<vector<double> > C(K,vector<double>(D,0.0));
	
	//First one's free.
	int rn = randint(gen)%N;
	for(int d=0;d<D;d++)
		C[0][d] = X[rn][d] + 1e-5*randn(gen);
	
	//Minmax (maxmin?) sampling... get the point in X with the farthest distance to C.
	for(int k=1;k<K;k++)
	{
		double k_dist = -1e15;
		int r = 0;
		for(int n=0;n<N;n++)
		{
			double r_dist = 1e15;
			for(int l=0;l<k;l++)
			{
				double a_dist = l2_dist(C[l],X[n]);
				if(a_dist < r_dist)
					r_dist = a_dist;
			}
			
			if(r_dist > k_dist)
			{
				k_dist = r_dist; 
				r = n;
			}
		}
		
		for(int d=0;d<D;d++)
			C[k][d] = X[r][d] + 1e-5*randn(gen);
	}
	
	vector<vector<double> > sum_wX(K,vector<double>(D,0.0));
	vector<double>          sum_w (K,0.0);

	vector<double> dist_S(N,1e15);

	vector<int> S(N,0);
	for(int n=0;n<N;n++)
	{
		double min_dist = 1e15;
		for(int k=0;k<K;k++)
		{
			double aux_dist = l2_dist(C[k],X[n]);
			if(min_dist > aux_dist)
			{
				min_dist = aux_dist;
				S[n] = k;
			}
		}
		dist_S[n] = min_dist;
	}

	double ssE_next = 0.0;
	double ssE_last = 1e15;
	for(int n=0;n<N;n++)
			ssE_next += dist_S[n];
	
	do
	{	
		//Calculate pairwise distances between clusters
		for(int k=0;k<K;k++)
			for(int l=k+1;l<K;l++)
			{
				dist[k][l] = l2_dist(C[k],C[l]);
				dist[l][k] = dist[k][l];
			}
		
		//Calculate the permutation from by increasing distances.
		for(int k=0;k<K;k++)
		{
			for(int l=0;l<K;l++)
				perm[k][l] = l;
				
			sort(perm[k].begin(),perm[k].end(),
				[&,k](int a,int b)->bool
				{
					return dist[k][a] < dist[k][b];
				});
		}
		
		//For each x in X find the closest c in C
		for(int n=0;n<N;n++)
		{
			int k = S[n];
			
			double min_dist = l2_dist(C[k],X[n]);
			double prev_dist = min_dist;
		
			for(int j=1;j<K;j++)
			{
				int l = perm[k][j];
				
				if(dist[k][l] >= 4.0*prev_dist)
					break;
				
				double aux_dist = l2_dist(C[l],X[n]);
				if(aux_dist <= min_dist)
				{
					min_dist = aux_dist;
					S[n] = l;
				}
			}
			dist_S[n] = min_dist;
			
		}
		
		//Clear sums
		for(int k=0;k<K;k++)
		{
			sum_w[k] = 0.0;
			for(int d=0;d<D;d++)
				sum_wX[k][d] = 0.0;
		}
		
		//Sum (w,x) depending on S
		for(int n=0;n<N;n++)
		{
			int k = S[n];
			
			sum_w[k] += w[n];
			for(int d=0;d<D;d++)
				sum_wX[k][d] += w[n]*X[n][d];
		}
		
		for(int k=0;k<K;k++)
			for(int d=0;d<D;d++)
				C[k][d] = sum_wX[k][d]/(1e-15+sum_w[k]);

		ssE_last = ssE_next;
		ssE_next = 0.0;		
		for(int n=0;n<N;n++)
			ssE_next += dist_S[n];
		
	}while((ssE_last - ssE_next)/(1e-15 + ssE_next) <= 1e-4);
	
	return C;
}
