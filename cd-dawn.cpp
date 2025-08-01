#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <set>
#include <fstream>
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <ios>
#include <string>
#include <set>
#include <map>
#include <list>
#include <numeric>
#include<unordered_map>
using namespace std;
class Graph
    {
    public:
        set<int> V; //only vertices with positive degree are stored
        map<int,set<int>>in_neighbors, out_neighbors,neighbors;
        map<int, float> indegree,outdegree,Centnode1,Centnode;
        map<int, map<int, float> > weight;
        void read_edgelist(string&, bool, bool);
        inline int order(){ return V.size(); }
        int ecount();
        float get_weight(int, int);
        void print_graph();
        float proximity(int, int);
        float get_sn_nbrs_proximity(int, map<int, float>&);
        friend float proximity(Graph&, int, int);
        friend float overlap(Graph&, set<int>&, set<int>&);
        friend float overlapping_weighted_modularity(Graph&, map<int, set<int> >&, int, map<int, set<int> >& );
    };
//other function declarations
void merge_communities(Graph&, map<int, set<int> >&, map<int, set<int> >&, float, int);
void write_seeds(set<int>&, string&);
float overlap(Graph& g, set<int>&, set<int>&);
void usage(void);
void read_communities(string& comfile, map<int, set<int> >&, map<int, set<int> >&);
void print_communities(map<int, set<int> >&);
float overlapping_weighted_modularity(Graph&, map<int, set<int> >&, map<int, set<int> >& );
void print_set(set<int>&);
void print_vector(vector<int>& );
map<int, set<int> > coms;
map<int, set<int> > memberships;



void Graph::read_edgelist(string& edgefile, bool weighted, bool directed)
{
    ifstream fin;
    fin.open(edgefile);
    if(!fin.is_open())
    {
        cout<<"The file containing the edgelist could not be opened."<<endl;
        exit(1);
    }
    string line;
    while ( getline(fin, line ) )
    {
        istringstream is( line );
        int u, v;
        is>>u;
        is>>v;

        if(u == v || weight[u].find(v) != weight[u].end())
            continue;

        float w;
        if(weighted == true)
        {
            if(is.str() == " " || is.eof())
            {   cout<<endl<<"The edge list has missing weights."<<endl;
                exit(1);
            }
            is>>w;
        }
        else
            w = 1;
       if (directed==true)
        {
             weight[u][v] = w;
            if(outdegree.find(u) == outdegree.end())
                outdegree[u] = w;
            else
                 outdegree[u] = outdegree[u]+w;
            if(indegree.find(v) == indegree.end())
                indegree[v] = w;
            else
                 indegree[v] = indegree[v]+w;
        }

       in_neighbors[v].insert(u);
       out_neighbors[u].insert(v);
       neighbors[u].insert(v);
       neighbors[v].insert(u);
        V.insert(u);
        V.insert(v);
    }
 for(auto st=V.begin();st!=V.end();++st)
    {
        //cout<<"strength of node "<<*st<<" = "<< indegree[*st]+outdegree[*st]<<endl;
        if((indegree[*st]+outdegree[*st])>0)
            {
              int i=*st;
              float s=0;
              for(auto si=neighbors[i].begin();si!=neighbors[i].end();++si)
                 {
                   for(auto sj=neighbors[i].begin();sj!=neighbors[i].end();++sj)
                    {
                       if(si!=sj && in_neighbors[*si].find(*sj)!=in_neighbors[*si].end())
                       s=s+weight[*sj][*si];
                       if(si!=sj && out_neighbors[*si].find(*sj)!=out_neighbors[*si].end())
                        s=s+weight[*si][*sj];
                    }
                  }
              Centnode.insert({i,(indegree[i]+outdegree[i])+(s/2)});
           }
     }

}

int Graph::ecount()
{
    int degree_sum = 0;
     for(auto mi = weight.begin(); mi != weight.end(); ++mi)
         {
           for(auto ik = mi->second.begin(); ik != mi->second.end(); ++ik)
               if(weight[mi->first][ik->first] != 0)
                  degree_sum++;
          }
    return degree_sum;
}

void print_set(set<int>& s)
{
    set<int>::iterator sitr1;
    for(sitr1 = s.begin(); sitr1 != s.end(); ++sitr1)
            cout<<*sitr1<<" ";
    cout<<endl;
}

float Graph::get_weight(int u, int v)
{
    if(weight[u].find(v) != weight[u].end())
        return weight[u][v];
    else
        return 0;
}


void Graph::print_graph()
{
    cout<<"No. of vertices: "<<order()<<endl;
    cout<<"No. of edges: "<<ecount()<<endl;
    cout<<"\n"<<"vertex"<<setw(8)<<"degree"<<endl;

    cout<<"\n"<<"node"<<setw(10)<<"nbrs"<<setw(10)<<"weight"<<endl;
    map<int, map<int, float> >::iterator it;
    for(it=weight.begin(); it != weight.end(); ++it)
    {
        cout<<it->first<<setw(5)<<"--- "<<endl;
        map<int, float>::iterator jt;
        for(jt = it->second.begin(); jt != it->second.end(); ++jt)
            cout<<setw(10)<<jt->first<<setw(10)<<jt->second<<endl;
    }

     cout<<endl<<"centnode"<<endl;
   map<int,float>::iterator ih;
   for(ih=Centnode.begin();ih!=Centnode.end();++ih){
    cout<<ih->first<<"    "<<ih->second<<endl;
    }
    cout<<"centnode size"<<Centnode.size()<<endl;

/*cout<<"in_neighbors"<<endl;

     for (auto it = in_neighbors.begin(); it != in_neighbors.end(); ++it)
     {
         cout<<it->first<<setw(5);
        print_set(in_neighbors[it->first]);
      }
    cout<<"out_neighbors"<<endl;

     for (auto it = out_neighbors.begin(); it != out_neighbors.end(); ++it)
     {
         cout<<it->first<<setw(5)<<"===";
        print_set(out_neighbors[it->first]);
      }
/*cout<<"indegree "<<endl;
      for(auto it=indegree.begin();it!=indegree.end();++it)
      {
          cout<<it->first<<"--"<<it->second<<endl;
      }
cout<<"outdegree "<<endl;
      for(auto it=outdegree.begin();it!=outdegree.end();++it)
      {
          cout<<it->first<<"--"<<it->second<<endl;
      }*/


}



float Graph::proximity(int u, int v)
{
    float prox;
    float weight_sum1 = 0;
    float weight_sum2 = 0;
    float common_nbrs_strength1 = 0, common_nbrs_strength2=0;
    float A = 0, a = 0, B=0, b=0;

   // First loop for node u
    for (auto i = out_neighbors[u].begin(); i!= out_neighbors[u].end();++i)
         {
         if (in_neighbors[v].find(*i) !=in_neighbors[v].end())
            {
             weight_sum1 += weight[u][*i];
             weight_sum2 += weight[*i][v];
             common_nbrs_strength1 += outdegree[*i];
             }
        }

    if (common_nbrs_strength1 != 0 && out_neighbors[u].find(v)!=out_neighbors[u].end())
        A = 0.5 * (((weight_sum1 + weight[u][v]) / outdegree[u]) + (weight_sum2 / common_nbrs_strength1));
            else if (common_nbrs_strength1 != 0 )
                A = 0.5 * (((weight_sum1) / outdegree[u]) + (weight_sum2 / common_nbrs_strength1));
                    else if (common_nbrs_strength1==0 && out_neighbors[u].find(v)!=out_neighbors[u].end())
                            a = weight[u][v] / outdegree[u];


   weight_sum1 = 0;
    weight_sum2 = 0;

    // Second loop for node v
    for (auto i = out_neighbors[v].begin(); i!= out_neighbors[v].end();++i)
         {
         if (in_neighbors[u].find(*i) !=in_neighbors[u].end())
            {
             weight_sum1 += weight[v][*i];
             weight_sum2 += weight[*i][u];
             common_nbrs_strength2 += outdegree[*i];
             }
        }

    if (common_nbrs_strength2 != 0 && out_neighbors[v].find(u)!=out_neighbors[v].end())
        B = 0.5 * (((weight_sum1 + weight[v][u]) / outdegree[v]) + (weight_sum2 / common_nbrs_strength2));
            else if (common_nbrs_strength2 != 0 )
                B = 0.5 * (((weight_sum1) / outdegree[v]) + (weight_sum2 / common_nbrs_strength2));
                    else if (common_nbrs_strength2==0 && out_neighbors[v].find(u)!=out_neighbors[v].end())
                            b = weight[v][u] / outdegree[v];


    if(common_nbrs_strength1==0 && common_nbrs_strength2==0)
      prox=max(a,b);
      else
       prox=max(A,B);


    return prox;

}

float Graph::get_sn_nbrs_proximity(int u, map<int, float>& prox)
{
    //set<int>::iterator mi, mj;
    float max_prox = 0;
   // cout<<"size of neighbour u == "<<outdegree[u].size()<<endl;
    for(auto mi = neighbors[u].begin(); mi != neighbors[u].end(); ++mi)
    {
        //cout<<"node prox 1 ="<<*mi<<endl;
        if(prox.find(*mi) == prox.end())
        {
            prox[*mi] = proximity(u, *mi);
                 if(prox[*mi] > max_prox)
                  max_prox = prox[*mi];
        }

      for(auto mj = neighbors[*mi].begin(); mj != neighbors[*mi].end(); ++mj)
            if(*mj!= u && prox.find(*mj) == prox.end())
            {
                 //cout<<"node prox 2 ="<<*mj<<endl;
                prox[*mj] = proximity(u, *mj);
                if(prox[*mj] > max_prox)
                    max_prox = prox[*mj];
            }
    }
    //cout<<"max prox = "<<max_prox<<endl;
    return max_prox;
}

void merge_communities(Graph& g, map<int, set<int> >& coms, map<int, set<int> >& members, float given_max_ov, int minc)
{
    //merges communities with overlap >= given_max_ov
    map<int, set<int> >::iterator ci, cj;
    set<int>::iterator si, sj, sk;
    map<int, float>::iterator mi;
    set<int> labels;
    for(ci = coms.begin(); ci != coms.end(); ++ci)
        labels.insert(ci->first);
    do
    {
        //cout<<"all coms"<<endl; print_community(coms);
        si = labels.begin();
        int r = rand()%labels.size();
        for(int i = 1; i<=r; ++i)
            ++si;
        //cout<<"curr com"<<endl<<*si<<"->"; print_set(coms[*si]);
        set<int> neighboring_com_labels;
        for(sj = coms[*si].begin(); sj != coms[*si].end(); ++sj)
            for(mi = g.weight[*sj].begin(); mi != g.weight[*sj].end(); ++mi)
                if(coms[*si].find(mi->first) == coms[*si].end())
                {
                    //cout<<"members["<<*sk<<"]->"; print_set(members[*sk]);
                    for(sk = members[mi->first].begin(); sk != members[mi->first].end(); ++sk)
                        if(coms[*sk].size() >= coms[*si].size())  //neighboring coms that are bigger than coms[*si]
                            neighboring_com_labels.insert(*sk);
                }
        //cout<<"neighboring coms "<<endl; print_set(neighboring_com_labels);
        float max_ov = 0;
        set<int> high_overlapping_coms;
        for(sj = neighboring_com_labels.begin(); sj != neighboring_com_labels.end(); ++sj)
        {
            float curr_ov = overlap(g, coms[*si], coms[*sj]);
            //cout<<"curr ov = "<<curr_ov<<endl;
            if(coms[*si].size() < minc)
            {
                //cout<<"This comm has size smaller than min"<<endl;
                if(curr_ov > max_ov)
                {
                    max_ov = curr_ov;
                    high_overlapping_coms.clear();
                    high_overlapping_coms.insert(*sj);
                }
                else
                    if(curr_ov == max_ov)
                        high_overlapping_coms.insert(*sj);
            }
            else
            {
                if( curr_ov >= given_max_ov)
                   high_overlapping_coms.insert(*sj);
            }
        }
        /*cout<<"high overlapping coms"<<endl;
        for(sj = high_overlapping_coms.begin(); sj != high_overlapping_coms.end(); ++sj)
        {
            cout<<*sj<<"->";
            print_set(coms[*sj]);
        } */
        for(sj = high_overlapping_coms.begin(); sj != high_overlapping_coms.end(); ++sj)
        {   //cout<<"curr ov = "<<curr_ov<<endl;
            for(sk = coms[*si].begin(); sk != coms[*si].end(); ++sk)
            {
                coms[*sj].insert(*sk);
                members[*sk].erase(*si);
                members[*sk].insert(*sj);
            }
        }
        /* cout<<"modified coms"<<endl;
        for(sj = high_overlapping_coms.begin(); sj != high_overlapping_coms.end(); ++sj)
        {
            cout<<*sj<<"->";
            print_set(coms[*sj]);
        } */
        if(!high_overlapping_coms.empty())
            coms.erase(*si);
        labels.erase(*si);
    }while(!labels.empty());
}

float overlap(Graph& g, set<int>& C1, set<int>& C2)
{

    float intercom_weight = 0, deg1 = 0, deg2 = 0;

    for(auto si = C1.begin(); si != C1.end(); ++si)
    {
        for(auto sj = g.out_neighbors[*si].begin(); sj != g.out_neighbors[*si].end(); ++sj)
            if(C2.find(*sj)!=C2.end())
                intercom_weight += g.weight[*si][*sj];
        for(auto sj = g.in_neighbors[*si].begin(); sj != g.in_neighbors[*si].end(); ++sj)
            if((C2.find(*sj)!= C2.end()))
                intercom_weight += g.weight[*sj][*si];
         deg1 += g.outdegree[*si] + g.indegree[*si];

    }
    for(auto sj = C2.begin(); sj != C2.end(); ++sj)
        deg2 += g.outdegree[*sj]+g.indegree[*sj];
    return intercom_weight/min(deg1, deg2);


}




void read_communities(string& comfile, map<int, set<int> >& coms, map<int, set<int> >& memberships)
{
    std::ifstream fin(comfile);
    string line;
    if(!fin.is_open())
    {   cout<<"Community file could not be opened."<<endl;
        exit(1);
    }
    int i = 1;
    while ( std::getline(fin, line ) )
    {
        istringstream is( line );
        int v;
        while(is>>v)
        {
            coms[i].insert(v);
            memberships[v].insert(i);
        }
        i++;
   }
}

pair<int, int>LargestValue(map<int, float> sampleMap)
{
	pair<int, float> MaxValue = *sampleMap.begin();
	map<int, float>::iterator currentEntry;
	for (currentEntry = sampleMap.begin();currentEntry != sampleMap.end();++currentEntry) {
		if (currentEntry->second> MaxValue.second) {
             MaxValue= make_pair(currentEntry->first,currentEntry->second);
		   }
	    }


    return MaxValue;
}
void print_communities(map<int, set<int> >& C)
{
    cout<<"community"<<endl;
    map<int, set<int> >::iterator sitr;
    for(sitr = C.begin(); sitr != C.end(); ++sitr)
    {
        copy((*sitr).second.begin(), (*sitr).second.end(), ostream_iterator<int>(cout, " "));
        cout<<endl;
    }
}




int main(int argc, char* argv[])
{
    std::chrono::time_point<std::chrono::system_clock> start_time, end_time;
    start_time = std::chrono::system_clock::now();  //time starts
    float max_overlap = 0.40, rho_0 = 0.30;
    bool weighted = false;
    bool directed =true;
    Graph g;
    map<int, set<int> > coms;
    string network_file;

     if(argc < 2 || argc > 7)
     {
	cout<<"Please see README file."<<endl;
	exit(1);
     }
    if(argv[1][0] == '-')
    {
	cout<<"Please see README file."<<endl;
	exit(1);
     }
    else
    {
        network_file = string(argv[1]);
    }
    int i=2;
    while(i < argc)
    {
        string arg = string(argv[i]);
        if(arg == "-w")
        {
            weighted = true;
            i++;
        }
        else
            if(arg == "-ov")
            {
                istringstream is(argv[i+1]);
                is>>max_overlap;
                if( max_overlap < 0 || max_overlap > 0.5)
		{
		    cout<<"Please see README file."<<endl;
		    exit(1);
		}
                i += 2;
            }
            else
                if(arg == "-rh")
                {
                    istringstream is(argv[i+1]);
                    is>>rho_0;
                    if( rho_0 < 0 || rho_0 >= 1)
		    {    
			cout<<"Please see README file."<<endl;
			exit(1);
		    }
                    i += 2;
                }
                else
                    {
		    cout<<"Please see README file."<<endl;
		    exit(1);
		    }
}
    g.read_edgelist(network_file,weighted, directed);
    cout<<"-------------------------------------------------------------------"<<endl;
    cout<<"Community Detection in Directed And Weighted Networks (CD-DAWN) "<<endl;
    cout<<"Authors: Abhinav Kumar, Pawan Kumar and Ravins Dohare"<<endl;
    cout<<"Email: abhinavkumar080395@gmail.com, pkumariitd@gmail.com, ravinsdohare@gmail.com"<<endl;
    cout<<"-------------------------------------------------------------------"<<endl;
 
 
    map<int, set<int> >::iterator ci;
    set<int> uncovered = g.V;
    //cout<<"initializing communities..."<<endl;
    set<int> seeds;
    map<int, set<int> > membership;
    set<int>::iterator si, sj;
    //cout<<"expanding communities..."<<endl;
    int com_count = 1;
      while(!g.Centnode.empty())
    {
        //cout<<"g.centnode size = "<<g.Centnode.size()<<endl;
        pair<int,int> yup=LargestValue(g.Centnode);
        int r=yup.first;
        //cout<<"curr uncovered = "<<r<<endl;
        map<int, float> prox;
        float max_prox;
        max_prox = g.get_sn_nbrs_proximity(r, prox);
        map<int, float>::iterator mi;
        //cout<<"nbrs proximities"<<endl;
       //for(mi = prox.begin(); mi != prox.end(); ++mi)
       // cout<<mi->first<<"--"<<mi->second<<endl;
       // cout<<"max prox = "<<max_prox<<endl;
        for(mi = prox.begin(); mi != prox.end(); ++mi)
        {
            if(mi->second > rho_0 || mi->second >= max_prox)
            {
                membership[r].insert(membership[mi->first].begin(), membership[mi->first].end());
                g.Centnode.erase(mi->first);

            }
        }
        if(membership[r].empty())
        {
            membership[r].insert(com_count);
            com_count++;
        }
        for(mi = prox.begin(); mi != prox.end(); ++mi)
             if(mi->second > rho_0 && membership[mi->first].empty())
             {
                 membership[mi->first].insert(membership[r].begin(), membership[r].end());

             }
          g.Centnode.erase(r);
     }


      for(si = g.V.begin(); si != g.V.end(); ++si)
    {
        for(sj = membership[*si].begin(); sj != membership[*si].end(); ++sj)
            coms[*sj].insert(*si);
    }

    //cout<<"communities before  merging = "<<coms.size()<<endl;
    int avg_csize = 0;
    for(ci = coms.begin(); ci != coms.end(); ++ci)
        {   //cout<<ci->second.size()<<", ";
            avg_csize += ci->second.size();
        }
    avg_csize /= coms.size();
    //cout<<"average size of initial communities = "<<avg_csize<<endl;
    //cout<<"merging communities..."<<endl;
   merge_communities(g, coms, membership, max_overlap, avg_csize);
    //cout<<"writing final communities..."<<endl;
    cout<<endl;
   ostringstream str_o;
    str_o.setf(ios::fixed, ios::floatfield);
    str_o.precision(2);
    str_o << max_overlap;

    ostringstream comfile;
    comfile<<"./cd-dawn-coms.txt";
    ofstream fout(comfile.str());
    if(!fout.is_open())
    {
        cout<<"Destination file for communities could not be opened.";
        exit(1);
    }
    for(ci = coms.begin(); ci != coms.end(); ++ci)
    {
        copy(ci->second.begin(), ci->second.end(), ostream_iterator<int>(fout, " "));
        fout<<endl;
    }
    fout.close();
      end_time = std::chrono::system_clock::now();  //time ends
    std::chrono::duration<double> elapsed_seconds = end_time-start_time;
    cout.setf(ios::left, ios::adjustfield);
    cout<<"-------------------------------------------------------------------"<<endl;
    cout<<setw(33)<<"Network file"<<"= "<<setw(20)<<network_file<<endl;
    cout<<setw(33)<<"Network order"<<"= "<<setw(20)<<g.order()<<endl;
    cout<<setw(33)<<"No. of edges"<<"= "<<setw(20)<<g.ecount()<<endl;
    cout<<setw(33)<<"Total non-singleton communities"<<"= "<<coms.size()<<endl;
    cout<<setw(33)<<"Community file"<<"= "<<setw(20)<<comfile.str()<<endl;
    cout<<setw(33)<<"Time elapsed"<<"= "<<elapsed_seconds.count()<<"s\n";
    cout<<"-------------------------------------------------------------------"<<endl;
    return 0;
}
