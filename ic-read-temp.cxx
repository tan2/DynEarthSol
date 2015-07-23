#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "array2d.hpp"

const int NDIMS = 2, NODES_PER_ELEM = 3;
typedef std::vector<double> double_vec;
typedef Array2D<double,NDIMS> array_t;
typedef Array2D<int,NODES_PER_ELEM> conn_t;

int main () {

    /* Read the eapsi_thermal file(x-coord y-coord temperature). The order of nodes in this file != the order of nodes_per_elem in EAPSI_Coord file.*/
    std::vector<double> xs, ys, Ts;
    std::vector<int> is;
    std::ifstream inputFile("/home/dt/Dropbox/Ubuntu/Sumatra_Subduction_boundary_positions/EAPSI_MAY_2015/eapsi_thermal_HD.dat"); 
    std::string line;
    int i = 0;

    while(getline(inputFile, line)) {
        if (!line.length() || line[0] == '%')
            continue;
        std::istringstream iss(line);

        double x = 0., y = 0., T = 0.;
        iss>>x>>y>>T;
        is.push_back(i);
        xs.push_back(x);
        ys.push_back(y);
        Ts.push_back(T);

        ++i;
    }
    /* Read the EAPSI-Coord file(x-coord  y-coord). The order of nodes in this file != the order of nodes in eapsi_thermal file.*/ 
    std::vector<double>nxs, nys, nTs;
    std::vector<int>nis;
    std::ifstream ninputFile("/home/dt/Dropbox/Ubuntu/Sumatra_Subduction_boundary_positions/EAPSI_MAY_2015/EAPSI_Coord_HD.dat");
    std::string nline;
    int ni = 0;

    while(getline(ninputFile, nline)) {
        if (!nline.length() || nline[0] == '#')
            continue;
        std::istringstream iss(nline);

        double x = 0., y = 0.;
        iss>>x>>y;
        nis.push_back(ni);
        nxs.push_back(x);
        nys.push_back(y);

        /* Assign temperature to the nodes according to the order in node-coord profile.*/
        int j = 0;
        for (j=0; j<is.size();++j) {
            if (abs(nxs[ni]-xs[j])<0.001 && abs(nys[ni]-ys[j])<0.001)
                nTs.push_back(Ts[j]);
        }
        ++ni;
    }

    /* Write x&y-coord to array_t coord(nnodes). Write temperature to double_vec temperature(nnodes).*/
    int n;
    int nnodes = nis.size();
    array_t coord(nnodes);	// coord[node#][dim#];
    double_vec temperature(nnodes);

    for (n=0; n<nis.size(); ++n) {
  	coord[n][0] = nxs[n];
   	coord[n][1] = nys[n];
   	temperature[n] = nTs[n];
    }

    /* Read the EAPSI_Connectivity file(n0 n1 n2).*/
    std::vector<long> es, n1s, n2s, n3s;
    std::ifstream einputFile("/home/dt/Dropbox/Ubuntu/Sumatra_Subduction_boundary_positions/EAPSI_MAY_2015/EAPSI_Connectivity_HD.dat"); 
    std::string eline;
    long e = 0;

    while(getline(einputFile, eline)) {
        if (!eline.length() || eline[0] == '#')
            continue;
        std::istringstream iss(eline);

        long n1, n2, n3;
        iss>>n1>>n2>>n3;
        n1s.push_back(n1);
        n2s.push_back(n2);
        n3s.push_back(n3);
        es.push_back(e);

        ++e;
    }

    /* Write nodes to conn_t connectivity(nelem).*/
    int m;
    int nelem = es.size();
    conn_t connectivity(nelem);	// connectivity[elem#][0-NODES_PER_ELEM-1]

    for (m=0; m<es.size(); ++m) {
  	connectivity[m][0] = n1s[m];
   	connectivity[m][1] = n2s[m];
   	connectivity[m][2] = n3s[m];
    }


    /* checking */
    for (m=0; m<5; ++m) {
   	int u = rand() % nis.size()/5+m*7000;
   	std::cout<<"The Temp @ point(node # "<<nis[u]<<"): ("<<coord[u][0]<<", "<<coord[u][1]<<") is "<<temperature[u]<<"C."<<std::endl;
    }
    std::cout<<"There are "<<nxs.size()<<" of nodes."<<std::endl;
    for (m=0; m<7; ++m) {
   	int v = rand() % es.size()/7+m*10000;
   	std::cout<<"The nodes for element "<<es[v]<<" are: "<<connectivity[v][0]<<", "<<connectivity[v][1]<<" and "<<connectivity[v][2]<<"."<<std::endl;
    }

    std::cout<<"There are "<<es.size()<<" of elements."<<std::endl;
    return 0;
}
