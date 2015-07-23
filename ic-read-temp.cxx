#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "array2d.hpp"

#include "parameters.hpp"
#include "ic-read-temp.hpp"

void read_external_temperature_from_comsol(const Param &param,
                                           const Variables &var,
                                           double_vec &temperature)
{

    /* Read the eapsi_thermal file(x-coord y-coord temperature). The order of nodes in this file != the order of nodes_per_elem in EAPSI_Coord file.*/
    std::vector<double> xs, ys, zs, Ts;
    std::ifstream inputFile(param.ic.Temp_filename.c_str());
    std::string line;
    int i = 0;

    while(getline(inputFile, line)) {
        if (!line.length() || line[0] == '%')
            continue;
        std::istringstream iss(line);

	double x = 0., y = 0., z = 0., T = 0.;
	if (NDIMS ==3)
	  iss>>x>>y>>z>>T;
	else
	  iss>>x>>y>>T;
	 
	xs.push_back(x);
	ys.push_back(y);
	zs.push_back(z);
	Ts.push_back(T);

        ++i;
    }
    /* Read the EAPSI-Coord file(x-coord  y-coord). The order of nodes in this file != the order of nodes in eapsi_thermal file.*/ 
    std::vector<double>nxs, nys, nzs, nTs;
    std::vector<int>nis;
    std::ifstream ninputFile(param.ic.Nodes_filename.c_str());
    std::string nline;
    int ni = 0;

    while(getline(ninputFile, nline)) {
        if (!nline.length() || nline[0] == '#')
            continue;
        std::istringstream iss(nline);

        double x = 0., y = 0., z = 0.;
	if (NDIMS == 3)
          iss>>x>>y>>z;
	else
	  iss>>x>>y;
        
	nis.push_back(ni);
        nxs.push_back(x);
        nys.push_back(y);
	nzs.push_back(z);
	
        /* Assign temperature to the nodes according to the order in node-coord profile.*/
        int j = 0;
        for (j=0; j<xs.size();++j) {
            if (abs(nxs[ni]-xs[j])<0.001 && abs(nys[ni]-ys[j])<0.001 && abs(nzs[ni]-zs[j])<0.001)
                nTs.push_back(Ts[j]);
        }
        ++ni;
    }

    /* Write x&y-coord to array_t coord(nnodes). Write temperature to double_vec temperature(nnodes).*/
    int n;
    int nnodes = nis.size();
    array_t coord(nnodes);	// coord[node#][dim#];
    double_vec inputtemperature(nnodes);

    for (n=0; n<nis.size(); ++n) {
  	coord[n][0] = nxs[n];
   	coord[n][1] = nys[n];
#ifdef THREED
	coord[n][2] = nzs[n];
#endif
   	inputtemperature[n] = nTs[n];
    }

    /* Read the EAPSI_Connectivity file(n0 n1 n2).*/
    std::vector<int> es, n0s, n1s, n2s, n3s;
    std::ifstream einputFile(param.ic.Connectivity_filename.c_str());
    std::string eline;
    int e = 0;

    while(getline(einputFile, eline)) {
        if (!eline.length() || eline[0] == '#')
            continue;
        std::istringstream iss(eline);

        int n0, n1, n2, n3 = 0;
	if (NDIMS == 3)
          iss>>n0>>n1>>n2>>n3;
	else
	  iss>>n0>>n1>>n2;
          n0s.push_back(n0);
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
  	connectivity[m][0] = n0s[m];
   	connectivity[m][1] = n1s[m];
   	connectivity[m][2] = n2s[m];
#ifdef THREED
	connectivity[m][3] = n3s[m];
#endif
    }

    if (1) {
    /* checking */
    std::cout << "# of nodes: " << coord.size() << '\n';
    std::cout << "# of elem:  " << connectivity.size() << '\n';
    for (m=0; m<5; ++m) {
   	//int u = rand() % nis.size()/5+m*7000;
   	int u = m*70;
   	std::cout<<"The Temp @ point(node # "<<nis[u]<<"): ("<<coord[u][0]<<", "<<coord[u][1]<<") is "<<inputtemperature[u]<<"C."<<std::endl;
    }
    std::cout<<"There are "<<nxs.size()<<" of nodes."<<std::endl;
    for (m=0; m<7; ++m) {
   	int v = rand() % es.size()/7+m*10000;
   	std::cout<<"The nodes for element "<<es[v]<<" are: "<<connectivity[v][0]<<", "<<connectivity[v][1]<<" and "<<connectivity[v][2]<<"."<<std::endl;
    }

    std::cout<<"There are "<<es.size()<<" of elements."<<std::endl;
    }
    std::exit(1);
    return;
}
