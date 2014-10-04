/****************************************************************************
Logiciel initial: MMG3D Version 4.0
Co-auteurs : Cecile Dobrzynski et Pascal Frey.
Propriétaires :IPB - UPMC -INRIA.

Copyright © 2004-2005-2006-2007-2008-2009-2010-2011,
diffusé sous les termes et conditions de la licence publique générale de GNU
Version 3 ou toute version ultérieure.

Ce fichier est une partie de MMG3D.
MMG3D est un logiciel libre ; vous pouvez le redistribuer et/ou le modifier
suivant les termes de la licence publique générale de GNU
Version 3 ou toute version ultérieure.
MMG3D est distribué dans l'espoir qu'il sera utile, mais SANS
AUCUNE GARANTIE ; sans même garantie de valeur marchande.
Voir la licence publique générale de GNU pour plus de détails.
MMG3D est diffusé en espérant qu’il sera utile,
mais SANS AUCUNE GARANTIE, ni explicite ni implicite,
y compris les garanties de commercialisation ou
d’adaptation dans un but spécifique.
Reportez-vous à la licence publique générale de GNU pour plus de détails.
Vous devez avoir reçu une copie de la licence publique générale de GNU
en même temps que ce document.
Si ce n’est pas le cas, aller voir <http://www.gnu.org/licenses/>.
/****************************************************************************
Initial software: MMG3D Version 4.0
Co-authors: Cecile Dobrzynski et Pascal Frey.
Owners: IPB - UPMC -INRIA.

Copyright © 2004-2005-2006-2007-2008-2009-2010-2011,
spread under the terms and conditions of the license GNU General Public License
as published Version 3, or (at your option) any later version.

This file is part of MMG3D
MMG3D is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.
MMG3D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with MMG3D. If not, see <http://www.gnu.org/licenses/>.
****************************************************************************/
#include "mesh.h"

int MMG_cutadd(pMesh mesh,pHedge hed,int icas,int k,int p0,int p1,int p2,int p3,int p4,int p5, int ref) {
  pTetra                pt;
  pPoint                ppt,pp0,pp1,pp2,pp3,pp4,pp5;

  pp0 = &mesh->point[p0];
  pp1 = &mesh->point[p1];
  pp2 = &mesh->point[p2];
  pp3 = &mesh->point[p3];
  pp4 = &mesh->point[p4];
  pp5 = &mesh->point[p5];
  mesh->np++;
  ppt = &mesh->point[mesh->np];
  ppt->c[2] = (1./6.)*(pp0->c[2]+pp1->c[2]+pp2->c[2]+pp3->c[2]+pp4->c[2]+pp5->c[2]);
  ppt->c[1] = (1./6.)*(pp0->c[1]+pp1->c[1]+pp2->c[1]+pp3->c[1]+pp4->c[1]+pp5->c[1]);
  ppt->c[0] = (1./6.)*(pp0->c[0]+pp1->c[0]+pp2->c[0]+pp3->c[0]+pp4->c[0]+pp5->c[0]);
  ppt->ref = pp0->ref;
  //printf("prisme : %d %d %d %d %d %d + %d\n",p0,p1,p2,p3,p4,p5,mesh->np);
  if(icas & 1) {
    //printf("icas %d --> 0 ?\n",icas);
    pt = &mesh->tetra[k + 1];
    pt->v[0] = p0;
    pt->v[1] = p4;
    pt->v[2] = p3;
    pt->v[3] = mesh->np;
    pt->ref  = ref;
    pt = &mesh->tetra[k + 2];
    pt->v[0] = p0;
    pt->v[1] = p1;
    pt->v[2] = p4;
    pt->v[3] = mesh->np;
    pt->ref  = ref;
  } else {
    if(icas & 4) {
      //printf("icas %d --> 2 ?\n",icas);

    } else {
      MMG_edgePut(hed,p1,p3,2);
    }
    pt = &mesh->tetra[k + 1];
    pt->v[0] = p0;
    pt->v[1] = p1;
    pt->v[2] = p3;
    pt->v[3] = mesh->np;
    pt->ref  = ref;
    pt = &mesh->tetra[k + 2];
    pt->v[0] = p1;
    pt->v[1] = p4;
    pt->v[2] = p3;
    pt->v[3] = mesh->np;
    pt->ref  = ref;
  }

  if(icas & 8) {
    //printf("icas %d --> 3 ?\n",icas);
    pt = &mesh->tetra[k + 3];
    pt->v[0] = p1;
    pt->v[1] = p2;
    pt->v[2] = p5;
    pt->v[3] = mesh->np;
    pt->ref  = ref;
    pt = &mesh->tetra[k + 4];
    pt->v[0] = p1;
    pt->v[1] = p5;
    pt->v[2] = p4;
    pt->v[3] = mesh->np;
    pt->ref  = ref;
  } else {
    if(icas & 32) {
      //printf("icas %d --> 5 ?\n",icas);

    } else {
      MMG_edgePut(hed,p2,p4,2);
    }
    pt = &mesh->tetra[k + 3];
    pt->v[0] = p1;
    pt->v[1] = p2;
    pt->v[2] = p4;
    pt->v[3] = mesh->np;
    pt->ref  = ref;
    pt = &mesh->tetra[k + 4];
    pt->v[0] = p4;
    pt->v[1] = p2;
    pt->v[2] = p5;
    pt->v[3] = mesh->np;
    pt->ref  = ref;
  }

  if(icas & 2) {
    //printf("icas %d --> 1 ?\n",icas);
    pt = &mesh->tetra[k + 5];
    pt->v[0] = p0;
    pt->v[1] = p5;
    pt->v[2] = p3;
    pt->v[3] = mesh->np;
    pt->ref  = ref;
    pt = &mesh->tetra[k + 6];
    pt->v[0] = p0;
    pt->v[1] = p5;
    pt->v[2] = p2;
    pt->v[3] = mesh->np;
    pt->ref  = ref;
  } else {
    if(icas & 16) {
      //printf("icas %d --> 4 ?\n",icas);

    } else {
      MMG_edgePut(hed,p2,p3,2);
    }
    pt = &mesh->tetra[k + 5];
    pt->v[0] = p0;
    pt->v[1] = p2;
    pt->v[2] = p3;
    pt->v[3] = mesh->np;
    pt->ref  = ref;
    pt = &mesh->tetra[k + 6];
    pt->v[0] = p2;
    pt->v[1] = p3;
    pt->v[2] = p5;
    pt->v[3] = mesh->np;
    pt->ref  = ref;
  }

  pt = &mesh->tetra[k + 7];
  pt->v[0] = p3;
  pt->v[1] = p4;
  pt->v[2] = p5;
  pt->v[3] = mesh->np;
  pt->ref  = ref;
  pt = &mesh->tetra[k + 8];
  pt->v[0] = p0;
  pt->v[1] = p1;
  pt->v[2] = p2;
  pt->v[3] = mesh->np;
  pt->ref  = ref;

  return(1);
}

//cut the hex in 6 tetras
// int MMG_cuthex(pMesh mesh,pHedge hed,int k,int p0,int p1,int p2,int p3,int p4,int p5,int p6,int p7, int ref) {
//      pTetra          pt;
//      int                             i,nu1,nu2;
//      int       p[8],mini,minnumber,nhex;
//      double    volhex;
//
//      /*permutation to have compatibility*/
//      /*1) the first face have to contain the smallest vertex number*/
//      p[0] = p0; p[1] = p1; p[2] = p2 ; p[3] = p3; p[4] = p4; p[5] = p5; p[6] = p6; p[7] = p7;
//      printf("points : %d %d %d %d %d %d %d %d \n",p0,p1,p2,p3,p4,p5,p6,p7);
//      minnumber = p[0];
//      mini = 0;
//      for(i=1 ; i<8 ; i++) {
//              if(p[i] < minnumber) {
//                      minnumber = p[i];
//                      mini = i;
//              }
//      }
//      if(mini > 3) {
//              printf("on doit echanger la face dessous et la face dessus\n");
//              exit(0);
//      } else {
//              if(mini) printf("mini %d --> %d\n",mini,k/6+1);
//              nhex = k/6+1;
//              //if(nhex==84841 || nhex==84842) {
//                      volhex = MMG_quickvol(mesh->point[p0].c,mesh->point[p1].c,mesh->point[p3].c,mesh->point[p4].c);
//                      printf("orientation of %d : %e\n",nhex,volhex);
//                      printf("points : %d %d %d %d %d %d %d %d -- %d\n",p0,p1,p2,p3,p4,p5,p6,p7,mini);
//              //}
//              switch(mini){
//                case(0):
//                        p[0] = p3; p[4] = p7;
//                        p[1] = p0; p[5] = p4;
//                        p[2] = p1; p[6] = p5;
//                        p[3] = p2; p[7] = p6;
//                break;
//                case(2):
//                      exit(0);
//                      case(3):
//                        p[0] = p2; p[4] = p6;
//                        p[1] = p3; p[5] = p7;
//                        p[2] = p0; p[6] = p4;
//                        p[3] = p1; p[7] = p5;
//                      break;
//                      default:
//                      break;
//              }
//      }
//
//      pt = &mesh->tetra[k+1];
//   pt->v[0] = p[0];
//   pt->v[1] = p[1];
//   pt->v[2] = p[3];
//   pt->v[3] = p[7];
//   pt->ref  = ref;
//      for(i=0 ; i<6 ; i++) {
//              nu1 = pt->v[MMG_iare[i][0]];
//              nu2 = pt->v[MMG_iare[i][1]];
//              MMG_edgePut(hed,nu1,nu2,2);
//      }
//   //if((netmp+(k-1)*6+1 ) == 2924) printf("i) %d tet %d %d %d %d\n",k,pt->v[0],pt->v[1],pt->v[2],pt->v[3]);
//   pt = &mesh->tetra[k+2];
//   pt->v[0] = p[7];
//   pt->v[1] = p[2];
//   pt->v[2] = p[1];
//   pt->v[3] = p[6];
//   pt->ref  = ref;
//      for(i=0 ; i<6 ; i++) {
//              nu1 = pt->v[MMG_iare[i][0]];
//              nu2 = pt->v[MMG_iare[i][1]];
//              MMG_edgePut(hed,nu1,nu2,2);
//      }
//   //if((netmp+(k-1)*6+1 + 1) == 2924) printf("ii) tet %d %d %d %d\n",pt->v[0],pt->v[1],pt->v[2],pt->v[3]);
//   pt = &mesh->tetra[k+3];
//   pt->v[0] = p[1];
//   pt->v[1] = p[4];
//   pt->v[2] = p[5];
//   pt->v[3] = p[7];
//      for(i=0 ; i<6 ; i++) {
//              nu1 = pt->v[MMG_iare[i][0]];
//              nu2 = pt->v[MMG_iare[i][1]];
//              MMG_edgePut(hed,nu1,nu2,2);
//      }
//   //if((netmp+(k-1)*6+1 + 2) == 2924) printf("iii) tet %d %d %d %d\n",pt->v[0],pt->v[1],pt->v[2],pt->v[3]);
//   pt->ref  = ref;
//   pt = &mesh->tetra[k+4];
//   pt->v[0] = p[7];
//   pt->v[1] = p[4];
//   pt->v[2] = p[0];
//   pt->v[3] = p[1];
//   pt->ref  = ref;
//      for(i=0 ; i<6 ; i++) {
//              nu1 = pt->v[MMG_iare[i][0]];
//              nu2 = pt->v[MMG_iare[i][1]];
//              MMG_edgePut(hed,nu1,nu2,2);
//      }
//   //if((netmp+(k-1)*6+1 + 3) == 2924) printf("iv) tet %d %d %d %d\n",pt->v[0],pt->v[1],pt->v[2],pt->v[3]);
//   pt = &mesh->tetra[k+5];
//   pt->v[0] = p[1];
//   pt->v[1] = p[6];
//   pt->v[2] = p[7];
//   pt->v[3] = p[5];
//   pt->ref  = ref;
//      for(i=0 ; i<6 ; i++) {
//              nu1 = pt->v[MMG_iare[i][0]];
//              nu2 = pt->v[MMG_iare[i][1]];
//              MMG_edgePut(hed,nu1,nu2,2);
//      }
//   pt = &mesh->tetra[k+6];
//   pt->v[0] = p[1];
//   pt->v[1] = p[3];
//   pt->v[2] = p[2];
//   pt->v[3] = p[7];
//   pt->ref  = ref;
//      for(i=0 ; i<6 ; i++) {
//              nu1 = pt->v[MMG_iare[i][0]];
//              nu2 = pt->v[MMG_iare[i][1]];
//              MMG_edgePut(hed,nu1,nu2,2);
//      }
//      return(1);
// }
int ddebug=0;

int MMG_decouphex(pMesh mesh, pHedge hed,int k,int* p,int ref) {
  pTetra  pt;
  int     i,nu1,nu2;

  pt = &mesh->tetra[k+1];
  pt->v[0] = p[0];
  pt->v[1] = p[1];
  pt->v[2] = p[3];
  pt->v[3] = p[7];
  pt->ref  = ref;
  // for(i=0 ; i<6 ; i++) {
  //    nu1 = pt->v[MMG_iare[i][0]];
  //    nu2 = pt->v[MMG_iare[i][1]];
  //    MMG_edgePut(hed,nu1,nu2,2);
  //    if(ddebug) printf("1) on rajoute %d %d\n",nu1,nu2);
  // }
  //if((netmp+(k-1)*6+1 ) == 2924) printf("i) %d tet %d %d %d %d\n",k,pt->v[0],pt->v[1],pt->v[2],pt->v[3]);
  pt = &mesh->tetra[k+2];
  pt->v[0] = p[7];
  pt->v[1] = p[2];
  pt->v[2] = p[1];
  pt->v[3] = p[6];
  pt->ref  = ref;
  // for(i=0 ; i<6 ; i++) {
  //    nu1 = pt->v[MMG_iare[i][0]];
  //    nu2 = pt->v[MMG_iare[i][1]];
  //    MMG_edgePut(hed,nu1,nu2,2);
  //    if(ddebug) printf("2) on rajoute %d %d\n",nu1,nu2);
  // }
  //if((netmp+(k-1)*6+1 + 1) == 2924) printf("ii) tet %d %d %d %d\n",pt->v[0],pt->v[1],pt->v[2],pt->v[3]);
  pt = &mesh->tetra[k+3];
  pt->v[0] = p[1];
  pt->v[1] = p[4];
  pt->v[2] = p[5];
  pt->v[3] = p[7];
  // for(i=0 ; i<6 ; i++) {
  //    nu1 = pt->v[MMG_iare[i][0]];
  //    nu2 = pt->v[MMG_iare[i][1]];
  //    MMG_edgePut(hed,nu1,nu2,2);
  //    if(ddebug) printf("3) on rajoute %d %d\n",nu1,nu2);
  // }
  //if((netmp+(k-1)*6+1 + 2) == 2924) printf("iii) tet %d %d %d %d\n",pt->v[0],pt->v[1],pt->v[2],pt->v[3]);
  pt->ref  = ref;
  pt = &mesh->tetra[k+4];
  pt->v[0] = p[7];
  pt->v[1] = p[4];
  pt->v[2] = p[0];
  pt->v[3] = p[1];
  pt->ref  = ref;
  // for(i=0 ; i<6 ; i++) {
  //    nu1 = pt->v[MMG_iare[i][0]];
  //    nu2 = pt->v[MMG_iare[i][1]];
  //    MMG_edgePut(hed,nu1,nu2,2);
  //    if(ddebug) printf("4) on rajoute %d %d\n",nu1,nu2);
  // }
  //if((netmp+(k-1)*6+1 + 3) == 2924) printf("iv) tet %d %d %d %d\n",pt->v[0],pt->v[1],pt->v[2],pt->v[3]);
  pt = &mesh->tetra[k+5];
  pt->v[0] = p[1];
  pt->v[1] = p[6];
  pt->v[2] = p[7];
  pt->v[3] = p[5];
  pt->ref  = ref;
  // for(i=0 ; i<6 ; i++) {
  //    nu1 = pt->v[MMG_iare[i][0]];
  //    nu2 = pt->v[MMG_iare[i][1]];
  //    MMG_edgePut(hed,nu1,nu2,2);
  //    if(ddebug) printf("5) on rajoute %d %d\n",nu1,nu2);
  // }
  pt = &mesh->tetra[k+6];
  pt->v[0] = p[1];
  pt->v[1] = p[3];
  pt->v[2] = p[2];
  pt->v[3] = p[7];
  pt->ref  = ref;
  // for(i=0 ; i<6 ; i++) {
  //    nu1 = pt->v[MMG_iare[i][0]];
  //    nu2 = pt->v[MMG_iare[i][1]];
  //    MMG_edgePut(hed,nu1,nu2,2);
  //    if(ddebug) printf("6) on rajoute %d %d\n",nu1,nu2);
  // }

  //add edges at hashtable
  MMG_edgePut(hed,p[0],p[7],2);
  MMG_edgePut(hed,p[1],p[3],2);
  MMG_edgePut(hed,p[2],p[7],2);
  MMG_edgePut(hed,p[1],p[6],2);
  MMG_edgePut(hed,p[1],p[4],2);
  MMG_edgePut(hed,p[5],p[7],2);
  return(1);
}

int MMG_checkcaseopp(pHexa ph,int nu1,int nu2,pHedge hed) {
  int i,nu3,nu4;
  unsigned char MMG_hied[8][3] = { {2,5,7}, {3,4,6}, {0,5,7}, {1,4,6}, {1,3,6}, {0,2,7}, {1,3,4}, {0,2,5} };
  unsigned char MMG_hop[8] = { 6,7,4,5,2,3,0,1 };
  nu4 = MMG_hop[nu1];
  for(i=0;i<3;i++) {
    nu3 = MMG_hied[nu4][i];
    if(nu3==MMG_hop[nu2]) continue;
    if(MMG_edgePoint(hed,ph->v[nu4],ph->v[nu3])) break;
  }
  if(i<3) return(1);
  else return(0);
}

int MMG_checkcase(pHexa ph,int nu1,int nu2,pHedge hed) {
  int i,nu3;
  unsigned char MMG_hied[8][3] = { {2,5,7}, {3,4,6}, {0,5,7}, {1,4,6}, {1,3,6}, {0,2,7}, {1,3,4}, {0,2,5} };

  for(i=0;i<3;i++) {
    nu3 = MMG_hied[nu1][i];
    if(nu3==nu2) continue;
    if(MMG_edgePoint(hed,ph->v[nu1],ph->v[nu3])) break;
  }
  if(i<3) return(1);
  else return(0);

}
int MMG_cuthex(pMesh mesh,pHedge hed,pHexa listhexa,int* adjahex,int nhex,int netmp) {
  pHexa     ph;
  pTetra    pt;
  pPoint    ppt;
  int                           i,k,nu1,nu2,nu3,nu4,adj,icas0,icasopp,nncut;
  int       *list,p[8],mini,minnumber,ipil,icurc,iface,iadr;
  int       iel,ip;
  unsigned char MMG_hidir[6][4] = { {0,3,2,1}, {0,4,7,3}, {0,1,5,4}, {4,5,6,7}, {1,2,6,5}, {2,3,7,6} };
  unsigned char MMG_hidirop[6][4] = { {7,4,5,6}, {5,1,2,6}, {7,3,2,6}, {1,0,3,2}, {3,0,4,7}, {0,1,5,4} };
  unsigned char MMG_hied[8][3] = { {2,5,7}, {3,4,6}, {0,5,7}, {1,4,6}, {1,3,6}, {0,2,7}, {1,3,4}, {0,2,5} };
  unsigned char MMG_hop[8] = { 6,7,4,5,2,3,0,1 };
  double    volhex,c[3];

  /*alloc*/
	list = (int*)M_calloc(10*nhex,sizeof(int),"alloclist");
  assert(adjahex);                             

  /*init pile*/
  ph = &listhexa[1];
  ph->mark = -1;
  for(i=0;i<8;i++) p[i] = ph->v[i];
  MMG_decouphex(mesh,hed,netmp,p,ph->ref);
  icurc=0;ipil=0;
  for(i=0;i<6;i++) {
    iadr = 1;
    adj = adjahex[iadr + i];
    if(!adj) continue;
    list[ipil++] = adj;
  }
  while(icurc++ < ipil) {
    k = list[icurc-1]/6;
    if(!k) continue;
    ph = &listhexa[k];
    if(ph->mark < 0) continue;
    ph->mark = -1;
    iface = list[icurc-1]%6;
    ddebug=0;
    if(ddebug) {
      printf("** eh hop on traite %d iface %d\n",k,iface);
      printf("face %d %d %d %d\n",ph->v[MMG_hidir[iface][0]],ph->v[MMG_hidir[iface][1]],
	     ph->v[MMG_hidir[iface][2]],ph->v[MMG_hidir[iface][3]]);
    }
    nu1 = MMG_hidir[iface][0];
    nu2 = MMG_hidir[iface][2];
    if(MMG_edgePoint(hed,ph->v[nu1],ph->v[nu2])) {
      if(ddebug)  printf("on a trouve l'arete %d %d\n",ph->v[MMG_hidir[iface][0]],ph->v[MMG_hidir[iface][2]]);
      if(ddebug) printf("iface %d : %d %d %d %d %d %d %d %d\n",iface,ph->v[0],ph->v[1],ph->v[2],ph->v[3],ph->v[4],ph->v[5],ph->v[6],ph->v[7]);
      //if edge opp sur face opp exist, on degage tout de suite
      nu3 = MMG_hidirop[iface][0];
      nu4 = MMG_hidirop[iface][2];
      if(MMG_edgePoint(hed,ph->v[nu3],ph->v[nu4])) {
	ph->mark = -10;
	continue;
      }
      if(iface==1 || iface==5) {
	//find if other edge with ph->v[MMG_hidir[iface][0]], if yes->renum
	icas0 = MMG_checkcase(ph,nu1,nu2,hed);
	icasopp = MMG_checkcaseopp(ph,nu1,nu2,hed);
	if(!icas0 && !icasopp) {
	  icas0 = 0;
	} else {
	  icas0 = MMG_checkcase(ph,nu2,nu1,hed);
	  icasopp = MMG_checkcaseopp(ph,nu2,nu1,hed);
	  if(icas0 || icasopp) {
	    ph->mark = -10;
	    continue;
	  }
	  icas0 = 1;
	}
	if(icas0) {
	  //debug check
	  for(i=0;i<3;i++) {
	    nu3 = MMG_hied[nu2][i];
	    if(nu3==nu1) continue;
	    if(MMG_edgePoint(hed,ph->v[nu2],ph->v[nu3])) break;
	  }
	  assert(i==3);
	  //printf("iface %d on a trouve une autre arete---> renum\n",iface);
	  if(iface==1) {
	    p[0] = ph->v[3]; p[1] = ph->v[0]; p[2] = ph->v[1]; p[3] = ph->v[2];
	    p[4] = ph->v[7]; p[5] = ph->v[4]; p[6] = ph->v[5]; p[7] = ph->v[6];
	  } else {
	    p[0] = ph->v[1]; p[1] = ph->v[2]; p[2] = ph->v[3]; p[3] = ph->v[0];
	    p[4] = ph->v[5]; p[5] = ph->v[6]; p[6] = ph->v[7]; p[7] = ph->v[4];
	  }
	  MMG_decouphex(mesh,hed,netmp+6*(k-1),p,ph->ref);
	} else {
	  MMG_decouphex(mesh,hed,netmp+6*(k-1),ph->v,ph->ref);
	}
      } else if (iface==4) {
	icas0 = MMG_checkcase(ph,nu2,nu1,hed);
	icasopp = MMG_checkcaseopp(ph,nu2,nu1,hed);
	if(!icas0 && !icasopp) {
	  icas0 = 0;
	} else {
	  icas0 = MMG_checkcase(ph,nu1,nu2,hed);
	  icasopp = MMG_checkcaseopp(ph,nu1,nu2,hed);
	  if(icas0 || icasopp) {
	    ph->mark = -10;
	    continue;
	  }
	  icas0 = 1;
	}
	if(icas0) {
	  //check debug
	  for(i=0;i<3;i++) {
	    nu3 = MMG_hied[nu1][i];
	    if(nu3==nu2) continue;
	    if(MMG_edgePoint(hed,ph->v[nu1],ph->v[nu3])) break;
	  }
	  assert(i==3);
	  //printf("iface 4 on a trouve une autre arete---> renum\n");
	  p[0] = ph->v[3]; p[1] = ph->v[0]; p[2] = ph->v[1]; p[3] = ph->v[2];
	  p[4] = ph->v[7]; p[5] = ph->v[4]; p[6] = ph->v[5]; p[7] = ph->v[6];
	  MMG_decouphex(mesh,hed,netmp+6*(k-1),p,ph->ref);
	} else {
	  MMG_decouphex(mesh,hed,netmp+6*(k-1),ph->v,ph->ref);
	}
      } else {
	if(ddebug)  printf("il faut renum iface %d\n",iface);//iface 0,2,3
	icas0 = MMG_checkcase(ph,nu1,nu2,hed);
	icasopp = MMG_checkcaseopp(ph,nu1,nu2,hed);
	if(!icas0 && !icasopp) {
	  icas0 = 0;
	} else {
	  icas0 = MMG_checkcase(ph,nu2,nu1,hed);
	  icasopp = MMG_checkcaseopp(ph,nu2,nu1,hed);
	  if(icas0 || icasopp) {
	    ph->mark = -10;
	    continue;
	  }
	  icas0 = 1;
	}
	if(icas0) {
	  //check debug
	  for(i=0;i<3;i++) {
	    nu3 = MMG_hied[nu2][i];
	    if(nu3==nu1) continue;
	    if(MMG_edgePoint(hed,ph->v[nu2],ph->v[nu3])) break;
	  }
	  assert(i==3);
	  icas0=1;
	}
	if(ddebug) printf("icas %d\n",icas0);
	switch(iface) {
	case(0):
	  if(icas0) {
	    p[0] = ph->v[3]; p[1] = ph->v[0]; p[2] = ph->v[1]; p[3] = ph->v[2];
	    p[4] = ph->v[7]; p[5] = ph->v[4]; p[6] = ph->v[5]; p[7] = ph->v[6];
	  } else {
	    p[0] = ph->v[1]; p[1] = ph->v[2]; p[2] = ph->v[3]; p[3] = ph->v[0];
	    p[4] = ph->v[5]; p[5] = ph->v[6]; p[6] = ph->v[7]; p[7] = ph->v[4];
	  }
	  break;
	case(2):
	  if(icas0) {
	    p[0] = ph->v[3]; p[1] = ph->v[0]; p[2] = ph->v[1]; p[3] = ph->v[2];
	    p[4] = ph->v[7]; p[5] = ph->v[4]; p[6] = ph->v[5]; p[7] = ph->v[6];
	  } else {
	    p[0] = ph->v[2]; p[1] = ph->v[3]; p[2] = ph->v[0]; p[3] = ph->v[1];
	    p[4] = ph->v[6]; p[5] = ph->v[7]; p[6] = ph->v[4]; p[7] = ph->v[5];
	  }
	  break;
	case(3):
	  if(icas0){// ph->v[MMG_hidir[iface][0]]<ph->v[MMG_hidir[iface][2]]) {
	    p[0] = ph->v[1]; p[1] = ph->v[2]; p[2] = ph->v[3]; p[3] = ph->v[0];
	    p[4] = ph->v[5]; p[5] = ph->v[6]; p[6] = ph->v[7]; p[7] = ph->v[4];
	  } else {
	    p[0] = ph->v[3]; p[1] = ph->v[0]; p[2] = ph->v[1]; p[3] = ph->v[2];
	    p[4] = ph->v[7]; p[5] = ph->v[4]; p[6] = ph->v[5]; p[7] = ph->v[6];
	  }
	  break;
	}
	if(ddebug) printf("after renum : %d %d %d %d %d %d %d %d\n",p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]);
	MMG_decouphex(mesh,hed,netmp+6*(k-1),p,ph->ref);
      }
    }  else if (MMG_edgePoint(hed,ph->v[MMG_hidir[iface][1]],ph->v[MMG_hidir[iface][3]])) {
      nu1 = MMG_hidir[iface][1];
      nu2 = MMG_hidir[iface][3];
      if(ddebug)  printf("-- on a trouve l'arete %d %d\n",ph->v[MMG_hidir[iface][1]],ph->v[MMG_hidir[iface][3]]);
      if(ddebug) printf("iface %d : %d %d %d %d %d %d %d %d\n",iface,ph->v[0],ph->v[1],ph->v[2],ph->v[3],ph->v[4],ph->v[5],ph->v[6],ph->v[7]);
      //if edge opp face opp exist on degage
      nu3 = MMG_hidirop[iface][1];
      nu4 = MMG_hidirop[iface][3];
      if(MMG_edgePoint(hed,ph->v[nu3],ph->v[nu4])) {
	ph->mark = -10;
	continue;
      }
      if(iface==0 || iface==3) {
	icas0 = MMG_checkcase(ph,nu1,nu2,hed);
	icasopp = MMG_checkcaseopp(ph,nu1,nu2,hed);
	if(ddebug) printf("oppnu1 %d  : %d %d\n",ph->v[MMG_hop[nu1]],icas0,icasopp);
	if(!icas0 && !icasopp) {
	  icas0 = 0;
	} else {
	  icas0 = MMG_checkcase(ph,nu2,nu1,hed);
	  icasopp = MMG_checkcaseopp(ph,nu2,nu1,hed);
	  if(icas0 || icasopp) {
	    ph->mark = -10;
	    continue;
	  }
	  icas0 = 1;
	}
	if(ddebug) printf("icas0 %d\n",icas0);

	if(icas0)       {
	  //check debug
	  for(i=0;i<3;i++) {
	    nu3 = MMG_hied[nu2][i];
	    if(nu3==nu1) continue;
	    if(MMG_edgePoint(hed,ph->v[nu2],ph->v[nu3])) {printf("on trouve arg %d %d\n",ph->v[nu2],ph->v[nu3]);break;    }
	  }
	  assert(i==3);
	  if(iface==0) {
	    p[0] = ph->v[2]; p[1] = ph->v[3]; p[2] = ph->v[0]; p[3] = ph->v[1];
	    p[4] = ph->v[6]; p[5] = ph->v[7]; p[6] = ph->v[4]; p[7] = ph->v[5];
	  } else {
	    p[0] = ph->v[2]; p[1] = ph->v[3]; p[2] = ph->v[0]; p[3] = ph->v[1];
	    p[4] = ph->v[6]; p[5] = ph->v[7]; p[6] = ph->v[4]; p[7] = ph->v[5];
	  }
	  MMG_decouphex(mesh,hed,netmp+6*(k-1),p,ph->ref);
	} else {
	  MMG_decouphex(mesh,hed,netmp+6*(k-1),ph->v,ph->ref);
	}
      } else if(iface==2){
	icas0 = MMG_checkcase(ph,nu2,nu1,hed);
	icasopp = MMG_checkcaseopp(ph,nu2,nu1,hed);
	if(!icas0 && !icasopp) {
	  icas0 = 0;
	} else {
	  icas0 = MMG_checkcase(ph,nu1,nu2,hed);
	  icasopp = MMG_checkcaseopp(ph,nu1,nu2,hed);
	  if(icas0 || icasopp) {
	    ph->mark = -10;
	    continue;
	  }
	  icas0 = 1;
	}
	if(icas0) {
	  //check debug
	  for(i=0;i<3;i++) {
	    nu3 = MMG_hied[nu1][i];
	    if(nu3==nu2) continue;
	    if(MMG_edgePoint(hed,ph->v[nu1],ph->v[nu3])) break;
	  }
	  assert(i==3);
	  p[0] = ph->v[1]; p[1] = ph->v[2]; p[2] = ph->v[3]; p[3] = ph->v[0];
	  p[4] = ph->v[5]; p[5] = ph->v[6]; p[6] = ph->v[7]; p[7] = ph->v[4];
	  MMG_decouphex(mesh,hed,netmp+6*(k-1),p,ph->ref);
	} else {
	  MMG_decouphex(mesh,hed,netmp+6*(k-1),ph->v,ph->ref);
	}
      }
      else {
	if(ddebug)  printf("il faut renum iface %d\n",iface);//iface 1,4,5
	icas0 = MMG_checkcase(ph,nu1,nu2,hed);
	icasopp = MMG_checkcaseopp(ph,nu1,nu2,hed);
	if(!icas0 && !icasopp) {
	  icas0 = 0;
	} else {
	  icas0 = MMG_checkcase(ph,nu2,nu1,hed);
	  icasopp = MMG_checkcaseopp(ph,nu2,nu1,hed);
	  if(icas0 || icasopp) {
	    ph->mark = -10;
	    continue;
	  }
	  icas0 = 1;
	}
	switch(iface) {
	case(1):
	  if(icas0){//ph->v[MMG_hidir[iface][1]]<ph->v[MMG_hidir[iface][3]]) {
	    p[0] = ph->v[1]; p[1] = ph->v[2]; p[2] = ph->v[3]; p[3] = ph->v[0];
	    p[4] = ph->v[5]; p[5] = ph->v[6]; p[6] = ph->v[7]; p[7] = ph->v[4];
	  } else {
	    p[0] = ph->v[2]; p[1] = ph->v[3]; p[2] = ph->v[0]; p[3] = ph->v[1];
	    p[4] = ph->v[6]; p[5] = ph->v[7]; p[6] = ph->v[4]; p[7] = ph->v[5];
	  }
	  break;
	case(4):
	  if(ddebug) printf("at the beginning %d : %d %d %d %d %d %d %d %d\n",k,ph->v[0],ph->v[1],ph->v[2],ph->v[3],ph->v[4],ph->v[5],ph->v[6],ph->v[7]);
	  if(icas0) {
	    p[0] = ph->v[1]; p[1] = ph->v[2]; p[2] = ph->v[3]; p[3] = ph->v[0];
	    p[4] = ph->v[5]; p[5] = ph->v[6]; p[6] = ph->v[7]; p[7] = ph->v[4];
	  } else {
	    p[0] = ph->v[2]; p[1] = ph->v[3]; p[2] = ph->v[0]; p[3] = ph->v[1];
	    p[4] = ph->v[6]; p[5] = ph->v[7]; p[6] = ph->v[4]; p[7] = ph->v[5];
	  }
	  if(ddebug)  printf("at the end %d : %d %d %d %d %d %d %d %d\n",k,p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]);
	  break;
	case(5):
	  if(icas0) {
	    p[0] = ph->v[2]; p[1] = ph->v[3]; p[2] = ph->v[0]; p[3] = ph->v[1];
	    p[4] = ph->v[6]; p[5] = ph->v[7]; p[6] = ph->v[4]; p[7] = ph->v[5];
	  } else {
	    p[0] = ph->v[3]; p[1] = ph->v[0]; p[2] = ph->v[1]; p[3] = ph->v[2];
	    p[4] = ph->v[7]; p[5] = ph->v[4]; p[6] = ph->v[5]; p[7] = ph->v[6];
	  }
	  break;
	}
	MMG_decouphex(mesh,hed,netmp+6*(k-1),p,ph->ref);
      }
    } else {
      printf("ya un pbs!!!!! %d %d et %d %d\n",ph->v[MMG_hidir[iface][1]],ph->v[MMG_hidir[iface][3]],
	     ph->v[MMG_hidir[iface][0]],ph->v[MMG_hidir[iface][2]]);
      exit(0);
    }
    for(i=0;i<6;i++) {
      iadr = 6*(k-1)+1;
      adj = adjahex[iadr + i];
      if(!adj) continue;
      list[ipil++] = adj;
      //printf("k= %d on rajoute hexa %d (iface %d) en %d -- par face %d \n",k,adj/6,adj%6,ipil-1,i);
    }
    //printf("---- end ipil %d icurc %d\n",ipil,icurc);
  }

  //stay few hexa not cutted...
  nncut = 0;
  for(k=1 ; k<=nhex ; k++) {
    ph = &listhexa[k];
    if(ph->mark==-1) continue;
    nncut++;
    if(k==60478) ddebug = 1;
    //create new vertex
    c[0] = c[1] = c[2] = 0.;
    for(i=0 ; i<8 ; i++) {
      ppt = &mesh->point[ph->v[i]];
      c[0] += ppt->c[0];        c[1] += ppt->c[1]; c[2] += ppt->c[2];
    }
    c[0] /= 8.; c[1] /= 8.; c[2] /= 8.;
    ip   = MMG_newPt(mesh,c);
    if ( ip < 1 )  {
      fprintf(stdout,"ERROR VERTEX CREATION\n");
      return(0);
    }
    //create 2 tets per faces
    for(i=0 ; i<6 ; i++) {
      nu1 = MMG_hidir[i][0];
      nu2 = MMG_hidir[i][2];
      if(MMG_edgePoint(hed,ph->v[nu1],ph->v[nu2])) {
	iel = MMG_newElt(mesh);
	pt = &mesh->tetra[iel];
	pt->v[0] = ip;
	pt->v[1] = ph->v[nu1];
	pt->v[2] = ph->v[nu2];
	pt->v[3] = ph->v[MMG_hidir[i][1]];
	pt->ref  = ph->ref;
	if(MMG_voltet(mesh,iel)<0) {
	  pt->v[3] = ph->v[nu2];
	  pt->v[2] = ph->v[MMG_hidir[i][1]];
	}
	iel = MMG_newElt(mesh);
	pt = &mesh->tetra[iel];
	pt->v[0] = ip;
	pt->v[1] = ph->v[nu1];
	pt->v[2] = ph->v[MMG_hidir[i][3]];
	pt->v[3] = ph->v[nu2];
	pt->ref  = ph->ref;
	if(MMG_voltet(mesh,iel)<0) {
	  pt->v[3] = ph->v[MMG_hidir[i][3]];
	  pt->v[2] = ph->v[nu2];
	}
      } else {
	nu1 = MMG_hidir[i][1];
	nu2 = MMG_hidir[i][3];
	if(!MMG_edgePoint(hed,ph->v[nu1],ph->v[nu2])) MMG_edgePut(hed,ph->v[nu1],ph->v[nu2],2);
	iel = MMG_newElt(mesh);
	pt = &mesh->tetra[iel];
	pt->v[0] = ip;
	pt->v[1] = ph->v[nu1];
	pt->v[2] = ph->v[MMG_hidir[i][0]];
	pt->v[3] = ph->v[nu2];
	pt->ref  = ph->ref;
	if(MMG_voltet(mesh,iel)<0) {
	  pt->v[3] = ph->v[MMG_hidir[i][0]];
	  pt->v[2] = ph->v[nu2];
	}
	iel = MMG_newElt(mesh);
	pt = &mesh->tetra[iel];
	pt->v[0] = ip;
	pt->v[1] = ph->v[nu1];
	pt->v[2] = ph->v[nu2];
	pt->v[3] = ph->v[MMG_hidir[i][2]];
	pt->ref  = ph->ref;
	if(MMG_voltet(mesh,iel)<0) {
	  pt->v[3] = ph->v[MMG_hidir[i][2]];
	  pt->v[2] = ph->v[nu2];
	}
      }
    }
  }
  if(nncut) fprintf(stdout,"  $$ %8d ADDED VERTEX\n",nncut);
  return(1);
}

int MMG_cutprism(pMesh mesh,pHedge hed,int k,int p0,int p1,int p2,int p3,int p4,int p5,int ref) {
  pTetra        pt;
  double        vol;
  int                   t0,t1,t2;
  char          icas;
  int ddebug;

  vol= MMG_quickvol(mesh->point[p0].c,mesh->point[p1].c,mesh->point[p2].c,mesh->point[p3].c);
  if(vol<0) {
    printf("inversion");
    t0 = p0;
    t1 = p1;
    t2 = p2;
    p0 = p3;
    p1 = p4;
    p2 = p5;
    p3 = t0;
    p4 = t1;
    p5 = t2;
  }
  if(k==606 || k==605 || k==604 || k==609 || k==608 || k==607) ddebug=1;
  else ddebug=0;
  ddebug=0;
  if(ddebug) printf("k = %d : %d %d %d %d %d %d\n",k,p0,p1,p2,p3,p4,p5);

  icas = 0;

  //find edge 2 : p1-p3 then edge 0 : 0 4
  if(!MMG_edgePoint(hed,p1,p3)) {
    if(MMG_edgePoint(hed,p0,p4))
      icas |= 1;
  } else {
    icas |= 4;
  }
  //find edge 5 : p2-p4 then edge 3 : 1 5
  if(!MMG_edgePoint(hed,p2,p4)) {
    if(MMG_edgePoint(hed,p1,p5))
      icas |= 8;
  } else {
    icas |= 32;
  }
  //find edge 4 : p2-p3 then edge 1 : 0 5
  if(!MMG_edgePoint(hed,p2,p3)) {
    if(MMG_edgePoint(hed,p0,p5))
      icas |= 2;
  } else {
    icas |= 16;
  }
  if(icas > 55) {
    fprintf(stdout,"grosgros bug %d\n",icas);
    exit(0);
  }
  if(ddebug) printf("on trouve %d\n",icas);

  switch(icas) {
  case 0:
    if(ddebug) printf("on rajoute %d %d -- %d %d -- %d %d\n",p0,p4,p1,p5,p0,p5);
    MMG_edgePut(hed,p2,p4,2);
    MMG_edgePut(hed,p1,p3,2);
    MMG_edgePut(hed,p3,p2,2);//MMG_edgePut(hed,p2,p3,2);
    icas = 52;
    break;
  case 1:
    MMG_edgePut(hed,p1,p5,2);
    MMG_edgePut(hed,p0,p5,2);
    icas = 11;//25;
    break;
  case 2:
    MMG_edgePut(hed,p1,p5,2);
    MMG_edgePut(hed,p0,p4,2);
    icas = 11;//14
    break;
  case 3:
    MMG_edgePut(hed,p1,p5,2);
    icas = 11;//35;
    break;
  case 4:
    MMG_edgePut(hed,p2,p4,2);//MMG_edgePut(hed,p1,p5,2);
    MMG_edgePut(hed,p2,p3,2);//MMG_edgePut(hed,p0,p5,2);
    icas = 52;//14;
    break;
  case 6:
    MMG_edgePut(hed,p1,p5,2);
    icas = 14;
    break;
  case 8:
    MMG_edgePut(hed,p0,p5,2);
    MMG_edgePut(hed,p0,p4,2);
    icas = 11;//14;
    break;
  case 9:
    MMG_edgePut(hed,p0,p5,2);
    icas = 11;//25;
    break;
  case 10:
    MMG_edgePut(hed,p0,p4,2);
    icas = 11;//14;
    break;
  case 12:
    MMG_edgePut(hed,p0,p5,2);
    icas = 14;
    break;
  case 16:
    MMG_edgePut(hed,p2,p4,2);//MMG_edgePut(hed,p1,p5,2);
    MMG_edgePut(hed,p3,p1,2);//MMG_edgePut(hed,p1,p3,2);
    icas = 52;//28;
    break;
  case 17:
    MMG_edgePut(hed,p4,p2,2);
    icas = 49;//25;
    break;
  case 20:
    MMG_edgePut(hed,p2,p4,2);    //MMG_edgePut(hed,p1,p5,2);
    icas = 52;//28;
    break;
  case 24:
    MMG_edgePut(hed,p1,p3,2);
    icas = 28;//25;
    break;
  case 32:
    MMG_edgePut(hed,p1,p3,2);//MMG_edgePut(hed,p0,p4,2);
    MMG_edgePut(hed,p2,p3,2);//MMG_edgePut(hed,p0,p5,2);
    icas = 52;//35;
    break;
  case 33:
    MMG_edgePut(hed,p0,p5,2);
    icas = 35;
    break;
  case 34:
    MMG_edgePut(hed,p0,p4,2);
    icas = 35;
    break;
  case 36:
    MMG_edgePut(hed,p3,p2,2);
    icas = 52;
    break;
  case 48:
    MMG_edgePut(hed,p1,p3,2);//MMG_edgePut(hed,p0,p4,2);
    icas = 52;//49;
    break;
  default:
    //5,7,11,13,15,18,19,21,22,23,26,27,29,30,31,37,39,40,41,42,43,44,45,46,47,50,51,52,53,54,55:
    //printf("icas imposssss %d\n",icas);
    //exit(0);
    break;

  }
  if(ddebug) printf("du coup %d\n",icas);
  switch(icas) {
  case 14:
    pt = &mesh->tetra[k + 1];
    pt->v[0] = p5;
    pt->v[1] = p1;
    pt->v[2] = p2;
    pt->v[3] = p0;
    pt->ref  = ref;//1;//ref;
    pt = &mesh->tetra[k + 2];
    pt->v[0] = p3;
    pt->v[1] = p5;
    pt->v[2] = p1;
    pt->v[3] = p0;
    pt->ref  = ref;//1;//ref;
    pt = &mesh->tetra[k + 3];
    pt->v[0] = p3;
    pt->v[1] = p4;
    pt->v[2] = p1;
    pt->v[3] = p5;
    pt->ref  = ref;//1;//ref;
    break;
  case 11://25:    //D3  --> bug!
    if(ddebug) printf("on create %d %d %d %d -- %d %d %d %d -- %d %d %d %d\n",p0,p4,p3,p5,p0,p1,p4,p5,p5,p1,p2,p0);
    pt = &mesh->tetra[k + 1];
    pt->v[0] = p0;
    pt->v[1] = p4;
    pt->v[2] = p3;
    pt->v[3] = p5;
    pt->ref  = ref;//3;//ref;
    pt = &mesh->tetra[k + 2];
    pt->v[0] = p0;
    pt->v[1] = p1;
    pt->v[2] = p4;
    pt->v[3] = p5;
    pt->ref  = ref;//3;//ref;
    pt = &mesh->tetra[k + 3];
    pt->v[0] = p5;
    pt->v[1] = p1;
    pt->v[2] = p2;
    pt->v[3] = p0;
    pt->ref  = ref;//3;//ref;
    break;
  case 28:    //D2
    pt = &mesh->tetra[k + 1];
    pt->v[0] = p4;
    pt->v[1] = p5;
    pt->v[2] = p1;
    pt->v[3] = p3;
    pt->ref  = ref;//2;//ref;
    pt = &mesh->tetra[k + 2];
    pt->v[0] = p1;
    pt->v[1] = p2;
    pt->v[2] = p5;
    pt->v[3] = p3;
    pt->ref  = ref;//2;//ref;
    pt = &mesh->tetra[k + 3];
    pt->v[0] = p2;
    pt->v[1] = p3;
    pt->v[2] = p1;
    pt->v[3] = p0;
    pt->ref  = ref;//2;//ref;
    break;
  case 35:    //D4 --> ok
    pt = &mesh->tetra[k + 1];
    pt->v[0] = p0;
    pt->v[1] = p4;
    pt->v[2] = p3;
    pt->v[3] = p5;
    pt->ref  = ref;//4;//ref;
    pt = &mesh->tetra[k + 2];
    pt->v[0] = p0;
    pt->v[1] = p4;
    pt->v[2] = p5;
    pt->v[3] = p2;
    pt->ref  = ref;//4;//ref;
    pt = &mesh->tetra[k + 3];
    pt->v[0] = p0;
    pt->v[1] = p2;
    pt->v[2] = p4;
    pt->v[3] = p1;
    pt->ref  = ref;//4;//ref;
    break;
  case 52:
    pt = &mesh->tetra[k + 1];
    pt->v[0] = p2;
    pt->v[1] = p4;
    pt->v[2] = p5;
    pt->v[3] = p3;
    pt->ref  = ref;//6;//ref;
    pt = &mesh->tetra[k + 2];
    pt->v[0] = p2;
    pt->v[1] = p4;
    pt->v[2] = p1;
    pt->v[3] = p3;
    pt->ref  = ref;//6;//ref;
    pt = &mesh->tetra[k + 3];
    pt->v[0] = p3;
    pt->v[1] = p0;
    pt->v[2] = p1;
    pt->v[3] = p2;
    pt->ref  = ref;//6;//ref;
    break;
  case 49:    //D5
    pt = &mesh->tetra[k + 1];
    pt->v[0] = p0;
    pt->v[1] = p4;
    pt->v[2] = p3;
    pt->v[3] = p2;
    pt->ref  = ref;//5;//ref;
    pt = &mesh->tetra[k + 2];
    pt->v[0] = p3;
    pt->v[1] = p2;
    pt->v[2] = p4;
    pt->v[3] = p5;
    pt->ref  = ref;//5;//ref;
    pt = &mesh->tetra[k + 3];
    pt->v[0] = p0;
    pt->v[1] = p2;
    pt->v[2] = p1;
    pt->v[3] = p4;
    pt->ref  = ref;//5;//ref;
    break;
  default:
    //5,7,11,13,15,18,19,21,22,23,26,27,29,30,31,37,39,40,41,42,43,44,45,46,47,50,51,52,53,54,55:
    MMG_cutadd(mesh,hed,icas,k,p0,p1,p2,p3,p4,p5,ref);

    //printf("icas imposssss %d\n",icas);
    return(0);
    //exit(0);
    break;
  }
  return(1);
}
