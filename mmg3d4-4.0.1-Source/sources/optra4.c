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

int MMG_pretreat(pMesh ,pSol ,double ,int *);
int MMG_optlenbdry_iso(pMesh mesh,pSol sol,double declic,int base);

int MMG_optra4(pMesh mesh,pSol sol) {  
  pTetra  pt;  
  double	declicw;
  double	declic;
  int		base,nm,ns,it,maxtou,alert,nw,nwold,k;
  
  pTetra ptet;
  int ncut,nbdry,j,ip,iel[4],*adja,*adja1,iadr,pp[4],ne,i;
  double coor[3],*m,*mip;

  /* optim coquil */
  alert  = 0;
  maxtou = 20;//0;
  it     = 0;
  
  MMG_caltet = ( sol->offset==6 ) ? MMG_caltet_ani:MMG_caltet_iso;
  MMG_caltet2 = ( sol->offset==6 ) ? MMG_caltet2_ani:MMG_caltet2_iso;
    
  for (k=1; k<=mesh->ne; k++) {
    mesh->tetra[k].flag = mesh->flag;
    mesh->tetra[k].qual = MMG_caltet(mesh,sol,k);
  }
  for (k=1; k<=mesh->np; k++) mesh->point[k].flag = mesh->flag;
  declicw = 5./ALPHAD;
  declic  = 1.5/ ALPHAD;
            
  do {
    base = ++mesh->flag;
//MMG_ratio(mesh,sol,NULL); 
    ns = 0;
    if ( !alert && !mesh->info.noswap ) {       
       ns = MMG_cendel(mesh,sol,declic,base);

     if ( ns < 0 ) {
        alert = 1;
    	ns    = -ns;
      }
    }
    nw = 0;
    /*if (!mesh->info.noinsert)  nw = MMG_pretreat(mesh,sol,declicw,&alert);
    */
    /*sur des surfaces courbes, il existe des tetras ayant  4 points de peau avec Q=3*/
    if ( it < 10 )  {
      nwold = nw;
#ifdef REQUIRED
      nw += 0;
#else
      nw += MMG_opttyp(mesh,sol,declicw,&alert,mesh->flag-1);
#endif
	  declicw *= 1.05;
    }
	nm = 0;    
    if (!mesh->info.nomove) {          
      nm = MMG_optlen(mesh,sol,declic,base);  
    }
    
    //if(!mesh->info.nomove && it<2) MMG_optlap(mesh,sol);
    if ( mesh->info.imprim && nw+ns+nm )
      fprintf(stdout,"     %8d IMPROVED  %8d SWAPPED  %8d MOVED\n",nw,ns,nm);     
    }
  //while ( (ns && ((ns > 0.005*mesh->ne /*&& nwold > nw*/) || it < 5)) && ++it < maxtou );
  while ( ns+nm && ++it < maxtou );

   //for LES   
  if(0 &&mesh->info.optles) {
    fprintf(stdout,"--------------------------------------------------------\n");     
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if( !pt->v[0] ) {
	continue;
      }
      pt->qual = MMG_caltet_LES(mesh,sol,k); 
    }
	  
    declic = 1.5/ALPHAD;
    do {
      base = ++mesh->flag;
      ns = 0;
      if ( !alert && !mesh->info.noswap ) {
	//ns = MMG_cendel(mesh,sol,declic,base);
	if ( ns < 0 ) {
	  alert = 1;
	  ns    = -ns;
	}
      }
      nw = 0;
      nm = 0;    
      if (!mesh->info.nomove) {          
	nm = MMG_optlen_LES(mesh,sol,0.3,base);  
      }
      
      //if(!mesh->info.nomove && it<2) MMG_optlap(mesh,sol);
      if ( mesh->info.imprim && nw+ns+nm )
	fprintf(stdout,"     %8d IMPROVED  %8d SWAPPED  %8d MOVED\n",nw,ns,nm);     
    }
    //while ( (ns && ((ns > 0.005*mesh->ne /*&& nwold > nw*/) || it < 5)) && ++it < maxtou );
    while ( ns+nm && ++it < maxtou );
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if( !pt->v[0] ) {
	continue;
      }
      pt->qual = MMG_caltet(mesh,sol,k); 
    }
  }
  

  ////////////////////////////////////////////////////////////////////
  /* /\*split tet if more than one bdry face*\/ */
  /* ncut = 0; */
  /* ne = mesh->ne; */
  /* for(k=1 ; k<=ne ; k++) { */
  /*   pt = &mesh->tetra[k]; */
  /*   if(!pt->v[0]) continue; */
  /*   iadr = 4*(k-1)+1; */
  /*   adja = &mesh->adja[iadr]; */
  /*   nbdry = 0; */
  /*   for(j=0 ; j<4 ; j++) { */
  /*     if(!adja[j]) nbdry++; */
  /*   } */
  /*   if(nbdry>1) { */
  /*     if(k==1814)printf(" %d bdry face for  %d\n",nbdry,k); */
  /*     ncut++; */
      
  /*     for(i=0 ; i<4 ; i++) pp[i] = pt->v[i]; */
  /*     for(j=0 ; j<3 ; j++) { */
  /*       coor[j]=0; */
  /*       for(i=0 ; i<4 ; i++) { */
  /*         coor[j]+=mesh->point[pp[i]].c[j]; */
  /*       } */
  /*       coor[j]*=0.25; */
  /*     } */
  /*     ip = MMG_newPt(mesh,coor); */
  /*     if(sol->offset==1) { */
  /* 	sol->met[ip] = sol->met[pp[0]]; */
  /*     } else { */
  /* 	iadr = (pp[0]-1)*sol->offset + 1; */
  /* 	m  = &sol->met[iadr];	  */
  /* 	iadr = (ip-1)*sol->offset + 1; */
  /* 	mip  = &sol->met[iadr];	   */
  /* 	memcpy(mip,m,sol->offset*sizeof(double)); */
  /*     } */
  /*     if(k==1814) printf("on cree %d : %e %e %e\n",ip,coor[0],coor[1],coor[2]); */
  /*     //create 3 new tets */
  /*     for(j=0 ; j<3 ; j++) { */
  /*       iel[j] = MMG_newElt(mesh); */
  /*       ptet = &mesh->tetra[iel[j]]; */
  /*       memcpy(ptet,pt,sizeof(Tetra)); */
  /*       ptet->v[j]=ip; */
  /* 	for(i=0 ; i<6 ; i++) ptet->bdryinfo[i] = 0; */
  /* 	ptet->bdryref[0] = ptet->bdryref[1] = ptet->bdryref[2] = ptet->bdryref[3] = 0; */
  /* 	ptet->bdryref[j] = pt->bdryref[j]; */
  /* 	pt->bdryref[j] = 0; */
  /* 	if(k==1814) printf("on cree %d : %d %d %d %d\n",iel[j],ptet->v[0],ptet->v[1],ptet->v[2],ptet->v[3]); */
  /*     } */
  /*     iel[3] = k; */
  /*     pt->v[3] = ip; */
  /*     for(i=0 ; i<6 ; i++) pt->bdryinfo[i] = 0; */
  /*     //adja of iel[0] */
  /*     iadr = 4*(iel[0]-1)+1; */
  /*     adja1 = &mesh->adja[iadr]; */
  /*     adja1[0] = adja[0]; */
  /*     (&mesh->adja[4*(adja[0]/4-1)+1])[adja[0]%4] = 4*iel[0] + 0; */
  /*     for(i=1 ; i<4 ; i++) adja1[i]=4*iel[i]+0; */

  /*     //adja of iel[1] */
  /*     iadr = 4*(iel[1]-1)+1; */
  /*     adja1 = &mesh->adja[iadr]; */
  /*     adja1[1] = adja[1]; */
  /*     (&mesh->adja[4*(adja[1]/4-1)+1])[adja[1]%4] = 4*iel[1] + 1; */
  /*     for(i=0 ; i<4 ; i++) { */
  /* 	if(i==1) continue; */
  /* 	adja1[i]=4*iel[i]+1; */
  /*     } */

  /*     //adja of iel[2] */
  /*     iadr = 4*(iel[2]-1)+1; */
  /*     adja1 = &mesh->adja[iadr]; */
  /*     adja1[2] = adja[2]; */
  /*     (&mesh->adja[4*(adja[2]/4-1)+1])[adja[2]%4] = 4*iel[2] + 2; */
  /*     for(i=0 ; i<4 ; i++) { */
  /* 	if(i==2) continue; */
  /* 	adja1[i]=4*iel[i]+2; */
  /*     } */

  /*     //adj of k */
  /*     for(i=0 ; i<3 ; i++) adja[i]=4*iel[i]+3; */
  /*   } */
  /* }   */
  /* fprintf(stdout," %8d TET SPLITTED\n",ncut); */
  /* MMG_chkmsh(mesh,0,0); */
  ////////////////////////////////////////////////////////////////////////////////////

  return(1);
}

int MMG_optra4_LES(pMesh mesh,pSol sol) {  
  pTetra  pt;  
  double	declicw;
  double	declic;
  int		base,nm,ns,it,maxtou,alert,nw,nwold,k;
  
  pTetra ptet;
  int ncut,nbdry,j,ip,iel[4],*adja,*adja1,iadr,pp[4],ne,i;
  double coor[3],*m,*mip;

  if(sol->offset==6) {
    fprintf(stdout,"THIS OPTIMIZATION IS NOT SUITABLE FOR ANISOTROPIC MESHES\n");
    return(1);
  }
  /* optim coquil */
  alert  = 0;
  maxtou = 10;//0;
  it     = 0;
  
  MMG_caltet = MMG_caltet_LES;
  MMG_caltet2 = MMG_caltet2_LES;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if( !pt->v[0] ) {
      continue;
    }
    pt->flag = mesh->flag;
    pt->qual = MMG_caltet(mesh,sol,k);
  }
  
  for (k=1; k<=mesh->np; k++) mesh->point[k].flag = mesh->flag;
  declicw = 0.7;//5./ALPHAD; 0.95 not enough 
  
  //try to delete bad element
  do {
    base = ++mesh->flag;
#ifdef REQUIRED
    nw = 0;
#else
    nw = MMG_opttyp(mesh,sol,declicw,&alert,-1);
#endif
    fprintf(stdout,"     %8d IMPROVED \n",nw);
  } while ( nw && ++it < 4 );
   

  it     = 0;
  declic  = 0.5;
  declicw = 0.8;
	  
  do {
    ++mesh->flag;
    base = -1; //avec base=mesh->flag-3 apres la 1e iter on est bien plus rapide mais il reste qq mauvais..
    //la proc de bouge de point n'est pas tres performante, il faut 10 iter pour bouger le point....
   
    ns = 0;
    if ( !alert && !mesh->info.noswap  ) {       
      ns = MMG_cendel_LES(mesh,sol,declic,base);
      if ( ns < 0 ) {
        alert = 1;
    	ns    = -ns;
      }
    }
    nm = 0;    
    if (!mesh->info.nomove ) {          
      nm = MMG_optlen_LES(mesh,sol,declic,-1);  
    }
    if(nw && it<5) {
#ifdef REQUIRED
      nw = 0;
#else
      nw = MMG_opttyp(mesh,sol,declicw,&alert,-1);
#endif
      if(nw) fprintf(stdout,"     %8d IMPROVED \n",nw);
    }
    //optimise en moy mais peut degrader un peu le plus mauvais
    if(!mesh->info.nomove && it==2) MMG_optlap(mesh,sol); 
    if ( mesh->info.imprim && ns+nm )
      fprintf(stdout,"    %8d SWAPPED  %8d MOVED\n",ns,nm);    
    }
  while ( (ns+nw && ((ns+nw > 0.005*mesh->ne ) || it < 5)) && ++it < maxtou );
  //while ( ns+nm && ++it <maxtou );

  MMG_caltet = ( sol->offset==6 ) ? MMG_caltet_ani:MMG_caltet_iso;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if( !pt->v[0] ) {
      continue;
    }
    pt->qual = MMG_caltet(mesh,sol,k); 
  }


  return(1);
}

