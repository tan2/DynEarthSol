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


double MMG_rao(pMesh mesh,int k,int inm);
double MMG_caltetrao(pMesh mesh,pSol sol,int iel) {
	return(MMG_rao(mesh,iel,0));
}

double MMG_isoperimetric(pMesh mesh) {
  pTria ptri;
  pTetra ptet;
  pPoint    pa,pb,pc;
  double area, vol, dd, s;
  int ia, ib, ic, tetraIdx, triaIdx;
  double abx, aby, abz, acx, acy, acz, bcx, bcy, bcz, ab, bc, ac;

  vol = 0;
  /* Compute volume for all tetrahedra */
  for(tetraIdx = 1 ; tetraIdx < mesh->ne + 1 ; tetraIdx++) {
    ptet = &mesh->tetra[tetraIdx];

    /* Testing if the tetra exists */
    if (!ptet->v[0]) continue;

    vol += MMG_quickvol(mesh->point[ptet->v[0]].c,
        mesh->point[ptet->v[1]].c,
        mesh->point[ptet->v[2]].c,
        mesh->point[ptet->v[3]].c) / 6.0;
  }

  area = 0;
  for(triaIdx = 1 ; triaIdx < mesh->nt + 1 ; triaIdx++) {
    ptri = &mesh->tria[triaIdx];

    /* Testing if the triangle exists */
    if (!ptri->v[0]) continue;

    ia = ptri->v[0];
    ib = ptri->v[1];
    ic = ptri->v[2];
    pa = &mesh->point[ia];
    pb = &mesh->point[ib];
    pc = &mesh->point[ic];

    /* area */
    abx = pb->c[0] - pa->c[0]; 
    aby = pb->c[1] - pa->c[1]; 
    abz = pb->c[2] - pa->c[2]; 

    bcx = pc->c[0] - pb->c[0]; 
    bcy = pc->c[1] - pb->c[1]; 
    bcz = pc->c[2] - pb->c[2]; 

    acx = pc->c[0] - pa->c[0]; 
    acy = pc->c[1] - pa->c[1]; 
    acz = pc->c[2] - pa->c[2]; 

    ab = sqrt (abx * abx + aby * aby + abz * abz);
    bc = sqrt (bcx * bcx + bcy * bcy + bcz * bcz);
    ac = sqrt (acx * acx + acy * acy + acz * acz);

    s = (ab + bc + ac) / 2.0;
    //dd = aby*acz - abz*acy;
    //s = dd * dd;
    //dd = abz*acx - abx*acz;
    //s += dd * dd;
    //dd = abx*acy - aby*acx;
    //s += dd * dd;
    area += sqrt(s * (s - ab) * (s - ac) * (s -bc));
  }

  printf("volume: %lf, area: %lf\n, quotient: %lf", vol, area, (double) ((36.0 * M_PI * vol * vol) / (area * area * area)));
  return ((double) ((36.0 * M_PI * vol * vol) / (area * area * area)));
}

/*approximation of the final number of vertex*/
int MMG_count(pMesh mesh,pSol sol, double *weightelt, long *npcible) {
  pTetra pt;
  double *ca,*cb,*ma,*mb,len;
  int    k,ia,ipa,ipb,iad,lon,l,iadr;
  int    *pdel,lenint,loc,nedel,longen;
  int    npbdry,isbdry;
  double   dned,dnface,dnint,dnins,w;
  double   dnpdel,dnadd,leninv,dnaddloc,dnpdelloc;
  List   list;
  int ib;
  double lenavg,lent[6];
  long nptot;
  pdel = (int*) calloc(mesh->np + 1,sizeof(int));
  nptot = (long) mesh->np;
  npbdry = 0;

  //substraction of the half of the number of bdry vertex to avoid the surestimation due of the interface
  //for (k=1; k<=mesh->np; k++) {
  //  if(mesh->point[k].tag & M_BDRY) npbdry++;
  //}
  //nptot -= 0.5*npbdry;

  dnadd = dnpdel = 0;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !pt->v[0] )  continue;

    /*longueur moyenne*/
    lenavg = 0;
    for(ib=0 ; ib<6 ; ib++) {
      ipa = MMG_iare[ib][0];
      ipb = MMG_iare[ib][1];
      ca  = &mesh->point[pt->v[ipa]].c[0];
      cb  = &mesh->point[pt->v[ipb]].c[0];

      iadr = (pt->v[ipa]-1)*sol->offset + 1;
      ma   = &sol->met[iadr];
      iadr = (pt->v[ipb]-1)*sol->offset + 1;
      mb   = &sol->met[iadr];
      if(sol->offset==6)
	lent[ib] = MMG_long_ani_init(ca,cb,ma,mb);
      else
	lent[ib] = MMG_length(ca,cb,ma,mb);
      lenavg+=lent[ib];
    }
    lenavg /= 6.;
    
    w = 0;
    if(weightelt)
      weightelt[k] = 0;
    nedel = 0;
    for (ia=0; ia<6; ia++) {
      //lon = MMG_coquil(mesh,k,ia,&list);
      longen = MMG_coquilgen(mesh,k,ia,&list);
      lon = longen/2;
      isbdry = 0;//longen%2;
      if(!lon) continue;
      /* if ( isbdry )  { */
      /* 	assert(longen%2); */
      /* 	//printf("coquil %d\n",longen/2); */
      /* 	continue; */
      /* } */
      //assert(!(longen%2));
      for (l=2; l<=lon; l++) {
	if ( list.tetra[l] < 6*k )  break;
      }

      if ( l <= lon )  {
	loc = 1;
	//continue;
      } else {
	loc = 0;
      } 

      dnaddloc = 0;
      dnpdelloc = 0;

      len = lent[ia];
      
      if(len > 3) {
        loc = 0; //CECILE : on compte toutes les aretes et on pondere par lon
	//si la longueur est tres gde, on prend la valeur moyenne
	len = lenavg;
	lenint = ((int) len); 
	if(fabs(lenint -len) > 0.5) lenint++;
	//POURQUOI SURESTIMER ???lenint++;
	//nb de point a inserer sur l'arete : longueur - 1
	dned = lenint - 1;
	//nb de point a l'interieur de la face si toutes les aretes sont coupees le meme nb de fois
	dnface = (lenint+2)*(lenint+1) / 2. - 3 - 3*dned;
	//nb de point a l'interieur du tetra si toutes les aretes sont coupees le meme nb de fois
	dnint = (lenint+3)*(lenint+2)*(lenint+1) / 6. - 4 - 4*dnface - 6*dned;
	//nb de point a inserer pour cette arete de ce tetra : on divise par lon
	dnins = dned*(1./lon) + (dnface/3. + dnint/6.);//(dnface/12. + dnint/6.);
	//if(!loc) {
	  if(!isbdry) {
	    //nb points sur l'arete + 
	    //lon*(2/3 nb point sur la face (ie 1/3 de face et 2 faces adj a l'arete) + 1/6 nb de point interne)
	    dnaddloc = dned + lon*(2*dnface/3. + dnint/6.);
	  } else {
	    dnaddloc = 0.5*(dned + lon*(2*dnface/3. + dnint/6.));
	  }
        dnaddloc *= 1./lon;
	//if(dnint > 10 && (ALPHAD * pt->qual < 10))	printf("dnint %d %e %e\n",(int) dnint,len,ALPHAD * pt->qual);
        if(!loc) {
          if(ALPHAD * pt->qual < 2) /*on ne compte les points internes que pour les (tres) bons tetras*/
            dnaddloc = dnaddloc;
          else if(ALPHAD * pt->qual < 5) 
	    dnaddloc = dned / lon + 2*dnface/3.;
	  else
            dnaddloc = dned / lon ;
	  //rajout de 30% parce que 1) on vise des longueurs de 0.7 et 
	  //2) on ne tient pas compte du fait qu'on divise tjs par 2 dans la generation
	  if( (ALPHAD*pt->qual > 50) )
	    dnaddloc = 0;
	  else  if((ALPHAD*pt->qual > 10) )
	    dnaddloc =  0.2*dnaddloc; //CEDRIC : essayer 0.3 0.4
	  else if((len > 10) && (ALPHAD*pt->qual < 1.5) ) //on sous-estime uniquement pour les tres bons
	    dnaddloc = dnaddloc*0.3 + dnaddloc; //CEDRIC : essayer 0.3 ?
	  else if(len < 6 && len>3) //CEDRIC : essayer len < 3,4, 6,7 mais aussi en commentant le test sur len puis pour la qual > 3, 5,8
	    dnaddloc = 0.7*dnaddloc; //CEDRIC : essayer 0.9 0.7 0.6


	    dnadd += dnaddloc;
	}
	//printf("len %d on a ajoute %e : %e -- %e %e %e tot %e == %e\n",lenint,dnaddloc,dnadd,dned,dnface,dnint,dned+2*dnface/3+dnint/6,(6*dned+4*dnface+dnint)/6);
        //}
      } else if(len > 2.8) {
	if(!isbdry) {
	  dnaddloc = 2.;
	} else {
	  dnaddloc = 1;
	}
	if(!loc){
	  if(!isbdry) {
	    dnadd += 2.;
	  } else {
	    dnadd++;
	  }
	}
	dnins = 2;
      } else if(len > 1.41) {
	if(!isbdry)
	  dnaddloc = 1;
	if(!loc) {
	  if(!isbdry) dnadd += 1.;
	}
	dnins = 1;
      } else if(len < 0.6) {
	nedel = 1;
	  leninv = 1./len;
	  if(pt->v[ipa]<pt->v[ipb]) {
	    if(!pdel[pt->v[ipa]]) {
	      if(!isbdry) {
              dnpdelloc = (leninv - 1.)/leninv;
	      } else {
              dnpdelloc = 0.5*(leninv - 1.)/leninv;
	      }
            if(!loc) {

              dnpdel+=dnpdelloc;
	      pdel[pt->v[ipa]]=1;
            }
	    } else if(!pdel[pt->v[ipb]]) {
	      if(!isbdry) {
		dnpdelloc = (leninv - 1.)/leninv;
	      } else {
		dnpdelloc = 0.5*(leninv - 1.)/leninv;
	      }
              if (!loc) {
                dnpdel += dnpdelloc;
	      pdel[pt->v[ipb]]=1;	  
	    }
	    }
	  } else {
	    if(!pdel[pt->v[ipb]]) {
	      if(!isbdry) {
		dnpdelloc = (leninv - 1.)/leninv;
	      } else {
		dnpdelloc = 0.5*(leninv - 1.)/leninv;
	      }
              if (!loc) {
                dnpdel += dnpdelloc;
	      pdel[pt->v[ipb]]=1;
              }
	    } else if(!pdel[pt->v[ipa]]) {
	      if(!isbdry) {
		dnpdelloc = (leninv - 1.)/leninv;
	      } else {
		dnpdelloc = 0.5*(leninv - 1.)/leninv;
	      }
              if (!loc) {
                dnpdel += dnpdelloc;
	      pdel[pt->v[ipa]]=1;
	    }
	  }
	}
      }

      //pour cette arete de ce tetra :
      //PHASE 1 = dnaddloc + nedel (on compte un si arete trop petite)
      //PHASE 2 = dnaddloc
      w += (2*dnaddloc);//1./lon*(2*dnaddloc + dnpdelloc);
    }/*for ia*/

   // if(weightelt)
    //  weightelt[k] = dnins*7. + 6.*nedel;

    w += 0.5*nedel;
    if(weightelt)
      weightelt[k] = 10*w;
 
  } /*For k*/


  nptot += (long) dnadd - (long) dnpdel;
  *npcible = nptot;
  fprintf(stdout,"ESTIMATION OF THE FINAL NUMBER OF NODES : %ld  ADD %f  DEL %f\n",nptot,dnadd,dnpdel);
	
  free(pdel);
  return(1);
}
/* compute tetra volume */
double MMG_voltet(pMesh mesh,int iel) {
  pTetra   pt;
  pPoint   p1,p2,p3,p4;
  double   ax,ay,az,bx,by,bz,vol,len;
  int      s1,s2,s3,s4;

  pt = &mesh->tetra[iel];
  s1 = pt->v[0];
  s2 = pt->v[1];
  s3 = pt->v[2];
  s4 = pt->v[3];

  /* sort tetra vertices */
  if ( pt->v[0] < pt->v[1] ) {
    if ( pt->v[0] < pt->v[2] ) {
      s3 = pt->v[2];
      if ( pt->v[0] < pt->v[3] ) {
        s1 = pt->v[0];
	      s2 = pt->v[1];
	      s4 = pt->v[3];
      }
      else {
        s1 = pt->v[3];
	      s2 = pt->v[0];
	      s4 = pt->v[1];
      }
    }
    else {
      s2 = pt->v[0];
      if ( pt->v[2] < pt->v[3] ) {
        s1 = pt->v[2];
	      s3 = pt->v[1];
	      s4 = pt->v[3];
      }
      else {
        s1 = pt->v[3];
	      s3 = pt->v[2];
	      s4 = pt->v[1];
      }
    }
  }
  else {
    if ( pt->v[1] < pt->v[2] ) {
      if ( pt->v[1] < pt->v[3] ) {
        s1 = pt->v[1];
	      s2 = pt->v[2];
        s3 = pt->v[0];
	      s4 = pt->v[3];
      }
      else {
        s1 = pt->v[3];
	      s2 = pt->v[0];
        s3 = pt->v[2];
	      s4 = pt->v[1];
      }
    }
    else {
      s2 = pt->v[0];
      if ( pt->v[2] < pt->v[3] ) {
        s1 = pt->v[2];
        s3 = pt->v[1];
        s4 = pt->v[3];
      }
      else {
        s1 = pt->v[3];
        s3 = pt->v[2];
        s4 = pt->v[1];
      }
    }
  }

  p1 = &mesh->point[s1];
  p2 = &mesh->point[s4];
  p3 = &mesh->point[s2];
  p4 = &mesh->point[s3];
  if ( s2 < s3 ) {
    if ( s2 < s4 ) {
      p2 = &mesh->point[s2];
      p3 = &mesh->point[s3];
      p4 = &mesh->point[s4];
    }
  }
  else if ( s3 < s4 ) {
    p2 = &mesh->point[s3];
    p3 = &mesh->point[s4];
    p4 = &mesh->point[s2];
  }

  /* compute volume */
  ax = p3->c[0] - p1->c[0];
  ay = p3->c[1] - p1->c[1];
  az = p3->c[2] - p1->c[2];

  bx = p4->c[0] - p1->c[0];
  by = p4->c[1] - p1->c[1];
  bz = p4->c[2] - p1->c[2];
  //printf("a %e %e %e \n",(ay*bz - az*by),(az*bx - ax*bz),(ax*by - ay*bx));
  vol = (p2->c[0]-p1->c[0]) * (ay*bz - az*by) \
      + (p2->c[1]-p1->c[1]) * (az*bx - ax*bz) \
      + (p2->c[2]-p1->c[2]) * (ax*by - ay*bx);
  //printf("%e %e %e\n",(p2->c[0]-p1->c[0]),(p2->c[1]-p1->c[1]),(p2->c[2]-p1->c[2]));
  //printf("vol1 %e %e %e\n",(p2->c[0]-p1->c[0]) * (ay*bz - az*by),(p2->c[1]-p1->c[1]) * (az*bx - ax*bz),
  //  (p2->c[2]-p1->c[2]) * (ax*by - ay*bx));
  len = ax*ax+ay*ay+az*az;
  len += bx*bx+by*by+bz*bz;
  len += (p2->c[0]-p1->c[0])*(p2->c[0]-p1->c[0])+(p2->c[1]-p1->c[1])*(p2->c[1]-p1->c[1])
       + (p2->c[2]-p1->c[2])*(p2->c[2]-p1->c[2]);
  len += (p2->c[0]-p3->c[0])*(p2->c[0]-p3->c[0])+(p2->c[1]-p3->c[1])*(p2->c[1]-p3->c[1])
       + (p2->c[2]-p3->c[2])*(p2->c[2]-p3->c[2]);
  len += (p2->c[0]-p4->c[0])*(p2->c[0]-p4->c[0])+(p2->c[1]-p4->c[1])*(p2->c[1]-p4->c[1])
       + (p2->c[2]-p4->c[2])*(p2->c[2]-p4->c[2]);
  len += (p3->c[0]-p4->c[0])*(p3->c[0]-p4->c[0])+(p3->c[1]-p4->c[1])*(p3->c[1]-p4->c[1])
       + (p3->c[2]-p4->c[2])*(p3->c[2]-p4->c[2]);    
  len /= 6.;
  len = sqrt(len);
  len = len*len*len; 
  //printf("vol %e %e\n",vol,len);       
  vol /= len;
  return(vol);
}


/* quick volume check */
double MMG_quickvol(double *c1,double *c2,double *c3,double *c4) {
  double   ax,ay,az,bx,by,bz,vol;

  ax = c3[0] - c1[0];
  ay = c3[1] - c1[1];
  az = c3[2] - c1[2];
  
  bx = c4[0] - c1[0];
  by = c4[1] - c1[1];
  bz = c4[2] - c1[2];
  
  vol = (c2[0]-c1[0]) * (ay*bz - az*by) \
      + (c2[1]-c1[1]) * (az*bx - ax*bz) \
      + (c2[2]-c1[2]) * (ax*by - ay*bx);

  return(vol);
}

//Q = Vref-V/Vref
//Vref = vol tet equilateral compris dans sphere cirsconstrite a iel
//Vref = sqrt(2)/12 * a**3 avec a = 2Rsqrt(6)/3
double MMG_caltet_LES(pMesh mesh,pSol sol,int iel) {
  pTetra     pt;
  double    *a,*b,*c,*d,ct[12],cs[3],rad,lref,Vref,V,cal;
  int        ia,ib,ic,id,j,l;

  cal = CALLIM;
  pt = &mesh->tetra[iel];  
  if ( !pt->v[0] )  return(cal);

  ia = pt->v[0];
  ib = pt->v[1];
  ic = pt->v[2];
  id = pt->v[3];
  a  = mesh->point[ia].c;
  b  = mesh->point[ib].c;
  c  = mesh->point[ic].c;
  d  = mesh->point[id].c; 
  for (j=0,l=0; j<4; j++,l+=3) {
    memcpy(&ct[l],mesh->point[pt->v[j]].c,3*sizeof(double));  
  }

  if(!MMG_cenrad_iso(mesh,ct,cs,&rad)) {
    //printf("pbs calcul rayon sphere\n");
    return(cal);
  }

  if(rad < 0) return(2);
  rad = sqrt(rad);
  /* volume Vref*/
  lref = (2.*rad*sqrt(6.))/3.;
  Vref = lref*lref*lref*sqrt(2.)/12.;
  V = MMG_quickvol(a,b,c,d)/6.; 
  if(V<=0) return(2);
  cal = (Vref - V) / Vref; 
  return(cal);
} 


/* compute tetra quality aniso */
double MMG_caltet_ani(pMesh mesh,pSol sol,int iel) {
  pTetra     pt;
  double     cal,abx,aby,abz,acx,acy,acz,adx,ady,adz;
  double     bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz;
  double     h1,h2,h3,h4,h5,h6,det,vol,rap,v1,v2,v3,num;
  double    *a,*b,*c,*d;
  double     *ma,*mb,*mc,*md,mm[6];
  int        j,ia,ib,ic,id,iadr;

  cal = CALLIM;
  pt  = &mesh->tetra[iel];
  if ( !pt->v[0] )  return(cal);

  ia = pt->v[0];
  ib = pt->v[1];
  ic = pt->v[2];
  id = pt->v[3];

  /* average metric */
  memset(mm,0,6*sizeof(double));
  iadr = (ia-1)*sol->offset + 1;
  ma   = &sol->met[iadr];
  iadr = (ib-1)*sol->offset + 1;
  mb   = &sol->met[iadr];
  iadr = (ic-1)*sol->offset + 1;
  mc   = &sol->met[iadr];
  iadr = (id-1)*sol->offset + 1;
  md   = &sol->met[iadr];
  for (j=0; j<6; j++)
    mm[j] = 0.25 * (ma[j]+mb[j]+mc[j]+md[j]);
  a = mesh->point[ia].c;
  b = mesh->point[ib].c;
  c = mesh->point[ic].c;
  d = mesh->point[id].c;

  /* volume */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];

  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];
  
  adx = d[0] - a[0];
  ady = d[1] - a[1];
  adz = d[2] - a[2];

  v1  = acy*adz - acz*ady;
  v2  = acz*adx - acx*adz;
  v3  = acx*ady - acy*adx;
  vol = abx * v1 + aby * v2 + abz * v3;            
  if ( vol <= 0. )  return(cal);
  det = mm[0] * ( mm[3]*mm[5] - mm[4]*mm[4]) \
      - mm[1] * ( mm[1]*mm[5] - mm[2]*mm[4]) \
      + mm[2] * ( mm[1]*mm[4] - mm[2]*mm[3]);   
  if ( det < EPSOK )   {
    //printf("--- INVALID METRIC : DET (%d) %e\n",iel,det);
    return(cal);
  }
  det = sqrt(det) * vol;
  /* edge lengths */
  h1 =      mm[0]*abx*abx + mm[3]*aby*aby + mm[5]*abz*abz \
     + 2.0*(mm[1]*abx*aby + mm[2]*abx*abz + mm[4]*aby*abz);
  h2 =      mm[0]*acx*acx + mm[3]*acy*acy + mm[5]*acz*acz \
     + 2.0*(mm[1]*acx*acy + mm[2]*acx*acz + mm[4]*acy*acz);
  h3 =      mm[0]*adx*adx + mm[3]*ady*ady + mm[5]*adz*adz \
     + 2.0*(mm[1]*adx*ady + mm[2]*adx*adz + mm[4]*ady*adz);

  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];

  bdx = d[0] - b[0];
  bdy = d[1] - b[1];
  bdz = d[2] - b[2];

  cdx = d[0] - c[0];
  cdy = d[1] - c[1];
  cdz = d[2] - c[2];

  h4 =      mm[0]*bdx*bdx + mm[3]*bdy*bdy + mm[5]*bdz*bdz \
     + 2.0*(mm[1]*bdx*bdy + mm[2]*bdx*bdz + mm[4]*bdy*bdz);
  h5 =      mm[0]*cdx*cdx + mm[3]*cdy*cdy + mm[5]*cdz*cdz \
     + 2.0*(mm[1]*cdx*cdy + mm[2]*cdx*cdz + mm[4]*cdy*cdz);
  h6 =      mm[0]*bcx*bcx + mm[3]*bcy*bcy + mm[5]*bcz*bcz \
     + 2.0*(mm[1]*bcx*bcy + mm[2]*bcx*bcz + mm[4]*bcy*bcz);

  /* quality */
  rap = h1 + h2 + h3 + h4 + h5 + h6;
  num = sqrt(rap) * rap;  
  //if ( det > EPS1 * num )  completement arbitraire comme seuil!!!!
  cal = num / det;  
  if(cal >= CALLIM) printf(" %d %e %e %e %e\n",iel,cal,num,det,vol);  
  
  assert(cal < CALLIM);
  return(cal);
}

/* compute tetra quality iso */
double MMG_caltet5_iso(pMesh mesh,pSol sol,int iel) {
  pTetra     pt;
  double     cal,abx,aby,abz,acx,acy,acz,adx,ady,adz;
  double     bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz;
  double     h1,h2,h3,h4,h5,h6,vol,v1,v2,v3,rap,num;
  double    *a,*b,*c,*d;
  int        ia,ib,ic,id;
double h;

  cal = CALLIM;
  pt = &mesh->tetra[iel];
  if ( !pt->v[0] )  return(cal);

  ia = pt->v[0];
  ib = pt->v[1];
  ic = pt->v[2];
  id = pt->v[3];
  a  = mesh->point[ia].c;
  b  = mesh->point[ib].c;
  c  = mesh->point[ic].c;
  d  = mesh->point[id].c;
  
  h1 = sol->met[ia];
  h2 = sol->met[ib];
  h3 = sol->met[ic];
  h4 = sol->met[id];
  h  = 0.25*(h1 + h2 + h3 + h4);

  /* volume */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];

  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];
  
  adx = d[0] - a[0];
  ady = d[1] - a[1];
  adz = d[2] - a[2];

  v1  = acy*adz - acz*ady;
  v2  = acz*adx - acx*adz;
  v3  = acx*ady - acy*adx;
  vol = abx * v1 + aby * v2 + abz * v3;
  if ( vol <= 0. )  return(cal);

  /* max edge */
  h1 = abx*abx + aby*aby + abz*abz;
  h2 = acx*acx + acy*acy + acz*acz;
  h3 = adx*adx + ady*ady + adz*adz;

  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];

  bdx = d[0] - b[0];
  bdy = d[1] - b[1];
  bdz = d[2] - b[2];

  cdx = d[0] - c[0];
  cdy = d[1] - c[1];
  cdz = d[2] - c[2];

  h4 = bdx*bdx + bdy*bdy + bdz*bdz;
  h5 = cdx*cdx + cdy*cdy + cdz*cdz;
  h6 = bcx*bcx + bcy*bcy + bcz*bcz;

  /* quality */
  rap = h1*h1 + h2*h2 + h3*h3 + h4*h4 + h5*h5 + h6*h6;
  num = sqrt(rap) * rap;
  cal = num / vol;
  assert( cal < CALLIM );

  return(cal);
}          


/* compute tetra quality iso */
double MMG_caltet_iso(pMesh mesh,pSol sol,int iel) {
  pTetra     pt;
  double     cal,abx,aby,abz,acx,acy,acz,adx,ady,adz;
  double     bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz;
  double     h1,h2,h3,h4,h5,h6,vol,v1,v2,v3,rap,num;
  double    *a,*b,*c,*d;
  int        ia,ib,ic,id;

  cal = CALLIM;
  pt = &mesh->tetra[iel];  
  if ( !pt->v[0] )  return(cal);

  ia = pt->v[0];
  ib = pt->v[1];
  ic = pt->v[2];
  id = pt->v[3];
  a  = mesh->point[ia].c;
  b  = mesh->point[ib].c;
  c  = mesh->point[ic].c;
  d  = mesh->point[id].c;
  
  /* volume */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];

  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];
  
  adx = d[0] - a[0];
  ady = d[1] - a[1];
  adz = d[2] - a[2];

  v1  = acy*adz - acz*ady;
  v2  = acz*adx - acx*adz;
  v3  = acx*ady - acy*adx;
  vol = abx * v1 + aby * v2 + abz * v3;   
  if ( vol <= 0. )  {
      /*tmp = pt->v[2];
      pt->v[2] =pt->v[1];
      pt->v[1] = tmp;*/  
      //printf("arg on a un vol nul!!%8e pour %d\n ",vol,iel);
      return(cal);}

  /* max edge */
  h1 = abx*abx + aby*aby + abz*abz;
  h2 = acx*acx + acy*acy + acz*acz;
  h3 = adx*adx + ady*ady + adz*adz;

  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];

  bdx = d[0] - b[0];
  bdy = d[1] - b[1];
  bdz = d[2] - b[2];

  cdx = d[0] - c[0];
  cdy = d[1] - c[1];
  cdz = d[2] - c[2];

  h4 = bdx*bdx + bdy*bdy + bdz*bdz;
  h5 = cdx*cdx + cdy*cdy + cdz*cdz;
  h6 = bcx*bcx + bcy*bcy + bcz*bcz;

  /* quality */
  rap = h1 + h2 + h3 + h4 + h5 + h6;
  num = sqrt(rap) * rap;
  //if ( vol > EPS1 * num ) completement arbitraire comme seuil!!!!    
  cal = num / vol;
  assert( cal < CALLIM );

  return(cal);
}


static double determ(double ab[3],double ac[3],double *mm) {
  double    a1,a2,a3,a4,a5,a6,m0,m1,m2,m3;

  a1 = mm[0]*ab[0] + mm[1]*ab[1] + mm[2]*ab[2];
  a2 = mm[1]*ab[0] + mm[3]*ab[1] + mm[4]*ab[2];
  a3 = mm[2]*ab[0] + mm[4]*ab[1] + mm[5]*ab[2];

  a4 = mm[0]*ac[0] + mm[1]*ac[1] + mm[2]*ac[2];
  a5 = mm[1]*ac[0] + mm[3]*ac[1] + mm[4]*ac[2];
  a6 = mm[2]*ac[0] + mm[4]*ac[1] + mm[5]*ac[2];

  m0 = a1*ab[0] + a2*ab[1] + a3*ab[2];
  m1 = a4*ab[0] + a5*ab[1] + a6*ab[2];
  m2 = a1*ac[0] + a2*ac[1] + a3*ac[2];
  m3 = a4*ac[0] + a5*ac[1] + a6*ac[2];

  return(m0*m3 - m1*m2);
}


/* compute tetra quality for output (same as GHS3D)
   note: must be normalized by ALPHAC  */
double MMG_calte1_ani(pMesh mesh,pSol sol,int iel) {
  pTetra     pt;
  double     cal,ab[3],ac[3],ad[3],bc[3],bd[3],cd[3];
  double     h1,h2,h3,h4,h5,h6,rapmax,vol,det,v1,v2,v3;
  double     air,dd,num;
  double    *a,*b,*c,*d;
  double     *ma,*mb,*mc,*md,mm[6];;
  int        j,ia,ib,ic,id,iadr;

  pt = &mesh->tetra[iel];
  if ( !pt->v[0] )  return(0.0);
  cal = CALLIM;

  ia = pt->v[0];
  ib = pt->v[1];
  ic = pt->v[2];
  id = pt->v[3];

  /* average metric */
  memset(mm,0,6*sizeof(double));
  iadr = (ia-1)*sol->offset + 1;
  ma   = &sol->met[iadr];
  iadr = (ib-1)*sol->offset + 1;
  mb   = &sol->met[iadr];
  iadr = (ic-1)*sol->offset + 1;
  mc   = &sol->met[iadr];
  iadr = (id-1)*sol->offset + 1;
  md   = &sol->met[iadr];
  for (j=0; j<6; j++)
    mm[j] = 0.25 * (ma[j]+mb[j]+mc[j]+md[j]);

  a  = mesh->point[ia].c;
  b  = mesh->point[ib].c;
  c  = mesh->point[ic].c;
  d  = mesh->point[id].c;

  /* volume */
  ab[0] = b[0] - a[0];
  ab[1] = b[1] - a[1];
  ab[2] = b[2] - a[2];

  ac[0] = c[0] - a[0];
  ac[1] = c[1] - a[1];
  ac[2] = c[2] - a[2];
  
  ad[0] = d[0] - a[0];
  ad[1] = d[1] - a[1];
  ad[2] = d[2] - a[2];

  v1  = ac[1]*ad[2] - ac[2]*ad[1];
  v2  = ac[2]*ad[0] - ac[0]*ad[2];
  v3  = ac[0]*ad[1] - ac[1]*ad[0];
  vol = ab[0] * v1 + ab[1] * v2 + ab[2] * v3;
  if ( vol <= 0. )  return(cal);

  det = mm[0] * ( mm[3]*mm[5] - mm[4]*mm[4]) \
      - mm[1] * ( mm[1]*mm[5] - mm[2]*mm[4]) \
      + mm[2] * ( mm[1]*mm[4] - mm[2]*mm[3]);
  if ( det < EPSOK )  {
    if(iel) printf("INVALID METRIC : DET (%d) %e\n",iel,det);
    return(cal);
  }
  det = sqrt(det) * vol;

  /* edge lengths */
  rapmax = 0.0;
  h1 =      mm[0]*ab[0]*ab[0] + mm[3]*ab[1]*ab[1] + mm[5]*ab[2]*ab[2] \
     + 2.0*(mm[1]*ab[0]*ab[1] + mm[2]*ab[0]*ab[2] + mm[4]*ab[1]*ab[2]);
  rapmax = M_MAX(rapmax,h1); 
  h2 =      mm[0]*ac[0]*ac[0] + mm[3]*ac[1]*ac[1] + mm[5]*ac[2]*ac[2] \
     + 2.0*(mm[1]*ac[0]*ac[1] + mm[2]*ac[0]*ac[2] + mm[4]*ac[1]*ac[2]);
  rapmax = M_MAX(rapmax,h2); 
  h3 =      mm[0]*ad[0]*ad[0] + mm[3]*ad[1]*ad[1] + mm[5]*ad[2]*ad[2] \
     + 2.0*(mm[1]*ad[0]*ad[1] + mm[2]*ad[0]*ad[2] + mm[4]*ad[1]*ad[2]);
  rapmax = M_MAX(rapmax,h3); 

  bc[0] = c[0] - b[0];
  bc[1] = c[1] - b[1];
  bc[2] = c[2] - b[2];

  bd[0] = d[0] - b[0];
  bd[1] = d[1] - b[1];
  bd[2] = d[2] - b[2];

  cd[0] = d[0] - c[0];
  cd[1] = d[1] - c[1];
  cd[2] = d[2] - c[2];

  h4 =      mm[0]*bd[0]*bd[0] + mm[3]*bd[1]*bd[1] + mm[5]*bd[2]*bd[2] \
     + 2.0*(mm[1]*bd[0]*bd[1] + mm[2]*bd[0]*bd[2] + mm[4]*bd[1]*bd[2]);
  rapmax = M_MAX(rapmax,h4); 
  h5 =      mm[0]*cd[0]*cd[0] + mm[3]*cd[1]*cd[1] + mm[5]*cd[2]*cd[2] \
     + 2.0*(mm[1]*cd[0]*cd[1] + mm[2]*cd[0]*cd[2] + mm[4]*cd[1]*cd[2]);
  rapmax = M_MAX(rapmax,h5); 
  h6 =      mm[0]*bc[0]*bc[0] + mm[3]*bc[1]*bc[1] + mm[5]*bc[2]*bc[2] \
     + 2.0*(mm[1]*bc[0]*bc[1] + mm[2]*bc[0]*bc[2] + mm[4]*bc[1]*bc[2]);
  rapmax = M_MAX(rapmax,h6); 

  /* surfaces */
  dd  = determ(bd,bc,mm);
  air = sqrt(dd);

  dd   = determ(ac,ad,mm);
  air += sqrt(dd);

  dd   = determ(ad,ab,mm);
  air += sqrt(dd);

  dd   = determ(ab,ac,mm);
  air += sqrt(dd);
  
  /* quality */
  num = sqrt(rapmax) * air;
  cal = num / det;
  assert( cal < CALLIM );

  return(cal);
}


double MMG_calte1_iso(pMesh mesh,pSol sol,int iel) {
  pTetra     pt;
  double     cal,abx,aby,abz,acx,acy,acz,adx,ady,adz;
  double     bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz;
  double     h1,h2,h3,h4,h5,h6,rapmax,vol,v1,v2,v3;
  double     s1,s2,s3,s4,dd,d2,num;
  double    *a,*b,*c,*d;
  int        ia,ib,ic,id;

  cal = CALLIM;
  pt = &mesh->tetra[iel];
  if ( !pt->v[0] )  return(cal);

  ia = pt->v[0];
  ib = pt->v[1];
  ic = pt->v[2];
  id = pt->v[3];
  a  = mesh->point[ia].c;
  b  = mesh->point[ib].c;
  c  = mesh->point[ic].c;
  d  = mesh->point[id].c;

  /* volume */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];

  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];
 
  adx = d[0] - a[0];
  ady = d[1] - a[1];
  adz = d[2] - a[2];

  v1  = acy*adz - acz*ady;
  v2  = acz*adx - acx*adz;
  v3  = acx*ady - acy*adx;
  vol = abx * v1 + aby * v2 + abz * v3;
  if ( vol <= 0. )  return(cal);

  /* max edge */
  rapmax = 0.0;
  h1 = abx*abx + aby*aby + abz*abz;
  rapmax = M_MAX(rapmax,h1);
  h2 = acx*acx + acy*acy + acz*acz;
  rapmax = M_MAX(rapmax,h2);
  h3 = adx*adx + ady*ady + adz*adz;
  rapmax = M_MAX(rapmax,h3);

  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];

  bdx = d[0] - b[0];
  bdy = d[1] - b[1];
  bdz = d[2] - b[2];

  cdx = d[0] - c[0];
  cdy = d[1] - c[1];
  cdz = d[2] - c[2];

  h4 = bdx*bdx + bdy*bdy + bdz*bdz;
  rapmax = M_MAX(rapmax,h4);
  h5 = cdx*cdx + cdy*cdy + cdz*cdz;
  rapmax = M_MAX(rapmax,h5);
  h6 = bcx*bcx + bcy*bcy + bcz*bcz;
  rapmax = M_MAX(rapmax,h6);

  /* surface */
  dd = cdy*bdz - cdz*bdy;
  d2 = dd * dd;
  dd = cdz*bdx - cdx*bdz;
  d2 = d2 + dd * dd;
  dd = cdx*bdy - cdy*bdx;
  d2 = d2 + dd * dd;
  s1 = sqrt(d2);
  
  s2 = sqrt(v1*v1 + v2*v2 + v3*v3);
  
  dd = bdy*adz - bdz*ady;
  d2 = dd * dd;
  dd = bdz*adx - bdx*adz;
  d2 = d2 + dd * dd;
  dd = bdx*ady - bdy*adx;
  d2 = d2 + dd * dd;
  s3 = sqrt(d2);

  dd = aby*acz - abz*acy;
  d2 = dd * dd;
  dd = abz*acx - abx*acz;
  d2 = d2 + dd * dd;
  dd = abx*acy - aby*acx;
  d2 = d2 + dd * dd;
  s4 = sqrt(d2);

  /* quality */
  num = sqrt(rapmax) * (s1+s2+s3+s4);
  cal = num / vol;
  assert( cal < CALLIM );

  return(cal);
}

/* compute tetra quality aniso : (ia ib ic id) and (ie ib id ic) */
int MMG_caltet2_ani(pMesh mesh,pSol sol,int iel,int ie,double crit, double * caltab) {
  pTetra     pt;
  double     cal,acx,acy,acz,adx,ady,adz;
  double     bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz,bax,bay,baz;
  double     bex,bey,bez,ecx,ecy,ecz,edx,edy,edz;
  double     h1,h2,h3,h4,h5,h6,det,dete,vol,vole,rap,v1,v2,v3;
  double     h1e,h2e,h3e,h4e,h5e,h6e,num;
  double    *a,*b,*c,*d,*e;
  double     *ma,*mb,*mc,*md,mm[6],mme[6],*me;
  int        j,ia,ib,ic,id,iadr;

  cal       = CALLIM;
  caltab[0] = cal;
  caltab[1] = cal;
  pt  = &mesh->tetra[iel];

  ia = pt->v[0];
  ib = pt->v[1];
  ic = pt->v[2];
  id = pt->v[3];

  /* average metric */
  memset(mm,0,6*sizeof(double));
  iadr = (ia-1)*sol->offset + 1;
  ma   = &sol->met[iadr];
  
  iadr = (ie-1)*sol->offset + 1;
  me   = &sol->met[iadr];
  
  iadr = (ib-1)*sol->offset + 1;
  mb   = &sol->met[iadr];
  
  iadr = (ic-1)*sol->offset + 1;
  mc   = &sol->met[iadr];
  
  iadr = (id-1)*sol->offset + 1;
  md   = &sol->met[iadr];
  for (j=0; j<6; j++) {
    mm[j]  = 0.25 * (ma[j]+mb[j]+mc[j]+md[j]);
    mme[j] = 0.25 * (me[j]+mb[j]+mc[j]+md[j]);
  }

  a = mesh->point[ia].c;
  b = mesh->point[ib].c;
  c = mesh->point[ic].c;
  d = mesh->point[id].c;
  e = mesh->point[ie].c;
  
  /* volume */
  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];

  bax = a[0] - b[0];
  bay = a[1] - b[1];
  baz = a[2] - b[2];
  
  bex = e[0] - b[0];
  bey = e[1] - b[1];
  bez = e[2] - b[2];
  
  bdx = d[0] - b[0];
  bdy = d[1] - b[1];
  bdz = d[2] - b[2];


  v1  = bdy*bcz - bdz*bcy;
  v2  = bdz*bcx - bdx*bcz;
  v3  = bdx*bcy - bdy*bcx;
  
  vol = bax * v1 + bay * v2 + baz * v3;
  if ( vol <= 0. )  return(0);
  
  vole = -bex * v1 - bey * v2 - bez * v3;
  if ( vole <= 0. ) return(0);
  
  det = mm[0] * ( mm[3]*mm[5] - mm[4]*mm[4]) \
      - mm[1] * ( mm[1]*mm[5] - mm[2]*mm[4]) \
      + mm[2] * ( mm[1]*mm[4] - mm[2]*mm[3]);
  if ( det < EPSOK )  return(0);
  
  dete = mme[0] * ( mme[3]*mme[5] - mme[4]*mme[4]) \
       - mme[1] * ( mme[1]*mme[5] - mme[2]*mme[4]) \
       + mme[2] * ( mme[1]*mme[4] - mme[2]*mme[3]);  
  if ( dete < EPSOK )  return(0);
  
  det = sqrt(det) * vol;

  /* edge lengths */
  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];
 
  adx = d[0] - a[0];
  ady = d[1] - a[1];
  adz = d[2] - a[2];
  
  h1 =      mm[0]*bax*bax + mm[3]*bay*bay + mm[5]*baz*baz \
     + 2.0*(mm[1]*bax*bay + mm[2]*bax*baz + mm[4]*bay*baz);
  h2 =      mm[0]*acx*acx + mm[3]*acy*acy + mm[5]*acz*acz \
     + 2.0*(mm[1]*acx*acy + mm[2]*acx*acz + mm[4]*acy*acz);
  h3 =      mm[0]*adx*adx + mm[3]*ady*ady + mm[5]*adz*adz \
     + 2.0*(mm[1]*adx*ady + mm[2]*adx*adz + mm[4]*ady*adz);
     
  cdx = d[0] - c[0];
  cdy = d[1] - c[1];
  cdz = d[2] - c[2];

  h4 =      mm[0]*bdx*bdx + mm[3]*bdy*bdy + mm[5]*bdz*bdz \
     + 2.0*(mm[1]*bdx*bdy + mm[2]*bdx*bdz + mm[4]*bdy*bdz);
  h5 =      mm[0]*cdx*cdx + mm[3]*cdy*cdy + mm[5]*cdz*cdz \
     + 2.0*(mm[1]*cdx*cdy + mm[2]*cdx*cdz + mm[4]*cdy*cdz);
  h6 =      mm[0]*bcx*bcx + mm[3]*bcy*bcy + mm[5]*bcz*bcz \
     + 2.0*(mm[1]*bcx*bcy + mm[2]*bcx*bcz + mm[4]*bcy*bcz);

  /* quality */
  rap = h1 + h2 + h3 + h4 + h5 + h6;
  num = sqrt(rap) * rap;
  caltab[0] = num / det;  
  if ( caltab[0] > crit )  return(0);

  dete = sqrt(dete) * vole;
  ecx = c[0] - e[0];
  ecy = c[1] - e[1];
  ecz = c[2] - e[2];
 
  edx = d[0] - e[0];
  edy = d[1] - e[1];
  edz = d[2] - e[2];
 
  h1e =      mme[0]*bex*bex + mme[3]*bey*bey + mme[5]*bez*bez \
      + 2.0*(mme[1]*bex*bey + mme[2]*bex*bez + mme[4]*bey*bez);
  h2e =      mme[0]*ecx*ecx + mme[3]*ecy*ecy + mme[5]*ecz*ecz \
      + 2.0*(mme[1]*ecx*ecy + mme[2]*ecx*ecz + mme[4]*ecy*ecz);
  h3e =      mme[0]*edx*edx + mme[3]*edy*edy + mme[5]*edz*edz \
      + 2.0*(mme[1]*edx*edy + mme[2]*edx*edz + mme[4]*edy*edz);   

  h4e =      mme[0]*bdx*bdx + mme[3]*bdy*bdy + mme[5]*bdz*bdz \
      + 2.0*(mme[1]*bdx*bdy + mme[2]*bdx*bdz + mme[4]*bdy*bdz);
  h5e =      mme[0]*cdx*cdx + mme[3]*cdy*cdy + mme[5]*cdz*cdz \
      + 2.0*(mme[1]*cdx*cdy + mme[2]*cdx*cdz + mme[4]*cdy*cdz);
  h6e =      mme[0]*bcx*bcx + mme[3]*bcy*bcy + mme[5]*bcz*bcz \
      + 2.0*(mme[1]*bcx*bcy + mme[2]*bcx*bcz + mme[4]*bcy*bcz);   

  rap = h1e + h2e + h3e + h4e + h5e + h6e;
  num = sqrt(rap) * rap;
  caltab[1] = ( sqrt(rap) * rap ) / dete;
  
  if ( caltab[1] > crit )  return(0);
  
  return(1);
}


/* compute tetra quality iso : (ia ib ic id) and (ie ib id ic) */
int MMG_caltet2_iso(pMesh mesh,pSol sol,int iel,int ie,double crit,double * caltab) {
  pTetra     pt;
  double     cal,acx,acy,acz,adx,ady,adz;
  double     bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz,bax,bay,baz;
  double     bex,bey,bez,ecx,ecy,ecz,edx,edy,edz;
  double     h1,h2,h3,h4,h5,h6,vol,vole,rap,v1,v2,v3;
  double     h1e,h2e,h3e,num;
  double    *a,*b,*c,*d,*e;
  int        ia,ib,ic,id;

  cal       = CALLIM;
  caltab[0] = cal;
  caltab[1] = cal;
  pt  = &mesh->tetra[iel];  
  ia = pt->v[0];
  ib = pt->v[1];
  ic = pt->v[2];
  id = pt->v[3];

  a = mesh->point[ia].c;
  b = mesh->point[ib].c;
  c = mesh->point[ic].c;
  d = mesh->point[id].c;
  e = mesh->point[ie].c;
  
  /* volume */
  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];

  bax = a[0] - b[0];
  bay = a[1] - b[1];
  baz = a[2] - b[2];
  
  bex = e[0] - b[0];
  bey = e[1] - b[1];
  bez = e[2] - b[2];
  
  bdx = d[0] - b[0];
  bdy = d[1] - b[1];
  bdz = d[2] - b[2];


  v1  = bdy*bcz - bdz*bcy;
  v2  = bdz*bcx - bdx*bcz;
  v3  = bdx*bcy - bdy*bcx;
  vol = bax * v1 + bay * v2 + baz * v3;
  if ( vol <= 0. )  return(0);
  
  vole = -bex * v1 - bey * v2 - bez * v3;
  if ( vole <= 0. ) return(0);
/////////////////////    
/*
caltab[0] = MMG_caltetrao(mesh,sol,iel);   
if ( caltab[0] > crit )  {
	return(0);
}
pt = &mesh->tetra[0];
pt->v[0] = ie;
pt->v[1] = ib;
pt->v[2] = id;
pt->v[3] = ic;       
caltab[1] = MMG_caltetrao(mesh,sol,0);   
if ( caltab[1] > crit )  {
	return(0);
}
return(1);      */
////////////////////
  
  /* edge lengths */
  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];
 
  adx = d[0] - a[0];
  ady = d[1] - a[1];
  adz = d[2] - a[2];
  
  
  h1 = bax*bax + bay*bay + baz*baz;
  h2 = acx*acx + acy*acy + acz*acz;
  h3 = adx*adx + ady*ady + adz*adz;
     

  cdx = d[0] - c[0];
  cdy = d[1] - c[1];
  cdz = d[2] - c[2];

  h4 = bdx*bdx + bdy*bdy + bdz*bdz;
  h5 = cdx*cdx + cdy*cdy + cdz*cdz;                                             
  h6 = bcx*bcx + bcy*bcy + bcz*bcz;
     
  /* quality */
  rap = h1 + h2 + h3 + h4 + h5 + h6;
  num = sqrt(rap) * rap;
  caltab[0] = num / vol;
  
  if ( caltab[0] > crit )  return(0);
  
  ecx = c[0] - e[0];
  ecy = c[1] - e[1];
  ecz = c[2] - e[2];
 
  edx = d[0] - e[0];
  edy = d[1] - e[1];
  edz = d[2] - e[2];

  h1e = bex*bex + bey*bey + bez*bez;
  h2e = ecx*ecx + ecy*ecy + ecz*ecz;
  h3e = edx*edx + edy*edy + edz*edz;   
 
  rap = h1e + h2e + h3e + h4 + h5 + h6;
  num = sqrt(rap) * rap;
  caltab[1] = num / vole;
  
  if ( caltab[1] > crit )  return(0);

  return(1);
}

/* compute tetra quality : 60*max(1/lmin,lmax) + iare : (ia ib ic id) and (ie ib id ic) */
int MMG_caltet2long_ani(pMesh mesh,pSol sol,int iel,int ie,double crit, double * caltab) {
  pTetra     pt;
  double     cal,acx,acy,acz,adx,ady,adz;
  double     bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz,bax,bay,baz;
  double     bex,bey,bez,ecx,ecy,ecz,edx,edy,edz;
  double     h1,h2,h3,h4,h5,h6,det,dete,vol,vole,v1,v2,v3;
  double     h1e,h2e,h3e,h4e,h5e,h6e;
  double    *a,*b,*c,*d,*e;
  double     *ma,*mb,*mc,*md,mm[6],mme[6],*me; 
  double     lmin,lmax;
  int        j,ia,ib,ic,id,iadr,iarmin,iarmax;

  cal       = CALLIM;
  caltab[0] = cal;
  caltab[1] = cal;
  pt  = &mesh->tetra[iel];

  ia = pt->v[0];
  ib = pt->v[1];
  ic = pt->v[2];
  id = pt->v[3];

  /* average metric */
  memset(mm,0,6*sizeof(double));
  iadr = (ia-1)*sol->offset + 1;
  ma   = &sol->met[iadr];
  
  iadr = (ie-1)*sol->offset + 1;
  me   = &sol->met[iadr];
  
  iadr = (ib-1)*sol->offset + 1;
  mb   = &sol->met[iadr];
  
  iadr = (ic-1)*sol->offset + 1;
  mc   = &sol->met[iadr];
  
  iadr = (id-1)*sol->offset + 1;
  md   = &sol->met[iadr];
  for (j=0; j<6; j++) {
    mm[j]  = 0.25 * (ma[j]+mb[j]+mc[j]+md[j]);
    mme[j] = 0.25 * (me[j]+mb[j]+mc[j]+md[j]);
  }

  a = mesh->point[ia].c;
  b = mesh->point[ib].c;
  c = mesh->point[ic].c;
  d = mesh->point[id].c;
  e = mesh->point[ie].c;
  
  /* volume */
  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];

  bax = a[0] - b[0];
  bay = a[1] - b[1];
  baz = a[2] - b[2];
  
  bex = e[0] - b[0];
  bey = e[1] - b[1];
  bez = e[2] - b[2];
  
  bdx = d[0] - b[0];
  bdy = d[1] - b[1];
  bdz = d[2] - b[2];


  v1  = bdy*bcz - bdz*bcy;
  v2  = bdz*bcx - bdx*bcz;
  v3  = bdx*bcy - bdy*bcx;
  
  vol = bax * v1 + bay * v2 + baz * v3;
  if ( vol <= 0. )  return(0);
  
  vole = -bex * v1 - bey * v2 - bez * v3;
  if ( vole <= 0. ) return(0);
  
  det = mm[0] * ( mm[3]*mm[5] - mm[4]*mm[4]) \
      - mm[1] * ( mm[1]*mm[5] - mm[2]*mm[4]) \
      + mm[2] * ( mm[1]*mm[4] - mm[2]*mm[3]);
  if ( det < EPSOK )  return(0);
  
  dete = mme[0] * ( mme[3]*mme[5] - mme[4]*mme[4]) \
       - mme[1] * ( mme[1]*mme[5] - mme[2]*mme[4]) \
       + mme[2] * ( mme[1]*mme[4] - mme[2]*mme[3]);  
  if ( dete < EPSOK )  return(0);
  
  det = sqrt(det) * vol;

  /* edge lengths */
  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];
 
  adx = d[0] - a[0];
  ady = d[1] - a[1];
  adz = d[2] - a[2];
  
  h1 =      mm[0]*bax*bax + mm[3]*bay*bay + mm[5]*baz*baz \
     + 2.0*(mm[1]*bax*bay + mm[2]*bax*baz + mm[4]*bay*baz);
  h2 =      mm[0]*acx*acx + mm[3]*acy*acy + mm[5]*acz*acz \
     + 2.0*(mm[1]*acx*acy + mm[2]*acx*acz + mm[4]*acy*acz);
  h3 =      mm[0]*adx*adx + mm[3]*ady*ady + mm[5]*adz*adz \
     + 2.0*(mm[1]*adx*ady + mm[2]*adx*adz + mm[4]*ady*adz);
     
  cdx = d[0] - c[0];
  cdy = d[1] - c[1];
  cdz = d[2] - c[2];

  h4 =      mm[0]*bdx*bdx + mm[3]*bdy*bdy + mm[5]*bdz*bdz \
     + 2.0*(mm[1]*bdx*bdy + mm[2]*bdx*bdz + mm[4]*bdy*bdz);
  h5 =      mm[0]*cdx*cdx + mm[3]*cdy*cdy + mm[5]*cdz*cdz \
     + 2.0*(mm[1]*cdx*cdy + mm[2]*cdx*cdz + mm[4]*cdy*cdz);
  h6 =      mm[0]*bcx*bcx + mm[3]*bcy*bcy + mm[5]*bcz*bcz \
     + 2.0*(mm[1]*bcx*bcy + mm[2]*bcx*bcz + mm[4]*bcy*bcz);

  /* quality */ 
  if( h1 < h2 ) {
    lmin = h1;
    iarmin = 0;
    lmax = h2;
    iarmax = 1;
  } else {
    lmin = h2;
    iarmin = 1;    
    lmax = h1;
    iarmax = 0;
  }
  
  if( h3 < lmin) {
    lmin = h3;
    iarmin = 2;
  }
  if( h3 > lmax) {
    lmax = h3;
    iarmax = 2;
  }
  if( h4 < lmin) {
    lmin = h4;
    iarmin = 3;
  }
  if( h4 > lmax) {
    lmax = h4;
    iarmax = 3;
  }
  if( h5 < lmin) {
    lmin = h5;
    iarmin = 4;
  }
  if( h5 > lmax) {
    lmax = h5;
    iarmax = 4;
  }
  if( h6 < lmin) {
    lmin = h6;
    iarmin = 5;
  }
  if( h6 > lmax) {
    lmax = h6;
    iarmax = 5;
  }                 
  
  lmin = sqrt(lmin);
  lmax = sqrt(lmax);     
  lmin = (lmin > 1.) ? lmin : 1./lmin;
  lmax = (lmax > 1.) ? lmax : 1./lmax;
  if ( lmin > lmax) {
     caltab[0] = 60*lmin + iarmin;
  } else {
     caltab[0] = 60*lmax + iarmax;
  } 
  
  /*rap = h1 + h2 + h3 + h4 + h5 + h6;
  num = sqrt(rap) * rap;
  caltab[0] = num / det;  
  */
  if ( caltab[0] > crit )  return(0);

  dete = sqrt(dete) * vole;
  ecx = c[0] - e[0];
  ecy = c[1] - e[1];
  ecz = c[2] - e[2];
 
  edx = d[0] - e[0];
  edy = d[1] - e[1];
  edz = d[2] - e[2];
 
  h1e =      mme[0]*bex*bex + mme[3]*bey*bey + mme[5]*bez*bez \
      + 2.0*(mme[1]*bex*bey + mme[2]*bex*bez + mme[4]*bey*bez);
  h2e =      mme[0]*ecx*ecx + mme[3]*ecy*ecy + mme[5]*ecz*ecz \
      + 2.0*(mme[1]*ecx*ecy + mme[2]*ecx*ecz + mme[4]*ecy*ecz);
  h3e =      mme[0]*edx*edx + mme[3]*edy*edy + mme[5]*edz*edz \
      + 2.0*(mme[1]*edx*edy + mme[2]*edx*edz + mme[4]*edy*edz);   

  h4e =      mme[0]*bdx*bdx + mme[3]*bdy*bdy + mme[5]*bdz*bdz \
      + 2.0*(mme[1]*bdx*bdy + mme[2]*bdx*bdz + mme[4]*bdy*bdz);
  h5e =      mme[0]*cdx*cdx + mme[3]*cdy*cdy + mme[5]*cdz*cdz \
      + 2.0*(mme[1]*cdx*cdy + mme[2]*cdx*cdz + mme[4]*cdy*cdz);
  h6e =      mme[0]*bcx*bcx + mme[3]*bcy*bcy + mme[5]*bcz*bcz \
      + 2.0*(mme[1]*bcx*bcy + mme[2]*bcx*bcz + mme[4]*bcy*bcz);   

  if( h1e < h2e ) {
    lmin = h1e;
    iarmin = 0;
    lmax = h2e;
    iarmax = 1;
  } else {
    lmin = h2e;
    iarmin = 1;    
    lmax = h1e;
    iarmax = 0;
  }
  
  if( h3e < lmin) {
    lmin = h3e;
    iarmin = 2;
  }
  if( h3e > lmax) {
    lmax = h3;
    iarmax = 2;
  }
  if( h4 < lmin) {
    lmin = h4;
    iarmin = 3;
  }
  if( h4 > lmax) {
    lmax = h4;
    iarmax = 3;
  }
  if( h5 < lmin) {
    lmin = h5;
    iarmin = 4;
  }
  if( h5 > lmax) {
    lmax = h5;
    iarmax = 4;
  }
  if( h6 < lmin) {
    lmin = h6;
    iarmin = 5;
  }
  if( h6 > lmax) {
    lmax = h6;
    iarmax = 5;
  }

  
  lmin = sqrt(lmin);
  lmax = sqrt(lmax);
  lmin = (lmin > 1.) ? lmin : 1./lmin;
  lmax = (lmax > 1.) ? lmax : 1./lmax;
  if ( lmin > lmax) {
     caltab[1] = 60*lmin + iarmin;
  } else {
     caltab[1] = 60*lmax + iarmax;
  } 
  /*rap = h1e + h2e + h3e + h4e + h5e + h6e;
  num = sqrt(rap) * rap;
  caltab[1] = ( sqrt(rap) * rap ) / dete;
  */
  if ( caltab[1] > crit )  return(0);
  
  return(1);
}


/* compute tetra quality iso : (ia ib ic id) and (ie ib id ic) */
int MMG_caltet2long_iso(pMesh mesh,pSol sol,int iel,int ie,double crit,double * caltab) {
  pTetra     pt;
  double     cal,acx,acy,acz,adx,ady,adz;
  double     bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz,bax,bay,baz;
  double     bex,bey,bez,ecx,ecy,ecz,edx,edy,edz;
  double     h1,h2,h3,h4,h5,h6,vol,vole,rap,v1,v2,v3;
  double     h1e,h2e,h3e,num;
  double    *a,*b,*c,*d,*e; 
  double     lmin,lmax;
  int        ia,ib,ic,id,iarmin,iarmax;
  cal       = CALLIM;
  caltab[0] = cal;
  caltab[1] = cal;
  pt  = &mesh->tetra[iel];  
  ia = pt->v[0];
  ib = pt->v[1];
  ic = pt->v[2];
  id = pt->v[3];

  a = mesh->point[ia].c;
  b = mesh->point[ib].c;
  c = mesh->point[ic].c;
  d = mesh->point[id].c;
  e = mesh->point[ie].c;
  
  /* volume */
  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];

  bax = a[0] - b[0];
  bay = a[1] - b[1];
  baz = a[2] - b[2];
  
  bex = e[0] - b[0];
  bey = e[1] - b[1];
  bez = e[2] - b[2];
  
  bdx = d[0] - b[0];
  bdy = d[1] - b[1];
  bdz = d[2] - b[2];


  v1  = bdy*bcz - bdz*bcy;
  v2  = bdz*bcx - bdx*bcz;
  v3  = bdx*bcy - bdy*bcx;
  vol = bax * v1 + bay * v2 + baz * v3;
  if ( vol <= 0. )  return(0);
  
  vole = -bex * v1 - bey * v2 - bez * v3;
  if ( vole <= 0. ) return(0);
/////////////////////    
/*
caltab[0] = MMG_caltetrao(mesh,sol,iel);   
if ( caltab[0] > crit )  {
	return(0);
}
pt = &mesh->tetra[0];
pt->v[0] = ie;
pt->v[1] = ib;
pt->v[2] = id;
pt->v[3] = ic;       
caltab[1] = MMG_caltetrao(mesh,sol,0);   
if ( caltab[1] > crit )  {
	return(0);
}
return(1);      */
////////////////////
  
  /* edge lengths */
  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];
 
  adx = d[0] - a[0];
  ady = d[1] - a[1];
  adz = d[2] - a[2];
  
  
  h1 = bax*bax + bay*bay + baz*baz;
  h2 = acx*acx + acy*acy + acz*acz;
  h3 = adx*adx + ady*ady + adz*adz;
     

  cdx = d[0] - c[0];
  cdy = d[1] - c[1];
  cdz = d[2] - c[2];

  h4 = bdx*bdx + bdy*bdy + bdz*bdz;
  h5 = cdx*cdx + cdy*cdy + cdz*cdz;                                             
  h6 = bcx*bcx + bcy*bcy + bcz*bcz;
     
  /* quality */
  if( h1 < h2 ) {
    lmin = h1;
    iarmin = 0;
    lmax = h2;
    iarmax = 1;
  } else {
    lmin = h2;
    iarmin = 1;    
    lmax = h1;
    iarmax = 0;
  }
  
  if( h3 < lmin) {
    lmin = h3;
    iarmin = 2;
  }
  if( h3 > lmax) {
    lmax = h3;
    iarmax = 2;
  }
  if( h4 < lmin) {
    lmin = h4;
    iarmin = 3;
  }
  if( h4 > lmax) {
    lmax = h4;
    iarmax = 3;
  }
  if( h5 < lmin) {
    lmin = h5;
    iarmin = 4;
  }
  if( h5 > lmax) {
    lmax = h5;
    iarmax = 4;
  }
  if( h6 < lmin) {
    lmin = h6;
    iarmin = 5;
  }
  if( h6 > lmax) {
    lmax = h6;
    iarmax = 5;
  }
  lmin = sqrt(lmin);
  lmax = sqrt(lmax);   
  lmin = (lmin > 1.) ? lmin : 1./lmin;
  lmax = (lmax > 1.) ? lmax : 1./lmax; 
  //printf("long %e %e %e %e %e %e\n",sqrt(h1),sqrt(h2),sqrt(h3),sqrt(h4),sqrt(h5),sqrt(h6));
  //printf("--lmin %e %e\n",lmin,lmax);
  if ( lmin > lmax) {
     caltab[0] = 60*lmin + iarmin;
  } else {
     caltab[0] = 60*lmax + iarmax;
  }
  rap = h1 + h2 + h3 + h4 + h5 + h6;
  num = sqrt(rap) * rap;
  cal = num / vol;
              
  //printf("cal %e ? %e\n",caltab[0],crit);
  if ( cal > crit )  return(0);   
  
  ecx = c[0] - e[0];
  ecy = c[1] - e[1];
  ecz = c[2] - e[2];
 
  edx = d[0] - e[0];
  edy = d[1] - e[1];
  edz = d[2] - e[2];

  h1e = bex*bex + bey*bey + bez*bez;
  h2e = ecx*ecx + ecy*ecy + ecz*ecz;
  h3e = edx*edx + edy*edy + edz*edz;   
 
  if( h1e < h2e ) {
    lmin = h1e;
    iarmin = 0;
    lmax = h2e;
    iarmax = 1;
  } else {
    lmin = h2e;
    iarmin = 1;    
    lmax = h1e;
    iarmax = 0;
  }
  
  if( h3e < lmin) {
    lmin = h3e;
    iarmin = 2;
  }
  if( h3e > lmax) {
    lmax = h3;
    iarmax = 2;
  }
  if( h4 < lmin) {
    lmin = h4;
    iarmin = 3;
  }
  if( h4 > lmax) {
    lmax = h4;
    iarmax = 3;
  }
  if( h5 < lmin) {
    lmin = h5;
    iarmin = 4;
  }
  if( h5 > lmax) {
    lmax = h5;
    iarmax = 4;
  }
  if( h6 < lmin) {
    lmin = h6;
    iarmin = 5;
  }
  if( h6 > lmax) {
    lmax = h6;
    iarmax = 5;
  }
  
  lmin = sqrt(lmin);
  lmax = sqrt(lmax);     
  lmin = (lmin > 1.) ? lmin : 1./lmin;
  lmax = (lmax > 1.) ? lmax : 1./lmax;   
  //printf("lmin %e %e\n",lmin,lmax);
  if ( lmin > lmax) {
     caltab[1] = 60*lmin + iarmin;
  } else {
     caltab[1] = 60*lmax + iarmax;
  }
  rap = h1e + h2e + h3e + h4 + h5 + h6;
  num = sqrt(rap) * rap;
  cal = num / vole;
  
  //printf("cal %e ? %e\n",caltab[1],crit);
  if ( cal > crit )  return(0);
//  puts("je passe");exit(0);

  return(1);
}

double MMG_calte3_ani(pMesh mesh,pSol sol,int iel) {
  pTetra     pt;
  double     cal,ab[3],ac[3],ad[3],bc[3],bd[3],cd[3];
  double     vol,det,v1,v2,v3,air,dd;
  double    *a,*b,*c,*d;
  double     *ma,*mb,*mc,*md,mm[6];;
  int        j,ia,ib,ic,id,iadr;
  static double id3[6] ={1.0,0.,0.,1.,0.,1.};

  pt = &mesh->tetra[iel];
  if ( !pt->v[0] )  return(0.0);
  cal = CALLIM;

  ia = pt->v[0];
  ib = pt->v[1];
  ic = pt->v[2];
  id = pt->v[3];

  /* average metric */
  memset(mm,0,6*sizeof(double));
  iadr = (ia-1)*sol->offset + 1;
  ma   = &sol->met[iadr];
  iadr = (ib-1)*sol->offset + 1;
  mb   = &sol->met[iadr];
  iadr = (ic-1)*sol->offset + 1;
  mc   = &sol->met[iadr];
  iadr = (id-1)*sol->offset + 1;
  md   = &sol->met[iadr];
  for (j=0; j<6; j++)
    mm[j] = 0.25 * (ma[j]+mb[j]+mc[j]+md[j]);

  a  = mesh->point[ia].c;
  b  = mesh->point[ib].c;
  c  = mesh->point[ic].c;
  d  = mesh->point[id].c;

  /* volume */
  ab[0] = b[0] - a[0];
  ab[1] = b[1] - a[1];
  ab[2] = b[2] - a[2];

  ac[0] = c[0] - a[0];
  ac[1] = c[1] - a[1];
  ac[2] = c[2] - a[2];
  
  ad[0] = d[0] - a[0];
  ad[1] = d[1] - a[1];
  ad[2] = d[2] - a[2];

  v1  = ac[1]*ad[2] - ac[2]*ad[1];
  v2  = ac[2]*ad[0] - ac[0]*ad[2];
  v3  = ac[0]*ad[1] - ac[1]*ad[0];
  vol = ab[0] * v1 + ab[1] * v2 + ab[2] * v3;
  if ( vol <= 0. )  return(cal);

  det = mm[0] * ( mm[3]*mm[5] - mm[4]*mm[4]) \
      - mm[1] * ( mm[1]*mm[5] - mm[2]*mm[4]) \
      + mm[2] * ( mm[1]*mm[4] - mm[2]*mm[3]);
  if ( det < EPSOK )  return(cal);
  det = sqrt(det) * vol;

  bc[0] = c[0] - b[0];
  bc[1] = c[1] - b[1];
  bc[2] = c[2] - b[2];

  bd[0] = d[0] - b[0];
  bd[1] = d[1] - b[1];
  bd[2] = d[2] - b[2];

  cd[0] = d[0] - c[0];
  cd[1] = d[1] - c[1];
  cd[2] = d[2] - c[2];

  /* surfaces */
  memcpy(mm,id3,6*sizeof(double));
  dd  = determ(bd,bc,mm);
  air = dd;
if ( iel ==17460 )  printf("aire1 %E\n",sqrt(dd));

  dd   = determ(ac,ad,mm);
if ( iel ==17460 )  printf("aire2 %E\n",sqrt(dd));
  air  = M_MAX(air,dd);

  dd   = determ(ad,ab,mm);
if ( iel ==17460 )  printf("aire3 %E\n",sqrt(dd));
  air  = M_MAX(air,dd);

  dd   = determ(ab,ac,mm);
if ( iel ==17460 )  printf("aire4 %E\n",sqrt(dd));
  air  = M_MAX(air,dd);
  
  /* quality */
  if ( iel == 20496 ) {
    printf("vol %E  \n",vol);
    printf("a %d: %f %f %f\n",ia,a[0],a[1],a[2]);
    printf("b %d: %f %f %f\n",ib,b[0],b[1],b[2]);
    printf("c %d: %f %f %f\n",ic,c[0],c[1],c[2]);
    printf("d %d: %f %f %f\n",id,d[0],d[1],d[2]);
  }
  cal = 3.0*vol / sqrt(air);

  return(cal);
}

/* compute tetra quality to be the cubic mean ration : 15552 * V^2/(somme l^2)^3 */
double MMG_caltetcubic(pMesh mesh,pSol sol,int iel) {
  pTetra     pt;
  double     cal,abx,aby,abz,acx,acy,acz,adx,ady,adz;
  double     bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz;
  double     h1,h2,h3,h4,h5,h6,vol,rap,v1,v2,v3;
  double    *a,*b,*c,*d;
  int        ia,ib,ic,id;

  cal = 0.;
  pt  = &mesh->tetra[iel];
  if ( !pt->v[0] )  return(cal);

  ia = pt->v[0];
  ib = pt->v[1];
  ic = pt->v[2];
  id = pt->v[3];

  a = mesh->point[ia].c;
  b = mesh->point[ib].c;
  c = mesh->point[ic].c;
  d = mesh->point[id].c;

  /* volume */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];

  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];
  
  adx = d[0] - a[0];
  ady = d[1] - a[1];
  adz = d[2] - a[2];

  v1  = acy*adz - acz*ady;
  v2  = acz*adx - acx*adz;
  v3  = acx*ady - acy*adx;
  vol = abx * v1 + aby * v2 + abz * v3;
  if ( vol <= 0. )  return(cal); 
  vol /= 6.;

  /*  edge */
  h1 = abx*abx + aby*aby + abz*abz;
  h2 = acx*acx + acy*acy + acz*acz;
  h3 = adx*adx + ady*ady + adz*adz;

  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];

  bdx = d[0] - b[0];
  bdy = d[1] - b[1];
  bdz = d[2] - b[2];

  cdx = d[0] - c[0];
  cdy = d[1] - c[1];
  cdz = d[2] - c[2];

  h4 = bdx*bdx + bdy*bdy + bdz*bdz;
  h5 = cdx*cdx + cdy*cdy + cdz*cdz;
  h6 = bcx*bcx + bcy*bcy + bcz*bcz;

  /* quality */
  rap = h1 + h2 + h3 + h4 + h5 + h6;
  cal = vol*vol;
  cal *= 15552.; 
  cal = exp(0.3333333333*log(cal));
  //cal = pow(cal,0.33333333333);
  cal = cal / rap; 

  return(cal);
}

/* compute tetra quality : 60*max(1/lmin,lmax) + iare*/
double MMG_callong(pMesh mesh,pSol sol,int iel) {
  pTetra     pt;
  double    *ca,*cb,len;
  double     *ma,*mb,lmin,lmax,cal; 
  double     ux,uy,uz,dd1,sa,dd2,sb,dd,rap;
  int        i,ia,ib,iadr,iarmin,iarmax,ipa,ipb;

  cal = CALLIM;
  pt  = &mesh->tetra[iel];
  if ( !pt->v[0] )  return(cal);
  
  lmax   = 0.;
  lmin   = CALLIM;
  iarmin = 0;
  iarmax = 0;
  for (i=0; i<6; i++) {
  
    /* edge length */
    ia  = MMG_iare[i][0];
    ib  = MMG_iare[i][1];
    ipa = pt->v[ia];
    ipb = pt->v[ib];

    ca  = &mesh->point[ipa].c[0];
    cb  = &mesh->point[ipb].c[0];
     
    iadr = (ipa-1)*sol->offset + 1;
    ma  = &sol->met[iadr];

    iadr = (ipb-1)*sol->offset + 1;
    mb  = &sol->met[iadr];        
    
    if(sol->offset==6) {
      ux = cb[0] - ca[0];
      uy = cb[1] - ca[1];
      uz = cb[2] - ca[2];

      dd1 =      ma[0]*ux*ux + ma[3]*uy*uy + ma[5]*uz*uz \
          + 2.0*(ma[1]*ux*uy + ma[2]*ux*uz + ma[4]*uy*uz);
      if ( dd1 <= 0.0 )  dd1 = 0.0;

      dd2 =     mb[0]*ux*ux + mb[3]*uy*uy + mb[5]*uz*uz \
          + 2.0*(mb[1]*ux*uy +mb[2]*ux*uz + mb[4]*uy*uz);
      if ( dd2 <= 0.0 )  dd2 = 0.0;

      /*longueur approchee*/   
      /*precision a 3.5 10e-3 pres*/
      if(fabs(dd1-dd2) < 0.05 ) {
        //printf("bonne precision %e \n",sqrt(0.5*(dd1+dd2)) - (sqrt(dd1)+sqrt(dd2)+4.0*sqrt(0.5*(dd1+dd2))) / 6.0 );
        len = sqrt(0.5*(dd1+dd2));
      } else{
        len = (sqrt(dd1)+sqrt(dd2)+4.0*sqrt(0.5*(dd1+dd2))) / 6.0;
      }
    } else {
      sa   = *ma;
      sb   = *mb;

      ux = cb[0] - ca[0];
      uy = cb[1] - ca[1];
      uz = cb[2] - ca[2];
      dd = sqrt(ux*ux + uy*uy + uz*uz);

      rap = (sb - sa) / sa;
      if ( fabs(rap) < EPS1 )
        /*len = dd * (2.0-EPS1) / (2.0*sa);*/
        len = dd / sa;
      else
        /*len = max(dd/sa,dd/sb);*/
        len = dd * (1.0/sa + 1.0/sb + 8.0 / (sa+sb)) / 6.0;
      
    }    
    if ( len < lmin ) {
      lmin = len;
      iarmin = i;
    } 
    if ( len > lmax ) {
      lmax = len;
      iarmax = i;
    }
  }  
  lmin = (lmin > 1.) ? lmin : 1./lmin;
  lmax = (lmax > 1.) ? lmax : 1./lmax;
  if ( lmin > lmax) {
    cal = 60*lmin + iarmin;
  } else {
    cal = 60*lmax + iarmax;
  } 
  return(cal);
}


/* compute tetra quality iso : (ia ib ic id) and (ie ib id ic) */
int MMG_caltet2_LES(pMesh mesh,pSol sol,int iel,int ie,double crit,double * caltab) {
  pTetra     pt,pt0;
  double     cal,acx,acy,acz,adx,ady,adz;
  double     bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz,bax,bay,baz;
  double     bex,bey,bez,ecx,ecy,ecz,edx,edy,edz;
  double     h1,h2,h3,h4,h5,h6,vol,vole,rap,v1,v2,v3;
  double     h1e,h2e,h3e,num;
  double    *a,*b,*c,*d,*e;
  int        ia,ib,ic,id;

  cal       = 1;
  caltab[0] = cal;
  caltab[1] = cal;
  pt  = &mesh->tetra[iel];  
  ia = pt->v[0];
  ib = pt->v[1];
  ic = pt->v[2];
  id = pt->v[3];

  pt0 = &mesh->tetra[0];
  pt->v[0] = ia;
  pt->v[1] = ib;
  pt->v[2] = ic;
  pt->v[3] = id;

  caltab[0] = MMG_caltet_LES(mesh,sol,0);
  
  if ( caltab[0] > crit )  return(0);
  //sol->met[0] = sol->met[iel];
  //caltab[0] = MMG_caltet(mesh,sol,0);

  pt->v[0] = ie;
  pt->v[1] = ib;
  pt->v[2] = id;
  pt->v[3] = ic;

  caltab[1] = MMG_caltet_LES(mesh,sol,0);
  
  if ( caltab[1] > crit )  return(0);
  // caltab[1] = MMG_caltet(mesh,sol,0);
  return(1);
}
