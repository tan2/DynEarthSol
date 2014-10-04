#include "mesh.h"

#define  HQCOEF    0.45
#define  HCRIT     0.8

#define QDEGRADBDRY     1e8//3075. /* *ALPHAD ~  150*/
#define EPSDEV          0.1
#define ALPHATRI ((double)sqrt(3.) / 6.)



/* find all points in the "boundary ball" of P
   in:  listtet  : list of tets sharing P
   list     : dynamic list structure (allocated outside)
   out: list  : list of points */
int MMG_boulebdry(pMesh mesh,int lontet,pList listtet,pList list) {
  pTetra        pt,ptk;
  pPoint        ppt;
  int           *adja,adj,i,j,indp,iel,iadr,base,ilist,nump,ip,l,kk,nk;
  int           depart,fdep,i0,i1,i2,kfa,kel,jel,nobj,nlast;
  int ik;
  i = kfa = nk = kk = -1;

  pt   = &mesh->tetra[listtet->tetra[1]/4];
  if ( !pt->v[0] )  return(0);
  ip   = listtet->tetra[1]%4;
  nump = pt->v[ip];
  ppt  = &mesh->point[nump];
  if (  ppt->tag & M_UNUSED )  return(0);
  if( !(ppt->tag & M_BDRY) )   return(0);

  //printf("******************* on cherche les tr contenant %d\n",nump);

  /*initialiser la liste*/
  for (l=1; l<=lontet; l++) {
    kk  = listtet->tetra[l] / 4;
    ptk = &mesh->tetra[kk];
    nk  = listtet->tetra[l] % 4;
    /*tts les faces sauf nk contiennent nump*/
    for(i=0 ; i<4 ; i++) {
      if(i==nk) continue;
      iadr = (kk-1)*4 + 1;
      adja = &mesh->adja[iadr];
      if(!adja[i] /*|| (ptk->ref != mesh->tetra[adja[i]/4].ref)*/) break;
    }
    if(i<4) break;
  }
  if(l==(lontet+1)) {
    printf("NO START\n");
    return(0);
  }
  depart = kk;
  fdep   = i;
  pt     = &mesh->tetra[kk];
  //printf("depart : face %d de %d : %d %d %d\n",fdep,depart,
  //                            pt->v[MMG_idir[fdep][0]],pt->v[MMG_idir[fdep][1]],pt->v[MMG_idir[fdep][2]]);

  ilist    = 0;
  /* NON ORIENTED
     for(i=0 ; i<3 ; i++) {
     if(MMG_idir[fdep][i] == nk) continue;
     list->tetra[++ilist] = pt->v[MMG_idir[fdep][i]];
     if(ilist == 1) kfa = MMG_idir[fdep][i];

     }*/

  for(i=0 ; i<3 ; i++) {
    if(MMG_idir[fdep][i] == nk) break;
  }
  assert(i!=3);
  ik = i;
  for(i=ik+1 ; i<ik+3 ; i++) {
    assert(MMG_idir[fdep][(i)%3] != nk);
    list->tetra[++ilist] = pt->v[MMG_idir[fdep][(i)%3]];
    if(ilist == 1) kfa = MMG_idir[fdep][(i)%3];

  }

  /*par adjacence remplir la liste*/
  nobj  = list->tetra[1];
  kel   = depart;
  nlast = list->tetra[ilist];
  while(nlast != nobj) {
    iadr  = (kel-1)*4 + 1;
    adja  = &mesh->adja[iadr];
    jel   = adja[kfa]/4;
    nk    = adja[kfa]%4;
    if(!jel) {
      /*pas de voisin => le tet a une 2e face de bord*/
      ptk = &mesh->tetra[kel];
      kk = -1;
      for(i=0 ; i<3 ; i++) {
	if(ptk->v[MMG_idir[kfa][i]] == nlast) {
	  kk = MMG_idir[kfa][i];
	  continue;
	}
	if(ptk->v[MMG_idir[kfa][i]] == nump) continue;
	list->tetra[++ilist] = ptk->v[MMG_idir[kfa][i]];
      }
      assert(kk!=-1);
      kfa   = kk;
      nlast = list->tetra[ilist];
      continue;
    }
    kel   = jel;
    ptk   = &mesh->tetra[kel];
    /*trouver kfa :  face list->tetra[2] nump ptk->v[nk]*/
    for(i = 0 ; i<4 ; i++) {
      if(i==nk) continue;
      i0 = MMG_idir[i][0];
      i1 = MMG_idir[i][1];
      i2 = MMG_idir[i][2];
      if((i0 != nk) && (i1 != nk) && (i2 != nk) ) continue;
      i0 = ptk->v[i0];
      i1 = ptk->v[i1];
      i2 = ptk->v[i2];
      if((i0 != nump) && (i1 != nump) && (i2 != nump) ) continue;
      if((i0 != nlast) && (i1 != nlast) && (i2 != nlast) ) continue;
      break;
    }
    assert(i<4);
    kfa = i;
    /*test si elle est frontiere*/
    iadr = (kel-1)*4 + 1;
    adja = &mesh->adja[iadr];
    if(!adja[kfa] || (ptk->ref != mesh->tetra[adja[kfa]/4].ref)) {
      list->tetra[++ilist] = ptk->v[nk];

      for(i=0 ; i<3 ; i++) {
	if(ptk->v[MMG_idir[kfa][i]] != nlast) continue;
	kfa = MMG_idir[kfa][i];
	break;
      }
      assert(i<3);
      /*       puts("------------------ apres rajout");
	       for (l=1; l<=ilist; l++) {
	       printf(" %d ",list->tetra[l]);
	       }
	       puts("------------------");
      */
      nlast = list->tetra[ilist];
    }






  } /*while(nlast != nobj);*/

  return(ilist);
}
int MMG_calnormal(pMesh mesh,int p1,int p2,int p3,double* normal) {
  double     ux,uy,uz,vx,vy,vz,wx,wy,wz;
  double     dd,dd1,dd2,dd3,peri,aire;
  double    *a,*b,*c;
  int        ia,ib,ic;

  ia = p1; ib = p2; ic = p3;
  //MMG_ord(&ic,&ib,&ia);
  //printf("----------------- point : %d %d %d\n",ia,ib,ic);

  a = mesh->point[ia].c;
  b = mesh->point[ib].c;
  c = mesh->point[ic].c;


  ux = b[0] - a[0];
  uy = b[1] - a[1];
  uz = b[2] - a[2];

  vx = c[0] - a[0];
  vy = c[1] - a[1];
  vz = c[2] - a[2];

  /* normal */
  normal[0]  = uy*vz - uz*vy;
  normal[1]  = uz*vx - ux*vz;
  normal[2]  = ux*vy - uy*vx;
  dd    = normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2];
  if ( dd <= 0. )  {
    puts("mauvaise orientation!!!!!");
    return(0);
  }
  dd    = 1.0 / sqrt(dd);
  normal[0] *= dd;
  normal[1] *= dd;
  normal[2] *= dd;

  /*aire*/
  wx = c[0] - b[0];
  wy = c[1] - b[1];
  wz = c[2] - b[2];

  dd1 = ux*ux + uy*uy + uz*uz;
  dd1 = sqrt(dd1);

  dd2 = vx*vx + vy*vy + vz*vz;
  dd2 = sqrt(dd2);

  dd3 = wx*wx + wy*wy + wz*wz;
  dd3 = sqrt(dd3);

  peri = 0.5 * (dd1 + dd2 + dd3);
  aire = peri * (peri-dd1) * (peri-dd2) * (peri-dd3);
  if ( aire <= 0.0 )  return(0);
  aire = sqrt(aire);
  normal[0] *= aire;
  normal[1] *= aire;
  normal[2] *= aire;
  //printf("normal %e %e %e\n",normal[0],normal[1],normal[2]);
  return(1);
}

int MMG_ord(int* i1,int* i2 , int* i3) {
  int it1,it2;

  if(*i1 > *i2) {
    it1 = *i2;
    *i2  = *i1;
    *i1  = it1;
  }
  if(i3 < i1) {
    it1 = *i1;
    *i1  = *i3;
    it2 = *i2;
    *i2  = it1;
    *i3  = it2;
  } else if(i3 < i2){
    it2  = *i2;
    *i2  = *i3;
    *i3  = it2;
  }

  return(1);
}

double MMG_caltri(pMesh mesh,pSol sol,int p1,int p2,int p3,double *aire) {
  double     cal,ux,uy,uz,vx,vy,vz,wx,wy,wz;
  double     dd,dd1,dd2,dd3,peri;
  double    *a,*b,*c,n[3],cotmax;
  int        ia,ib,ic;

  ia = p1; ib = p2; ic = p3;
  MMG_ord(&ic,&ib,&ia);
  //printf("----------------- point : %d %d %d\n",ia,ib,ic);

  cal = 1.0e20;

  a = mesh->point[ia].c;
  b = mesh->point[ib].c;
  c = mesh->point[ic].c;


  ux = b[0] - a[0];
  uy = b[1] - a[1];
  uz = b[2] - a[2];

  vx = c[0] - a[0];
  vy = c[1] - a[1];
  vz = c[2] - a[2];

  /* normal */
  n[0]  = uy*vz - uz*vy;
  n[1]  = uz*vx - ux*vz;
  n[2]  = ux*vy - uy*vx;
  dd    = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if ( dd <= 0. )  {
    //puts("mauvaise orientation!!!!!");
    return(cal);
  }
  /*dd    = 1.0 / sqrt(dd);
    n[0] *= dd;
    n[1] *= dd;
    n[2] *= dd;
  */

  /* edge lengths */
  wx = c[0] - b[0];
  wy = c[1] - b[1];
  wz = c[2] - b[2];

  dd1 = ux*ux + uy*uy + uz*uz;
  dd1 = sqrt(dd1);

  dd2 = vx*vx + vy*vy + vz*vz;
  dd2 = sqrt(dd2);

  dd3 = wx*wx + wy*wy + wz*wz;
  dd3 = sqrt(dd3);

  cotmax = dd1;
  if (dd2 > cotmax) cotmax = dd2;
  if (dd3 > cotmax) cotmax = dd3;

  /* quality */
  peri = 0.5 * (dd1 + dd2 + dd3);
  if ( peri < EPS1 )  return(cal);
  *aire = peri * (peri-dd1) * (peri-dd2) * (peri-dd3);
  if ( *aire <= 0.0 )  return(cal);
  *aire = sqrt(*aire);
  cal  = (cotmax * peri) / *aire;//sqrt(aire * peri) / peri;

  return(cal);
}



/* optimise point on boundary */
int MMG_optlenbdry_iso(pMesh mesh,pSol sol,double declic,int base) {
  pTetra        pt,pt1;
  pPoint        ppa,ppb,ppc;
  pDispl                pdisp;
  List                  listtet,listp;
  double        oldc[3],cal,calt,ctg,cx,cy,cz,ux,uy,uz,cpx,cpy,cpz,coe;
  double                dd,len,lmi,lma,normal[3],n[3],cop;
  double        hb,hp,*ca,*cb;
  double                maxdev,dev,air;
  int           i,j,k,l,iel,lontet,lon,nm;
  int           ipa,ipb,ipc,nb,nk,npp,iadr,iter,maxtou,it,maxit;

  /*for normal save*/
  maxtou = 5;
  maxit  = 30;
  it = 0;
  nm = 1;
  base = mesh->flag;

  pdisp = mesh->disp;

  while((it++ < maxit) && nm) {
    npp  = 0;
    for(k=1 ; k<=mesh->ne ; k++) {
      pt = &mesh->tetra[k];
      if(!pt->v[0]) continue;
      for(j=0 ; j<4 ; j++) {
	ipa  = pt->v[j];
	ppa = &mesh->point[ipa];

	if(ppa->geom) continue;
	if(!(ppa->tag & M_BDRY)) continue;
	if(ppa->flag != base) continue;
	/*boule du point*/
	lontet = MMG_bouleg(mesh,k,j,&listtet);
	lon    = MMG_boulebdry(mesh,lontet,&listtet,&listp);
	if(!lon) continue;
	npp++;

	iadr = (ipa-1)*sol->offset + 1;
	hp   = sol->met[iadr];

	cal  = 0.;

	/*optimal point*/
	nb = 0;
	cx   = 0.0;
	cy   = 0.0;
	cz   = 0.0;
	lmi  = LSHORT;
	lma  = LLONG;
	for (l=1; l<=lon; l++) {
	  ipb = listp.tetra[l];
	  ppb = &mesh->point[ipb];
	  if( l < lon ) {
	    calt = MMG_caltri(mesh,sol,ipa,ipb,listp.tetra[l+1],&air);
	    if ( calt > cal )  cal = calt;
	  }
	  iadr = (ipb-1)*sol->offset + 1;
	  hb   = sol->met[iadr];

	  ca  = &mesh->point[ipa].c[0];
	  cb  = &mesh->point[ipb].c[0];
	  len = MMG_length(ca,cb,&hp,&hb);

	  //len = MMG_length(mesh,sol,ipa,ipb);
	  if ( len < lmi )       lmi = len;
	  else if ( len > lma )  lma = len;

	  ux  = ppb->c[0] - ppa->c[0];
	  uy  = ppb->c[1] - ppa->c[1];
	  uz  = ppb->c[2] - ppa->c[2];

	  len = 1.0 / len;
	  cx += ppa->c[0] + ux * len;
	  cy += ppa->c[1] + uy * len;
	  cz += ppa->c[2] + uz * len;
	  nb++;
	}
	assert(nb==lon);
	dd  = 1.0 / (double)nb;
	cpx = cx*dd - ppa->c[0];
	cpy = cy*dd - ppa->c[1];
	cpz = cz*dd - ppa->c[2];

	/* adjust position */
	coe  = HQCOEF;
	iter = 1;
	if ( cal > 10.0 / ALPHATRI )
	  ctg = cal * HCRIT;
	else if ( cal > 5.0 / ALPHATRI )
	  ctg = cal * 0.95;
	else
	  ctg = cal * 0.9975;


	memcpy(oldc,ppa->c,3*sizeof(double));

	/*proj du point sur le plan tangent*/
	/*normales*/
	for(i=0 ; i<3 ; i++) n[i] = 0;

	maxdev = EPSDEV;
	for (l=1; l<lon; l++) {
	  ipb  = listp.tetra[l];
	  ipc  = listp.tetra[l+1];
	  iadr = (ipb-1)*sol->offset + 1;
	  hb   = sol->met[iadr];
	  ppb = &mesh->point[ipb];
	  ppc = &mesh->point[ipc];

	  /*normal*/
	  if(!MMG_calnormal(mesh,ipa,ipb,ipc,&normal[0])) printf("pbs normale\n");

	  /*deviation*/
	  ux     = ppb->c[0] - ppa->c[0];
	  uy     = ppb->c[1] - ppa->c[1];
	  uz     = ppb->c[2] - ppa->c[2];
	  dd     = ux*ux + uy*uy + uz*uz;
	  dev    = (normal[0]*ux + normal[1]*uy + normal[2]*uz)/dd;
	  //        printf("dev %e : %e %e %e -- %e -- %e %e %e\n",dev,normal[0],normal[1],normal[2],dd,ux,uy,uz);
	  // printf("dev2 %e %e %e -- %e\n",normal[0]*ux , normal[1]*uy,normal[2]*uz,normal[0]*ux + normal[1]*uy + normal[2]*uz);
	  // dd = 1./dd;
	  // printf("devv %e %e\n",dd,(normal[0]*ux + normal[1]*uy + normal[2]*uz)*dd);
	  maxdev = M_MAX(maxdev,dev);


	  for(i=0 ; i<3 ; i++) n[i] += normal[i];
	}

	// printf("dev max = %e %e\n",180*asin(maxdev)/M_PI,180*asin(EPSDEV)/M_PI);
	/*normalize*/
	dd = 0;
	for(i=0 ; i<3 ; i++) dd += n[i]*n[i];
	dd = sqrt(dd);
	for(i=0 ; i<3 ; i++) n[i] /= dd;
	for(i=0 ; i<3 ; i++) pdisp->mv[3*(ipa-1) + i] = n[i];
	do {
	  ppa->c[0] = oldc[0] + coe * cpx;
	  ppa->c[1] = oldc[1] + coe * cpy;
	  ppa->c[2] = oldc[2] + coe * cpz;

	  /*projection*/
	  cop = 0;
	  for(i=0 ; i<3 ; i++) cop += (ppa->c[i] * n[i]);
	  for(i=0 ; i<3 ; i++) cop -= (oldc[i] * n[i]);

	  for(i=0 ; i<3 ; i++)
	    ppa->c[i] -= cop * n[i];


	  for (l=1; l<lon; l++) {
	    ipb  = listp.tetra[l];
	    ipc  = listp.tetra[l+1];
	    iadr = (ipb-1)*sol->offset + 1;
	    hb   = sol->met[iadr];

	    /*check length : edge ipa ipb*/
	    ca  = &mesh->point[ipa].c[0];
	    cb  = &mesh->point[ipb].c[0];
	    len = MMG_length(ca,cb,&hp,&hb);
	    //len = MMG_length(mesh,sol,ipa,ipb);
	    if ( len < lmi || len > lma ) break;

	    /*check deviation :ipa ipb ipc*/
	    /*normal*/
	    if(!MMG_calnormal(mesh,ipa,ipb,ipc,&normal[0])) printf("pbs normale\n");

	    /*deviation*/
	    ppb    = &mesh->point[ipb];
	    ux     = ppb->c[0] - ppa->c[0];
	    uy     = ppb->c[1] - ppa->c[1];
	    uz     = ppb->c[2] - ppa->c[2];
	    dd     = ux*ux + uy*uy + uz*uz;
	    dev    = (normal[0]*ux + normal[1]*uy + normal[2]*uz)/dd;
	    if(dev > maxdev ) {
	      printf("dev incorrect %e > %e\n",180*asin(dev)/M_PI,180*asin(maxdev)/M_PI);
	    }


	    /*check triangle qual :ipa ipb ipc*/
	    cal = MMG_caltri(mesh,sol,ipa,ipb,ipc,&air);
	    if ( cal > ctg )  break;
	  }
	  if(l==lon) {
	    /*check last edge*/
	    ipb = listp.tetra[lon];
	    iadr = (ipb-1)*sol->offset + 1;
	    hb  = sol->met[iadr];

	    /*check length : edge ipa ipb*/
	    ca  = &mesh->point[ipa].c[0];
	    cb  = &mesh->point[ipb].c[0];
	    len = MMG_length(ca,cb,&hp,&hb);
	    //len = MMG_length(mesh,sol,ipa,ipb);
	    if ( !(len > lmi && len < lma) ){
	      coe*=0.5;
	      continue;
	    }
	  } else {
	    coe*=0.5;
	    continue;
	  }
	  /*check tet degrad : par tet*/
	  for (l=1; l<=lontet; l++) {
	    iel = listtet.tetra[l] >> 2;
	    nk  = listtet.tetra[l] % 4;
	    pt1 = &mesh->tetra[iel];

	    cal = MMG_caltet(mesh,sol,iel);
	    if ( cal > QDEGRADBDRY )  break;
	    listtet.qual[l] = cal;
	  }
	  if ( l > lontet )  break;

	  coe *= 0.5;

	}
	while ( ++iter <= maxtou );

	if ( iter > maxtou ) {
	  //printf("on bouge pas!\n");
	  memcpy(ppa->c,oldc,3*sizeof(double));
	  continue;
	}

	nm++;
	//printf("---------------------------- on bouge\n");
	//exit(0);
	/* update tetra */
	for (l=1; l<=lontet; l++) {
	  iel = listtet.tetra[l] >> 2;
	  nk  = listtet.tetra[l] % 4;
	  pt1 = &mesh->tetra[iel];
	  pt1->qual = listtet.qual[l];
	  pt1->flag = mesh->flag;
	}

	/* interpol metric */
	ppa->flag = base + 1;

      }/*end j*/
    }/*end k*/

    base++;
    if ( mesh->info.imprim < - 4 )
      fprintf(stdout,"BDRY       %7d PROPOSED  %7d MOVED\n",npp,nm);

    //saveMesh(mesh,"normal.mesh");
    //saveVect(mesh,"normal");
    //return(1);
  }/*end it*/

  if(it==maxit) puts("iter max atteint");



  return(--nm);
}
