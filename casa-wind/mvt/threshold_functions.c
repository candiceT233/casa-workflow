/* $Id: threshold_functions.c,v 0.9 2014-12-10 11:16:23 elyons Exp $ */
/* Copyright 2014 University of Massachusetts Amherst all rights reserved*/

/* threshold_functions.c contains a number of common functions used in the UMass suite of thresholding algorithms */

#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <math.h>
#include "threshold_functions.h"

void handle_error(int errid){
  if(errid == NC_NOERR) return;
  printf("%s\n", nc_strerror(errid));
  return;
}

void malloc_error(char *In)
{
  printf("Error allocating memory for %s\n",In);
  exit(1);
};

double degToRad(double degs) {
  double rad = 0;
  rad = degs * (M_PI / 180);
  return rad;
};

double radToDeg(double rads) {
  double deg = 0;
  deg = rads * (180 / M_PI);
  return deg;
}

int contained(int nps, struct index_pair pair_in_question, struct index_pair *array_of_pairs_in_question) {
  int k;
  int is_contained = 0;
  for (k=0; k<nps; ++k)
    if ((pair_in_question.xind == array_of_pairs_in_question[k].xind) && (pair_in_question.yind == array_of_pairs_in_question[k].yind))
      is_contained = 1;
  return is_contained;
}

void Contour(float **d,int ilb,int iub,int jlb,int jub,
             int *x,int *y,int nc,float *z, struct index_pair ***indp1, struct index_pair ***indp2, struct index_pair ***indp3, struct index_pair ***indp4) {
#define xsect(p1,p2) (h[p2]*xh[p1]-h[p1]*xh[p2])/(h[p2]-h[p1])
#define ysect(p1,p2) (h[p2]*yh[p1]-h[p1]*yh[p2])/(h[p2]-h[p1])
  int m1,m2,m3,case_value;
  float dmin,dmax;
  int x1=0,x2=0,y1=0,y2=0;
  int i,j,k,m;
  float h[5];
  int sh[5];
  float xh[5],yh[5];
  int im[4] = {0,1,1,0},jm[4]={0,0,1,1};
  int castab[3][3][3] = {
    { {0,0,8},{0,2,5},{7,6,9} },
    { {0,3,4},{1,3,1},{4,3,0} },
    { {9,6,7},{5,2,0},{8,0,0} }
  };
  float temp1,temp2;
  for (j=(jub-1);j>=jlb;j--) {
    for (i=ilb;i<=iub-1;i++) {
      temp1 = fminf(d[i][j],d[i][j+1]);
      temp2 = fminf(d[i+1][j],d[i+1][j+1]);
      dmin  = fminf(temp1,temp2);
      temp1 = fmaxf(d[i][j],d[i][j+1]);
      temp2 = fmaxf(d[i+1][j],d[i+1][j+1]);
      dmax  = fmaxf(temp1,temp2);
      //printf("mag %d %d %f\n", i, j, d[i][j]);                                                                                                                                                         
      //printf("dmin %f dmax %f\n", dmin, dmax);                                                                                                                                                         
      if (dmax < z[0] || dmin > z[nc-1])
	continue;
      for (k=0;k<nc;k++) {
        //printf("k %d zk %f\n", k, z[k]);                                                                                                                                                               
        if (z[k] < dmin || z[k] > dmax)
          continue;
        for (m=4;m>=0;m--) {
          if (m > 0) {
            h[m]  = d[i+im[m-1]][j+jm[m-1]]-z[k];
            xh[m] = x[i+im[m-1]];
            yh[m] = y[j+jm[m-1]];
          } else {
            h[0]  = 0.25 * (h[1]+h[2]+h[3]+h[4]);
            xh[0] = 0.50 * (x[i]+x[i+1]);
            yh[0] = 0.50 * (y[j]+y[j+1]);
          }
          if (h[m] > 0.0)
            sh[m] = 1;
          else if (h[m] < 0.0)
            sh[m] = -1;
          else
            sh[m] = 0;
        }
	for (m=1;m<=4;m++) {
          m1 = m;
          m2 = 0;
          if (m != 4)
            m3 = m + 1;
          else
            m3 = 1;
          if ((case_value = castab[sh[m1]+1][sh[m2]+1][sh[m3]+1]) == 0)
	    continue;
          switch (case_value) {
          case 1: /* Line between vertices 1 and 2 */
            x1 = xh[m1];
            y1 = yh[m1];
            x2 = xh[m2];
            y2 = yh[m2];
            break;
          case 2: /* Line between vertices 2 and 3 */
            x1 = xh[m2];
            y1 = yh[m2];
            x2 = xh[m3];
            y2 = yh[m3];
            break;
          case 3: /* Line between vertices 3 and 1 */
            x1 = xh[m3];
            y1 = yh[m3];
            x2 = xh[m1];
            y2 = yh[m1];
            break;
          case 4: /* Line between vertex 1 and side 2-3 */
            x1 = xh[m1];
            y1 = yh[m1];
            x2 = xsect(m2,m3);
            y2 = ysect(m2,m3);
            break;
	  case 5: /* Line between vertex 2 and side 3-1 */
            x1 = xh[m2];
            y1 = yh[m2];
            x2 = xsect(m3,m1);
	    y2 = ysect(m3,m1);
            break;
          case 6: /* Line between vertex 3 and side 1-2 */
            x1 = xh[m3];
            y1 = yh[m3];
            x2 = xsect(m1,m2);
            y2 = ysect(m1,m2);
            break;
          case 7: /* Line between sides 1-2 and 2-3 */
            x1 = xsect(m1,m2);
            y1 = ysect(m1,m2);
            x2 = xsect(m2,m3);
            y2 = ysect(m2,m3);
            break;
          case 8: /* Line between sides 2-3 and 3-1 */
            x1 = xsect(m2,m3);
            y1 = ysect(m2,m3);
            x2 = xsect(m3,m1);
            y2 = ysect(m3,m1);
            break;
          case 9: /* Line between sides 3-1 and 1-2 */
            x1 = xsect(m3,m1);
            y1 = ysect(m3,m1);
            x2 = xsect(m1,m2);
            y2 = ysect(m1,m2);
            break;
	  default:
            break;
          }
          //printf("i %d j %d x1 %d, y1 %d, x2 %d, y2 %d, zk %f\n", i, j, x1, y1, x2, y2, z[k]);
	  //ip1 is who points to me
          struct index_pair ip1;
          //ip2 is who I point to
          struct index_pair ip2;
          if ((x1 != x2) || (y1 != y2)) {
            if((i != x1) || (j != y1)) {
              ip2.xind = x1;
              ip2.yind = y1;
              ip1.xind = x2;
              ip1.yind = y2;
            }
            else {
              ip2.xind = x2;
              ip2.yind = y2;
              ip1.xind = x1;
              ip1.yind = y1;
            }
            if (indp2[k][i][j].xind == -1)
              indp2[k][i][j] = ip2;
            else
              indp4[k][i][j] = ip2;
            if (indp1[k][ip2.xind][ip2.yind].xind == -1)
              indp1[k][ip2.xind][ip2.yind] = ip1;
            else
              indp3[k][ip2.xind][ip2.yind] = ip1;
          }
	}
      }
    }
  }
}

int vincenty(double inLat, double inLon, double inDx, double inBrng, struct latLonBrng *llb_out) {
  /* WS-84 Ellipsoid */

  double n1 = 6378137;
  double n2 = 6356752.3142;
  double n3 = 1/298.257223563;

  double inBrngRad = degToRad(inBrng);
  double sinBrng = sin(inBrngRad);
  double cosBrng = cos(inBrngRad);

  double f1 = (1-n3) * tan(degToRad(inLat));
  double f2 = 1/(sqrt((1 + (f1 * f1))));
  double f3 = f1 * f2;
  double f4 = atan2(f1, cosBrng);
  double f5 = f2 * sinBrng;
  double f6 = 1 - (f5 * f5);
  double f7 = f6 * (((n1 * n1) - (n2 * n2))/ (n2 * n2));
  double f8 = 1 + ((f7/16384) * (4096 + (f7 * (-768 + (f7 * (320 - (175*f7)))))));
  double f9 = (f7/1024) * (256 + (f7 * (-128 + (f7 * (74 - (47 * f7))))));
  double f10 = inDx / (n2 * f8);
  double f10a = 2 * M_PI;
  double f11 = 0; double f12 = 0; double f13 = 0; double f14 = 0.0;
  while (fabs(f10 - f10a) > pow(10,-12)) {
    f11 = cos((2*f4) + f10);
    f12 = sin(f10);
    f13 = cos(f10);
    f14 = f9*f12*(f11 + f9/4*(f13*(-1+2*f11*f11)-f9/6*f11*(-3+4*f12*f12)*(-3+4*f11*f11)));
    f10a = f10;
    f10 = (inDx / (n2 * f8)) + f14;
  }
  double f15 = (f3 * f12) - (f2 * f13 * cosBrng);
  double outLat = atan2(((f3*f13) + (f2*f12*cosBrng)), ((1-n3) * sqrt((f5*f5) + (f15*f15))));
  double f16 = atan2((f12*sinBrng), ((f2*f13) - (f3 * f12 * cosBrng)));
  double f17 = n3/16*f6*(4+(n3*(4-(3*f6))));
  double f18 = f16 - ((1-f17) * n3 * f5 * (f10 + (f17 * f12 * (f11 + (f17 * f13 * (-1 + (2 * f11 * f11)))))));
  double outLon = fmod((degToRad(inLon) + f18 + (3 * M_PI)), (2 * M_PI)) - M_PI;
  double f19 = atan2(f5, -f15);

  llb_out->lat = radToDeg(outLat);
  llb_out->lon = radToDeg(outLon);
  llb_out->brng = radToDeg(f19);
  //printf("inlat: %f, inlon: %f, inDx: %f, inBrng: %f, outlat: %f, outlon: %f outBrng: %f\n", inLat, inLon, inDx, inBrng, llb_out->lat, llb_out->lon, llb_out->brng);                          

  return 0;
};

int inverse_vincenty(double inLat1, double inLon1, double inLat2, double inLon2, double *displacement)
{
  /* WS-84 Ellipsoid */
  double n1 = 6378137;
  double n2 = 6356752.3142;
  double n3 = 1/298.257223563;

  double diffLon = degToRad(inLon2 - inLon1);
  double f1 = atan(tan(degToRad(inLat1)) * (1-n3));
  double sinf1 = sin(f1); double cosf1 = cos(f1);
  double f2 = atan(tan(degToRad(inLat2)) * (1-n3));
  double sinf2 = sin(f2); double cosf2 = cos(f2);
  double f3 = diffLon;
  double f3new, sinf3, cosf3, f4, f5, f6, f7, f8, f9, f10;
  do {
    sinf3 = sin(f3); cosf3 = cos(f3);
    f4 = sqrt(((cosf2*sinf3) * (cosf2*sinf3)) + (((cosf1*sinf2)-(sinf1*cosf2*cosf3)) * ((cosf1*sinf2) - (sinf1*cosf2*cosf3))));
    if (f4 == 0) {
      (*displacement) = 0;
      return 0;
    }
    f5 = (sinf1*sinf2) + (cosf1*cosf2*cosf3);
    f6 = atan2(f4,f5);
    f7 = (cosf1 * cosf2 * sinf3)/f4;
    f8 = 1 - (f7 * f7);
    f9 = f5 - ((2 * sinf1 * sinf2)/f8);
    if (isnan(f9))
      f9 = 0;
    f10 = (n3/16)*f8*(4+(n3*(4 - 3*f8)));
    f3new = f3;
    f3 = diffLon + ((1-f10) * n3 * f7 * (f6 + (f10*f4*(f9+(f10*f5*(-1 + (2*f9*f9)))))));
  } while (abs(f3-f3new) > 1e-12);

  double f11 = (f8 * ((n1 * n1) - (n2 * n2))) / (n2 * n2);
  double f12 = 1 + (f11/16384) * (4096 + (f11 * (-768 + (f11 * (320 - (175*f11))))));
  double f13 = (f11/1024) * (256 + (f11 * (-128 + (f11 * (74 - (47 * f11))))));
  double f14 = f13 * f4 * (f9 + (f13/4) * ((f5 * (-1 + (2 * f9 * f9))) - ((f13/6) * f9 * (-3 + (4 * f4 * f4)) * (-3 + (4 * f9 * f9)))));
  double f15 = n2*f12*(f6-f14);
  (*displacement) = f15;
  return 0;
  // optional return bearing from point 1 to point 2 for the future
  // double bearingrad = atan2((cosf2*sinf3), ((cosf1*sinf2) - (sinf1*cosf2*cosf3)));
  // double bearingdeg = radToDeg(bearingrad); 
}

int find_min_dx(struct latLon *ll_arr, int ll_arr_point, int ll_arr_size, int *contour_arr, int contour_arr_size, int *min_dx_index, double *min_disp) {
  double min_dx = 999999999; //some big number to initialize;                                                                                                                                   
  double min_index = -1;
  int counter;
  int counter2;
  int contained = 0;
  double disp;
  for(counter = 0; counter < ll_arr_size; ++counter) {
    //printf("counter:%d\n", counter);
    if (counter == ll_arr_point) {
      //printf("comparing the same point... moving on\n");
      continue;
    }
    for (counter2 = 0; counter2 < contour_arr_size; ++counter2) {
      //printf("counter:%d contour_arr[%d]: %d\n", counter, counter2, contour_arr[counter2]);
      if (counter == contour_arr[counter2]) {
        //printf("this point already exists in the contour... moving on\n");
        contained = 1;
        break;
      }
    }
    if (contained == 1) {
      contained = 0;
      continue;
    }

    inverse_vincenty(ll_arr[ll_arr_point].lat, ll_arr[ll_arr_point].lon, ll_arr[counter].lat, ll_arr[counter].lon, &disp);
    //printf("x:%d, xlat:%f xlon:%f counter:%d counterlat:%f counterlon:%f disp:%f\n", ll_arr_point, ll_arr[ll_arr_point].lat, ll_arr[ll_arr_point].lon, counter, ll_arr[counter].lat, ll_arr[counter].lon, disp);
    if(disp < min_dx) {
      min_dx = disp;
      min_index = counter;
    }
    
  }
  (*min_dx_index) = min_index;
  (*min_disp) = min_dx;
  
  return 0;
}

void makeDiscreteContours(size_t nlats, float *latarr, size_t nlons, float *lonarr, int nc, struct index_pair ***indp1, struct index_pair ***indp2, struct index_pair ***indp3, struct index_pair ***indp4) {
  int h, i, j;  //counters
  
  //for each contour level that was produced... for example 35dbZ, then 45dbZ, and so on...
  for (h=0; h<nc; h++) {
    int k;
    int numpoints=0; //number of points in this particular contour
    int num_used_up_ips=0; //number of points total in all contours
    int numcontours=0; //number of discrete closed contours
    struct latLon *ll_array;
    ll_array = (struct latLon*)malloc(sizeof(struct latLon));
    int edgecount = 0;
    
    struct index_pair reference; //the i&j index reference of the cell that pointed to the current cell
    int ref_index; //the array that pointed the reference to the current_ip 0=indp1 1=indp2 2=indp3 3=indp4
    reference.xind = -1;
    reference.yind = -1;

    struct index_pair *used_up_ips; //an array of all the index points we've used for all contours
    struct index_pair *contour_ips; //an array of the index points in the current contour
    used_up_ips = (struct index_pair*)malloc(sizeof(struct index_pair));
    contour_ips = (struct index_pair*)malloc(sizeof(struct index_pair));

    for (i=0; i<nlats;i++) {
      for (j=0; j<nlons;j=j+0) {
	
	struct index_pair current_ip;
	current_ip.xind = i;
	current_ip.yind = j;

        //print out the arrays
	//printf("%d %d %d %d %d %d %d %d %d %d\n", i, j, indp1[h][i][j].xind, indp1[h][i][j].yind, indp3[h][i][j].xind, indp3[h][i][j].yind, indp2[h][i][j].xind, indp2[h][i][j].yind, indp4[h][i][j].xind, indp4[h][i][j].yind);

	//each point in the grid has a three possibilities, post conrec.

	//firstly, we may have already dealt with it... if so we move on
	if (contained(num_used_up_ips, current_ip, used_up_ips)) {
	  ++j;
	  //printf("i:%d j:%d is contained... moving on\n", i, j);
	  continue;
	}

	//secondly, it may not be part of a contour... if so, we add it to the list of points we've already dealt with and move on
	if ((indp1[h][i][j].xind == -1) && (indp2[h][i][j].xind == -1) && (indp3[h][i][j].xind == -1) && (indp4[h][i][j].xind == -1)) {
	  //printf("i:%d j:%d is not part of a contour\n", i, j);
	  used_up_ips[num_used_up_ips++] = current_ip;
	  num_used_up_ips++;
	  used_up_ips = (struct index_pair*)realloc(used_up_ips, (num_used_up_ips+1)*sizeof(struct index_pair));
	  ++j;
	  continue;
	}

	//thirdly, it may be both new AND part of a contour...

	//some points have multiple paths... these are harder to deal with... we want to eliminate them asap
	int multipath = 0;
	int kntr = 0;
	if (indp1[h][i][j].xind != -1)
	  ++kntr;
	if (indp3[h][i][j].xind != -1)
	  ++kntr;
	if (indp2[h][i][j].xind != -1)
	  ++kntr;
	if (indp4[h][i][j].xind != -1)
	  ++kntr;
	if (kntr > 2) {
	  multipath = 1;
	}

	//now that that's out of the way, let's see if we have a point that referred us here...
	//if not this must be the first point in a new contour
	if(reference.xind == -1) {
	  //we don't want the first point to be a multipath... if it is, toss it back and move on
	  if (multipath == 1) {
            //printf("i:%d j:%d is a multipoint... pass\n", i, j);
	    ++j;
	    continue;
	  }
	  
	  //otherwise start a new contour
	  //printf("i:%d j:%d is the first point in a new contour\n", i, j);

	  //increment the number of points in this contour
	  //numpoints++;
	  //set the first point in this contour to this index point

	  contour_ips[0] = current_ip;

	  //now we start down the path
	  if(indp3[h][i][j].xind != -1) {
	    i = indp3[h][i][j].xind;
	    j = indp3[h][current_ip.xind][j].yind;
	    ref_index = 2;
	  }
	  else if (indp1[h][i][j].xind != -1){
	    i = indp1[h][i][j].xind;
	    j = indp1[h][current_ip.xind][j].yind;
	    ref_index = 0;
	  }
	  else if (indp4[h][i][j].xind != -1){
	    i = indp4[h][i][j].xind;
	    j = indp4[h][current_ip.xind][j].yind;
	    ref_index = 3;
	  }
	  else if (indp2[h][i][j].xind != -1){
	    i = indp2[h][i][j].xind;
	    j = indp2[h][current_ip.xind][j].yind;
	    ref_index = 1;
	  }
            
	  //set the reference point and move on
	  reference = current_ip;
          //printf("reference set\n");
	  continue;
	}
	
	//ok so if we made it here we have a contour with at least one point already, and a reference
	//if this point is not already contained in the ip array of this contour, we are not closed yet
	//or it could be a multipath and we intentionally didn't add it to the array
	if (!contained(numpoints+1, current_ip, contour_ips)) {
	  //printf("new point with a reference\n");
	  //first let's figure out which index is the reference...
          //printf("reference.xind: %d reference.yind: %d\n", reference.xind, reference.yind);
	  int rfr;
	  if ((indp1[h][i][j].xind == reference.xind) && (indp1[h][i][j].yind == reference.yind)) {
	    rfr = 0;
	  }
	  else if ((indp2[h][i][j].xind == reference.xind) && (indp2[h][i][j].yind == reference.yind)) {
	    rfr = 1;
	  }
	  else if ((indp3[h][i][j].xind == reference.xind) && (indp3[h][i][j].yind == reference.yind)){
	    rfr = 2;
	  }
	  else if ((indp4[h][i][j].xind == reference.xind) && (indp4[h][i][j].yind == reference.yind)){
	    rfr = 3;
	  }

	  //if it isn't a multipath it is simpler
	  if(multipath == 0) {
	    //printf("this is a non-multipath point in a contour\n");
	    //increase the number of points
	    numpoints++;
	    //increase the size of the dynamic contour array
	    contour_ips = (struct index_pair*)realloc(contour_ips, (numpoints+1)*sizeof(struct index_pair));
	    //add this ip to our list of ips in this contour
	    contour_ips[numpoints] = current_ip;

	    //if this isn't a multipath, there is only one path in and out
	    if ((indp1[h][i][j].xind != -1) && (rfr != 0)) {
	      i = indp1[h][i][j].xind;
	      j = indp1[h][current_ip.xind][j].yind;
	      ref_index = 0;
	    }
	    else if ((indp2[h][i][j].xind != -1) && (rfr != 1)) {
	      i = indp2[h][i][j].xind;
	      j = indp2[h][current_ip.xind][j].yind;
	      ref_index = 1;
	    }
	    else if ((indp3[h][i][j].xind != -1) && (rfr != 2)) {
	      i = indp3[h][i][j].xind;
	      j = indp3[h][current_ip.xind][j].yind;
	      ref_index = 2;
	    }
	    else if ((indp4[h][i][j].xind != -1) && (rfr != 3)) {
	      i = indp4[h][i][j].xind;
	      j = indp4[h][current_ip.xind][j].yind;
	      ref_index = 3;
	    }
	    else
	      printf("point is not currently handled... skipping\n");
	    
	    reference = current_ip;
            //printf("non-multipath point processed\n");
	    continue;
	  }

	  //if we're here it's a multipath point in a contour
	  //printf("this is a multipath point in a contour\n");
	  //we have to toss this point and replace references to this point
	  //with references to another point... probably the closest point by proximity
          //printf("multipath point\n");
	  double dx2p;
	  double mndx = 999999999999.0; //Huge val
	  struct index_pair destination;
	  int which_index;  //0=indp1 1=indp2 2=indp3 3=indp4
	  if ((indp1[h][i][j].xind != -1) && (rfr != 0)) {
	    inverse_vincenty(latarr[reference.xind], lonarr[reference.yind], latarr[indp1[h][i][j].xind], lonarr[indp1[h][i][j].yind], &dx2p);
	    if (dx2p < mndx) {
	      mndx = dx2p;
	      destination.xind = indp1[h][i][j].xind;
	      destination.yind = indp1[h][i][j].yind;
	      which_index = 0;
	    }
	  }

	  if ((indp2[h][i][j].xind != -1) && (rfr != 1)) {
	    inverse_vincenty(latarr[reference.xind], lonarr[reference.yind], latarr[indp2[h][i][j].xind], lonarr[indp2[h][i][j].yind], &dx2p);
	    if (dx2p < mndx) {
	      mndx = dx2p;
	      destination.xind = indp2[h][i][j].xind;
	      destination.yind = indp2[h][i][j].yind;
	      which_index = 1;
	    }
	  }

	  if ((indp3[h][i][j].xind != -1) && (rfr != 2)) {
	    inverse_vincenty(latarr[reference.xind], lonarr[reference.yind], latarr[indp3[h][i][j].xind], lonarr[indp3[h][i][j].yind], &dx2p);
	    if (dx2p < mndx) {
	      mndx = dx2p;
	      destination.xind = indp3[h][i][j].xind;
	      destination.yind = indp3[h][i][j].yind;
	      which_index = 2;
	    }
	  }

	  if ((indp4[h][i][j].xind != -1) && (rfr != 3)) {
	    inverse_vincenty(latarr[reference.xind], lonarr[reference.yind], latarr[indp4[h][i][j].xind], lonarr[indp4[h][i][j].yind], &dx2p);
	    if (dx2p < mndx) {
	      mndx = dx2p;
	      destination.xind = indp4[h][i][j].xind;
	      destination.yind = indp4[h][i][j].yind;
	      which_index = 3;
	    }
	  }

	  //set the reference to the destination
	  switch (ref_index) {
	  case 0:
	    indp1[h][reference.xind][reference.yind] = destination;
	    break;
	  case 1:
	    indp2[h][reference.xind][reference.yind] = destination;
	    break;
	  case 2:
	    indp3[h][reference.xind][reference.yind] = destination;
	    break;
	  case 3:
	    indp4[h][reference.xind][reference.yind] = destination;
	    break;
	  }

	  //change the destination's reference to the current_ip over to the current ip's reference
	  if ((indp1[h][destination.xind][destination.yind].xind == current_ip.xind) && (indp1[h][destination.xind][destination.yind].yind == current_ip.yind))
	    indp1[h][destination.xind][destination.yind] = reference;
	  else if ((indp2[h][destination.xind][destination.yind].xind == current_ip.xind) && (indp2[h][destination.xind][destination.yind].yind == current_ip.yind))
	    indp2[h][destination.xind][destination.yind] = reference;
	  else if((indp3[h][destination.xind][destination.yind].xind == current_ip.xind) && (indp3[h][destination.xind][destination.yind].yind == current_ip.yind))
	    indp3[h][destination.xind][destination.yind] = reference;
	  else if((indp4[h][destination.xind][destination.yind].xind == current_ip.xind) && (indp4[h][destination.xind][destination.yind].yind == current_ip.yind))
	    indp4[h][destination.xind][destination.yind] = reference;
	  else
	    printf("An unhandled state has been reached... skipping\n");

	  //scrub the current_ip's path to the reference and the destination
	  switch (rfr) {
	  case 0:
	    indp1[h][i][j].xind = -1;
	    indp1[h][i][j].yind = -1;
	    break;
	  case 1:
	    indp2[h][i][j].xind = -1;
	    indp2[h][i][j].yind = -1;
	    break;
	  case 2:
	    indp3[h][i][j].xind = -1;
	    indp3[h][i][j].yind = -1;
	    break;
	  case 3:
	    indp4[h][i][j].xind = -1;
	    indp4[h][i][j].yind = -1;
	    break;
	  }
	  switch (which_index) {
	  case 0:
	    indp1[h][i][j].xind = -1;
	    indp1[h][i][j].yind = -1;
	    break;
	  case 1:
	    indp2[h][i][j].xind = -1;
	    indp2[h][i][j].yind = -1;
	    break;
	  case 2:
	    indp3[h][i][j].xind = -1;
	    indp3[h][i][j].yind = -1;
	    break;
	  case 3:
	    indp4[h][i][j].xind = -1;
	    indp4[h][i][j].yind = -1;
	    break;
	  }

	  //we may need to further modify the paths... tbd

	  //move to the destination
	  i = destination.xind;
	  j = destination.yind;
          //printf("multipoint handled\n");
	  continue;
	}

	if (contained(numpoints+1, current_ip, contour_ips)) {
	  //if the point is contained we are at the end.... maybe?
	  //increment the number of contours
	  ++numcontours;
	  
	  for (k=0; k<numpoints+1; ++k) {
	    //print out the contour
	    printf("contour %d pair %d i %d j %d\n", numcontours, k, contour_ips[k].xind, contour_ips[k].yind);
	    used_up_ips[num_used_up_ips] = contour_ips[k];
	    num_used_up_ips++;
	    used_up_ips = (struct index_pair*)realloc(used_up_ips, (num_used_up_ips+1)*sizeof(struct index_pair));
	  }

	  //dump contour to output
	  int a = 0;
	  //if (numpoints > 3) {
	  for (a=0; a < numpoints+1; a++) {
	    ll_array[edgecount].lat = latarr[contour_ips[a].xind];
	    ll_array[edgecount].lon = lonarr[contour_ips[a].yind];
	    ++edgecount;
	    ll_array = (struct latLon*)realloc(ll_array, (edgecount+1)*sizeof(struct latLon));
	  }
	  //now append the first point back to the end per some polygon format standards
	  ll_array[edgecount].lat = latarr[contour_ips[0].xind];
	  ll_array[edgecount].lon = lonarr[contour_ips[0].yind];
	  ++edgecount;
	  ll_array = (struct latLon*)realloc(ll_array, (edgecount+1)*sizeof(struct latLon));
	    
	  //print the contour
	  for (k=0; k<edgecount; k++)
	    printf("contour:%d point:%d lat: %f lon: %f\n", numcontours, k, ll_array[k].lat, ll_array[k].lon);
		  
	  //return to the start of the first contour
	  i = contour_ips[0].xind;
	  j = contour_ips[0].yind + 1;

	  //clear the reference
	  reference.xind = -1;
	  reference.yind = -1;

	  //reset the numpoints to zero
	  numpoints = 0;

	  //reset the edgecount
	  edgecount = 0;

	  //free up and malloc contours_ips and ll_array to the minimum
	  free(contour_ips);
	  contour_ips = (struct index_pair*)malloc(sizeof(struct index_pair));

	  free(ll_array);
	  ll_array = (struct latLon*)malloc(sizeof(struct latLon));

	  continue;
	}
      }
    }
    free(used_up_ips);
  }
}
