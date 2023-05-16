/* $Id: UM_VEL.c,v 0.9 2018-04-27 15:05:16 elyons Exp $ */
/* Copyright 2018 University of Massachusetts Amherst all rights reserved*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <netcdf.h>
#include <libconfig.h>
#include <libxml/parser.h>
#include <libxml/tree.h>

void handle_error(int errid){
  if(errid == NC_NOERR) return;
  printf("%s\n", nc_strerror(errid));
  exit(-1);
  return;
};

struct latLonBrng {
  double lat;
  double lon;
  double brng;
};

struct latLon {
  double lat;
  double lon;
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
};

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

int UM_VEL(size_t count, char* i_filenames[])
{
  //locate and verify the config file
  char *UM_home = getenv("UM_VEL");  
  config_t UM_config;
  char config_fn[128];
  sprintf(config_fn, "%s/UM_VEL_config.txt", UM_home);
  printf("config_fn is %s\n", config_fn);
  
  config_init(&UM_config);
  if (!config_read_file(&UM_config, config_fn)) {
    printf("config file not found or contains errors.  Should be in $UM_VEL directory and called UM_VEL_config.txt. exiting...\n");
    config_destroy(&UM_config);
    return -1;
  }

  //get output directory
  const char *ncdir;
  if (config_lookup_string(&UM_config, "ncdir", &ncdir)) {
    printf("ncdir: %s\n", ncdir);
    if (strlen(ncdir) > 128) 
      printf("ncdir has a 128 character limit;  Please enter a shorter path.  Exiting...\n");
  }
  else {
    printf("ncdir not defined.  Please check config file. Exiting...\n");
    exit(-1);
  }

  //get alg parameters from config file
  //average = 0, max = 1 
  int routine;
  if (config_lookup_int(&UM_config, "routine", &routine)) {
    if (routine == 1)
      printf("routine: max\n");
    else {
      routine = 0;
      printf("routine: mean\n");
    }
  }
  else {
    printf("routine not defined in config file.  Default is mean");
    routine = 0;
  }

  //truemap off = 0, on = 1
  int truemap;
  if (config_lookup_int(&UM_config, "truemap", &truemap)) {
    if (truemap == 1)
      printf("truemap: on\n");
    else {
      truemap = 0;
      printf("truemap: off\n");
    }
  }
  else {
    printf("truemap not defined in config file.  Default is off");
    truemap = 0;
  }

  //smoothing off = 0, on = 1
  int smoothing;
  //double smoothdx = 1000; //meters
  int smoothvoxels=12; 
  if (config_lookup_int(&UM_config, "smoothing", &smoothing)) {
    if (smoothing == 1)
      printf("smoothing: on\n");
    else {
      smoothing = 0;
      printf("smoothing: off\n");
    }
    
    if (!config_lookup_int(&UM_config, "smoothvoxels", &smoothvoxels)) {
      printf("smoothvoxels not defined in config file.  Using 12\n");
    }
    else 
      printf("smoothvoxels: %d\n", smoothvoxels);
    
  }
  else {
    printf("smoothing not defined in config file.  Default is off");
    smoothing = 0;
  }

  //get grid definition from config file
  double nlat, slat, wlon, elon, latspac, lonspac, minalt, maxalt, mintilt, maxtilt, minval, maxval;
  int numlats, numlons;
 
  if (config_lookup_float(&UM_config, "nlat", &nlat)) 
    printf("nlat: %f\n", nlat);
  else {
    printf("nlat not defined.  Please check config file.  Exiting....\n");
    exit(-1);
  }

  if (config_lookup_float(&UM_config, "wlon", &wlon))
    printf("wlon: %f\n", wlon);
  else {
    printf("wlon not defined.  Please check config file.  Exiting....\n");
    exit(-1);
  }
  
  if (config_lookup_float(&UM_config, "latspac", &latspac))
    printf("latspac: %f\n", latspac);
  else {
    printf("latspac not defined.  Please check config file.  Exiting....\n");
    exit(-1);
  }

  if (config_lookup_float(&UM_config, "lonspac", &lonspac))
    printf("lonspac: %f\n", lonspac);
  else {
    printf("lonspac not defined.  Please check config file.  Exiting....\n");
    exit(-1);
  }

  if (config_lookup_float(&UM_config, "minalt", &minalt))
    printf("minalt: %f\n", minalt);
  else {
    printf("minalt not defined.  Please check config file.  Exiting....\n");
    exit(-1);
  }

  if (config_lookup_float(&UM_config, "maxalt", &maxalt))
    printf("maxalt: %f\n", maxalt);
  else {
    printf("maxalt not defined.  Please check config file.  Exiting....\n");
    exit(-1);
  }

  if (config_lookup_float(&UM_config, "mintilt", &mintilt))
    printf("mintilt: %f\n", mintilt);
  else {
    printf("mintilt not defined.  Using 0.0 as the default\n");
    mintilt = 0.0;
  }

  if (config_lookup_float(&UM_config, "maxtilt", &maxtilt))
    printf("maxtilt: %f\n", maxtilt);
  else {
    printf("maxtilt not defined.  Using 2.1 as the default\n");
    maxtilt = 2.1;
  }
  
  if (config_lookup_int(&UM_config, "numlats", &numlats))
    printf("numlats: %d\n", numlats);
  else {
    printf("numlats not defined.  Please check config file.  Exiting....\n");
    exit(-1);
  }

  if (config_lookup_int(&UM_config, "numlons", &numlons))
    printf("numlons: %d\n", numlons);
  else {
    printf("numlons not defined.  Please check config file.  Exiting....\n");
    exit(-1);
  }
  
  slat = nlat - (numlats * latspac);
  elon = wlon + (numlons * lonspac);

  //get name and units of variable to grid
  const char *varname;
  if (config_lookup_string(&UM_config, "varname", &varname)) 
    printf("varname: %s\n", varname);
  else {
    printf("varname not defined.  Please check config file. Exiting...\n");
    exit(-1);
  }

  const char *units;
  if (config_lookup_string(&UM_config, "units", &units))
    printf("units: %s\n", units);
  else {
    printf("units not defined.  Please check config file. Exiting...\n");
    exit(-1);
  }
  
  //get min/max vals

  if (config_lookup_float(&UM_config, "minval", &minval))
    printf("minval: %f\n", minval);
  else {
    printf("minval not defined.  Please check config file.  Exiting....\n");
    exit(-1);
  }

  if (config_lookup_float(&UM_config, "maxval", &maxval))
    printf("maxval: %f\n", maxval);
  else {
    printf("maxval not defined.  Please check config file.  Exiting....\n");
    exit(-1);
  }

  //start the grid process
  float **Data_In; 
  float **SNR_In;
  float **NCP_In;
  float **RHO_In;
  float **REF_In;
  int **Flag_In;
  double *GateWidth_In, *Azimuth_In, *Elevation_In, *StartRange_In;
  int *Time_In;
  float grid[numlats][numlons];
  int sample[numlats][numlons];
  char starttime[16];
  double AntennaBeamwidth;
  int status, ncid, Radial_id, Gate_id, Data_id, SNR_id, NCP_id, RHO_id, REF_id, Time_id, GateWidth_id, Az_id, El_id, Sr_id, GF_id, i, j, k, x, y, ScanType, maxtime;
  //int Vel_id;
  //time_t tm_t_0;
  size_t numgates, numradials;

  //struct tm broketime;
  struct latLonBrng radar_llb, dest_llb;

  for (x = 0; x<numlats; x++) {
    for (y = 0; y<numlons; y++) {
      grid[x][y] = 0;
      sample[x][y] = 0;
    }
  }
  
  maxtime = 0;
  char *pch;
  char tmpname[256];
  memset(tmpname, '\0', sizeof(tmpname));
  strcpy(tmpname, i_filenames[count-1]);
  char ymd[9];
  char hms[7];
  
  pch = strtok (tmpname, "-");
  int spl = 0;
  while (pch != NULL) {
    if (spl == 1) {
      //printf ("%s\n", pch);
      sprintf(ymd, "%s", pch);
    }
    if (spl == 2) {
      char *pch2;
      pch2 = strtok (pch, ".");
      sprintf(hms, "%s", pch2);
    }
    pch = strtok (NULL, "-");
    spl++;
  }
  //printf("ymd: %s\n", ymd);
  //printf("hms: %s\n", hms);
  for (k=0;k<count;k++) {
    printf("%s\n", i_filenames[k]);

    status = nc_open(i_filenames[k], NC_NOWRITE, &ncid);
    handle_error(status);

    status = nc_inq_dimid(ncid, "Radial", &Radial_id);                      
    handle_error(status);
  
    status = nc_inq_dimlen(ncid, Radial_id, &numradials);              
    handle_error(status);
    
    status = nc_inq_dimid(ncid, "Gate", &Gate_id);                          
    handle_error(status);
  
    status = nc_inq_dimlen(ncid, Gate_id, &numgates);                  
    handle_error(status);
  
    Data_In = malloc(numradials*sizeof(float *));       
    Data_In[0]=malloc(numgates*numradials*sizeof(float));
    for(i = 0; i < numradials; i++)
      Data_In[i] = Data_In[0]+i*numgates;
    
    Azimuth_In = (double *)malloc(numradials*sizeof(double));
    Elevation_In = (double *)malloc(numradials*sizeof(double));
    GateWidth_In = (double *)malloc(numradials*sizeof(double));
    StartRange_In = (double *)malloc(numradials*sizeof(double));
    Time_In = (int *)malloc(numradials*sizeof(int));
  
    status = nc_inq_varid(ncid, varname, &Data_id);    
    handle_error(status);
   
    status = nc_inq_varid(ncid, "Azimuth", &Az_id);
    handle_error(status);

    status = nc_inq_varid(ncid, "Elevation", &El_id);
    handle_error(status);

    status = nc_inq_varid(ncid, "GateWidth", &GateWidth_id);
    handle_error(status);

    status = nc_inq_varid(ncid, "Time", &Time_id);
    handle_error(status);

    status = nc_inq_varid(ncid, "StartRange", &Sr_id);
    handle_error(status);
    
    status = nc_get_var_float(ncid, Data_id, &Data_In[0][0]);
    handle_error(status);
    
    for (i=0; i<numradials; i++) {
      for (j=0; j<numgates; j++) {
	Data_In[i][j] = fabs(Data_In[i][j]);
	if (Data_In[i][j] < 11.176)
	  Data_In[i][j] = 0;
      }
    }

    char *radstr;
    if ((radstr=strstr(i_filenames[k], "cleburne.tx"))!=NULL) {
      for (i=0; i<numradials; i++) {
	for (j=0; j<50; j++) {
	  Data_In[i][j] = 0;
	}
      }
    }
    else if ((radstr=strstr(i_filenames[k], "denton.tx"))!=NULL) {
      Flag_In = malloc(numradials*sizeof(int *));
      Flag_In[0]=malloc(numgates*numradials*sizeof(int));
      for(i = 0; i < numradials; i++)
        Flag_In[i] = Flag_In[0]+i*numgates;
      status = nc_inq_varid(ncid, "GateFlags", &GF_id);
      handle_error(status);
      status = nc_get_var_int(ncid, GF_id, &Flag_In[0][0]);
      handle_error(status);
      
      REF_In = malloc(numradials*sizeof(float *));
      REF_In[0]=malloc(numgates*numradials*sizeof(float));
      for(i = 0; i < numradials; i++)
        REF_In[i] = REF_In[0]+i*numgates;
      status = nc_inq_varid(ncid, "Reflectivity", &REF_id);
      handle_error(status);
      status = nc_get_var_float(ncid, REF_id, &REF_In[0][0]);
      handle_error(status);
      
      for (i=0; i<numradials; i++) {
        for (j=0; j<numgates; j++) {
          if ((REF_In[i][j] < 25) || (Flag_In[i][j] == 4)) {
            Data_In[i][j] = 0;
          }
        }
      }
    }
    else if (((radstr=strstr(i_filenames[k], "arlington.tx"))!=NULL) || ((radstr=strstr(i_filenames[k], "midlothian.tx"))!=NULL) || ((radstr=strstr(i_filenames[k], "ftworth.tx"))!=NULL) || ((radstr=strstr(i_filenames[k], "mesquite.tx"))!=NULL)) {
      SNR_In = malloc(numradials*sizeof(float *));
      SNR_In[0]=malloc(numgates*numradials*sizeof(float));
      for(i = 0; i < numradials; i++)
        SNR_In[i] = SNR_In[0]+i*numgates;
      status = nc_inq_varid(ncid, "SignalToNoiseRatio", &SNR_id);
      handle_error(status);
      status = nc_get_var_float(ncid, SNR_id, &SNR_In[0][0]);
      handle_error(status);

      NCP_In = malloc(numradials*sizeof(float *));
      NCP_In[0]=malloc(numgates*numradials*sizeof(float));
      for(i = 0; i < numradials; i++)
        NCP_In[i] = NCP_In[0]+i*numgates;
      status = nc_inq_varid(ncid, "NormalizedCoherentPower", &NCP_id);
      handle_error(status);
      status = nc_get_var_float(ncid, NCP_id, &NCP_In[0][0]);
      handle_error(status);

      RHO_In = malloc(numradials*sizeof(float *));
      RHO_In[0]=malloc(numgates*numradials*sizeof(float));
      for(i = 0; i < numradials; i++)
        RHO_In[i] = RHO_In[0]+i*numgates;
      status = nc_inq_varid(ncid, "CrossPolCorrelation", &RHO_id);
      handle_error(status);
      status = nc_get_var_float(ncid, RHO_id, &RHO_In[0][0]);
      handle_error(status);

      REF_In = malloc(numradials*sizeof(float *));
      REF_In[0]=malloc(numgates*numradials*sizeof(float));
      for(i = 0; i < numradials; i++)
        REF_In[i] = REF_In[0]+i*numgates;
      status = nc_inq_varid(ncid, "Reflectivity", &REF_id);
      handle_error(status);
      status = nc_get_var_float(ncid, REF_id, &REF_In[0][0]);
      handle_error(status);

      for (i=0; i<numradials; i++) {
	for (j=0; j<numgates; j++) {
	  if ((SNR_In[i][j] < 5) || (RHO_In[i][j] < .9) || (NCP_In[i][j] < .5) || (REF_In[i][j] < 25)) {
	    Data_In[i][j] = 0;
	  }
	}
      }
    }

    status = nc_get_var_double(ncid, Az_id, &Azimuth_In[0]);
    handle_error(status);

    status = nc_get_var_double(ncid, El_id, &Elevation_In[0]);
    handle_error(status);
    
    status = nc_get_var_double(ncid, GateWidth_id, &GateWidth_In[0]);
    handle_error(status);
    
    status = nc_get_var_int(ncid, Time_id, &Time_In[0]);
    handle_error(status);
    
    status = nc_get_var_double(ncid, Sr_id, &StartRange_In[0]);
    handle_error(status);
    
    status = nc_get_att_double(ncid, NC_GLOBAL, "Latitude", &radar_llb.lat); 
    handle_error(status);
    
    status = nc_get_att_double(ncid, NC_GLOBAL, "Longitude", &radar_llb.lon);      
    handle_error(status);
  
    status = nc_get_att_double(ncid, NC_GLOBAL, "AntennaBeamwidth", &AntennaBeamwidth);
    //account for a naming convention problem in the EEC radar data
    if (status != NC_NOERR) { 
      status = nc_get_att_double(ncid, NC_GLOBAL, "AntennaBeamWidth", &AntennaBeamwidth);
    }
    handle_error(status);
  
    status = nc_get_att_int(ncid, NC_GLOBAL, "ScanType", &ScanType);
    handle_error(status);
      
    if ((ScanType != 2) && (ScanType != 8)) {
      printf("Not a PPI. Exiting...");
      exit(1);
    }

    if (Time_In[0] > 0) {
      if (Time_In[0] > maxtime)
	maxtime = Time_In[0];
    }

    for (i=0; i<numradials; i++) {
      if (Elevation_In[i] < mintilt)
	continue;
      else if (Elevation_In[i] > maxtilt)
	continue;
 
      double sr = abs(StartRange_In[i]/1000);
      double gw = abs(GateWidth_In[i]/1000);

      for (j=0; j<numgates; j++) {
	
	/*
	if ((radstr=strstr(i_filenames[k], "denton.tx"))!=NULL) {
	  if(Flag_In[i][j] == 4)
	    continue;
	}
	*/
	/*
	if ((radstr=strstr(i_filenames[k], "XADD"))!=NULL) {
	  if((j < 600) && (Data_In[i][j] < 25))
	    continue;
	}
	*/
	double dlat, dlon;
	int latvox, lonvox;
	double dx = (j*gw) + sr;
	double hrise = dx * sin(degToRad(Elevation_In[i]));
	double hcurve = ((dx * cos(degToRad(Elevation_In[i])))/1000) * .12616097;
	
	if ((hcurve + hrise) > maxalt) {
	  //printf("radar: %s j: %d height: %f\n", i_filenames[k], j, (hcurve + hrise));
	  j = numgates - 1;
	}
	if (((hcurve + hrise) < minalt) || ((hcurve + hrise) > maxalt))
	  continue;
	
	//handle front center of beam
	vincenty(radar_llb.lat, radar_llb.lon, dx, Azimuth_In[i], &dest_llb);
	//if the center is outside the grid, we toss
	if ((dest_llb.lat > nlat) || (dest_llb.lat < slat) || (dest_llb.lon < wlon) || (dest_llb.lon > elon))
	  continue;

	dlat = nlat - dest_llb.lat;
	dlon = dest_llb.lon - wlon;
	latvox = (int) (dlat / latspac);
	lonvox = (int) (dlon / lonspac);

	if ((Data_In[i][j] > minval) && (Data_In[i][j] < maxval)) {
	  //check routine
	  if (routine == 0) {
	    grid[latvox][lonvox] = grid[latvox][lonvox] + Data_In[i][j];
	    sample[latvox][lonvox]++;
	  }
	  else if (routine == 1) {
	    if (Data_In[i][j] > grid[latvox][lonvox])
	      grid[latvox][lonvox] = Data_In[i][j];
	  }
	}

	if (truemap == 1) {
	  double tmpaz, tmpdx, tmprise, tmpcurve;
	  int tmplatvox, tmplonvox, tmlatvox, tmlonvox; 
	
	  //handle back center of beam
	  tmpdx = dx + (GateWidth_In[i]/1000);
	  tmprise = tmpdx * sin(degToRad(Elevation_In[i]));
	  tmpcurve = ((tmpdx * cos(degToRad(Elevation_In[i])))/1000) * .12616097;
	  if (((tmpcurve + tmprise) > minalt) && ((tmpcurve + tmprise) < maxalt)) {
	    vincenty(radar_llb.lat, radar_llb.lon, tmpdx, Azimuth_In[i], &dest_llb);
	    if ((dest_llb.lat < nlat) && (dest_llb.lat > slat) && (dest_llb.lon > wlon) && (dest_llb.lon < elon)) {
	      dlat = nlat - dest_llb.lat;
	      dlon = dest_llb.lon - wlon;
	      tmplatvox = (int) (dlat / latspac);
	      tmplonvox = (int) (dlon / lonspac);
	      //if the front of the beam and the back of the beam are the same voxel we avoid double counting
	      if ((tmplatvox != latvox) || (tmplonvox != lonvox)) {
		if ((Data_In[i][j] > minval) && (Data_In[i][j] < maxval)) {
		  if (routine == 0) {
		    grid[tmplatvox][tmplonvox] = grid[tmplatvox][tmplonvox] + Data_In[i][j];
		    sample[tmplatvox][tmplonvox]++;
		  }
		  else if (routine == 1) {
		    if (Data_In[i][j] > grid[tmplatvox][tmplonvox])
		      grid[tmplatvox][tmplonvox] = Data_In[i][j];
		  }
		}
	      }
	    }
	  }
	  
	  //handle ccw front side of beam
	  tmpaz = Azimuth_In[i] - (.5*AntennaBeamwidth);
	  if (tmpaz < 0)
	    tmpaz = tmpaz + 360;
	  vincenty(radar_llb.lat, radar_llb.lon, dx, tmpaz, &dest_llb);
	  if ((dest_llb.lat < nlat) && (dest_llb.lat > slat) && (dest_llb.lon > wlon) && (dest_llb.lon < elon)) {
	    dlat = nlat - dest_llb.lat;
	    dlon = dest_llb.lon - wlon;
	    tmplatvox = (int) (dlat / latspac);
	    tmplonvox = (int) (dlon / lonspac);
	    //if the center of the beam and the edge of the beam are the same voxel we move on to avoid double counting
	    if ((tmplatvox != latvox) || (tmplonvox != lonvox)) {
	      if ((Data_In[i][j] > minval) && (Data_In[i][j] < maxval)) {
		//check routine
		if (routine == 0) {
		  grid[tmplatvox][tmplonvox] = grid[tmplatvox][tmplonvox] + Data_In[i][j];
		  sample[tmplatvox][tmplonvox]++;
		}
		else if (routine == 1) {
		  if (Data_In[i][j] > grid[tmplatvox][tmplonvox])
		    grid[tmplatvox][tmplonvox] = Data_In[i][j];
		}
	      }
	    }
	  }
	  
	  //handle ccw back side of beam
	  tmpaz = Azimuth_In[i] - (.5*AntennaBeamwidth);
	  if (tmpaz < 0)
	    tmpaz = tmpaz + 360;
	  if (((tmpcurve + tmprise) > minalt) && ((tmpcurve + tmprise) < maxalt)) {
	    vincenty(radar_llb.lat, radar_llb.lon, tmpdx, tmpaz, &dest_llb);
	    if ((dest_llb.lat < nlat) && (dest_llb.lat > slat) && (dest_llb.lon > wlon) && (dest_llb.lon < elon)) {
	      dlat = nlat - dest_llb.lat;
	      dlon = dest_llb.lon - wlon;
	      tmlatvox = (int) (dlat / latspac);
	      tmlonvox = (int) (dlon / lonspac);
	      
	      if ((tmlatvox != tmplatvox) || (tmlonvox != tmplonvox)) {
		if ((Data_In[i][j] > minval) && (Data_In[i][j] < maxval)) {
		  if (routine == 0) {
		    grid[tmlatvox][tmlonvox] = grid[tmlatvox][tmlonvox] + Data_In[i][j];
		    sample[tmlatvox][tmlonvox]++;
		  }
		  else if (routine == 1) {
		    if (Data_In[i][j] > grid[tmlatvox][tmlonvox])
		      grid[tmlatvox][tmlonvox] = Data_In[i][j];
		  }
		}
	      }
	    }
	  }
	  
	  //handle cw front side of beam
	  tmpaz = Azimuth_In[i] + (.5*AntennaBeamwidth);
	  if (tmpaz > 359.999999)
	    tmpaz = tmpaz - 360;
	  vincenty(radar_llb.lat, radar_llb.lon, dx, tmpaz, &dest_llb);
	  if ((dest_llb.lat < nlat) && (dest_llb.lat > slat) && (dest_llb.lon > wlon) && (dest_llb.lon < elon)) {
	    dlat = nlat - dest_llb.lat;
	    dlon = dest_llb.lon - wlon;
	    tmplatvox = (int) (dlat / latspac);
	    tmplonvox = (int) (dlon / lonspac);
	    
	    if ((tmplatvox != latvox) || (tmplonvox != lonvox)) {
	      if ((Data_In[i][j] > minval) && (Data_In[i][j] < maxval)) {
		if (routine == 0) {
		  grid[tmplatvox][tmplonvox] = grid[tmplatvox][tmplonvox] + Data_In[i][j];
		  sample[tmplatvox][tmplonvox]++;
		}
		else if (routine == 1) {
		  if (Data_In[i][j] > grid[tmplatvox][tmplonvox])
		    grid[tmplatvox][tmplonvox] = Data_In[i][j];
		}
	      }
	    }
	  }
	  
	  //handle cw back of the beam
	  tmpaz = Azimuth_In[i] + (.5*AntennaBeamwidth);
	  if (tmpaz > 359.999999)
	    tmpaz = tmpaz - 360;
	  if (((tmpcurve + tmprise) > minalt) && ((tmpcurve + tmprise) < maxalt)) {
	    vincenty(radar_llb.lat, radar_llb.lon, tmpdx, tmpaz, &dest_llb);
	    if ((dest_llb.lat < nlat) && (dest_llb.lat > slat) && (dest_llb.lon > wlon) && (dest_llb.lon < elon)) {
	      dlat = nlat - dest_llb.lat;
	      dlon = dest_llb.lon - wlon;
	      tmlatvox = (int) (dlat / latspac);
	      tmlonvox = (int) (dlon / lonspac);
	      
	      if ((tmlatvox != tmplatvox) || (tmlonvox != tmplonvox)) {
		if ((Data_In[i][j] > minval) && (Data_In[i][j] < maxval)) {
		  if (routine == 0) {
		    grid[tmlatvox][tmlonvox] = grid[tmlatvox][tmlonvox] + Data_In[i][j];
		    sample[tmlatvox][tmlonvox]++;
		  }
		  else if (routine == 1) {
		    if (Data_In[i][j] > grid[tmlatvox][tmlonvox])
		      grid[tmlatvox][tmlonvox] = Data_In[i][j];
		  }
		}
	      }
	    }
	  }
	  
	  //handle in between areas
	  
          tmpaz = Azimuth_In[i] - (.25*AntennaBeamwidth);
          if (tmpaz < 0)
            tmpaz = tmpaz + 360;
          vincenty(radar_llb.lat, radar_llb.lon, dx, tmpaz, &dest_llb);
          if ((dest_llb.lat < nlat) && (dest_llb.lat > slat) && (dest_llb.lon > wlon) && (dest_llb.lon < elon)) {
            dlat = nlat - dest_llb.lat;
            dlon = dest_llb.lon - wlon;
            tmplatvox = (int) (dlat / latspac);
            tmplonvox = (int) (dlon / lonspac);

            if ((tmplatvox != latvox) || (tmplonvox != lonvox)) {
              if ((Data_In[i][j] > minval) && (Data_In[i][j] < maxval)) {
                //check routine                                                                                                                                                                
                if (routine == 0) {
                  grid[tmplatvox][tmplonvox] = grid[tmplatvox][tmplonvox] + Data_In[i][j];
                  sample[tmplatvox][tmplonvox]++;
                }
                else if (routine == 1) {
                  if (Data_In[i][j] > grid[tmplatvox][tmplonvox])
                    grid[tmplatvox][tmplonvox] = Data_In[i][j];
                }
              }
            }
          }

	  tmpaz = Azimuth_In[i] - (.25*AntennaBeamwidth);
          if (tmpaz < 0)
            tmpaz = tmpaz + 360;
          if (((tmpcurve + tmprise) > minalt) && ((tmpcurve + tmprise) < maxalt)) {
            vincenty(radar_llb.lat, radar_llb.lon, tmpdx, tmpaz, &dest_llb);
            if ((dest_llb.lat < nlat) && (dest_llb.lat > slat) && (dest_llb.lon > wlon) && (dest_llb.lon < elon)) {
              dlat = nlat - dest_llb.lat;
              dlon = dest_llb.lon - wlon;
              tmlatvox = (int) (dlat / latspac);
              tmlonvox = (int) (dlon / lonspac);

	      if ((tmlatvox != tmplatvox) || (tmlonvox != tmplonvox)) {
                if ((Data_In[i][j] > minval) && (Data_In[i][j] < maxval)) {
                  if (routine == 0) {
                    grid[tmlatvox][tmlonvox] = grid[tmlatvox][tmlonvox] + Data_In[i][j];
                    sample[tmlatvox][tmlonvox]++;
                  }
                  else if (routine == 1) {
                    if (Data_In[i][j] > grid[tmlatvox][tmlonvox])
                      grid[tmlatvox][tmlonvox] = Data_In[i][j];
                  }
                }
              }
            }
          }

	  tmpaz = Azimuth_In[i] + (.25*AntennaBeamwidth);
          if (tmpaz > 359.999999)
            tmpaz = tmpaz - 360;
          vincenty(radar_llb.lat, radar_llb.lon, dx, tmpaz, &dest_llb);
          if ((dest_llb.lat < nlat) && (dest_llb.lat > slat) && (dest_llb.lon > wlon) && (dest_llb.lon < elon)) {
            dlat = nlat - dest_llb.lat;
            dlon = dest_llb.lon - wlon;
            tmplatvox = (int) (dlat / latspac);
            tmplonvox = (int) (dlon / lonspac);

            if ((tmplatvox != latvox) || (tmplonvox != lonvox)) {
              if ((Data_In[i][j] > minval) && (Data_In[i][j] < maxval)) {
                if (routine == 0) {
                  grid[tmplatvox][tmplonvox] = grid[tmplatvox][tmplonvox] + Data_In[i][j];
                  sample[tmplatvox][tmplonvox]++;
                }
                else if (routine == 1) {
                  if (Data_In[i][j] > grid[tmplatvox][tmplonvox])
                    grid[tmplatvox][tmplonvox] = Data_In[i][j];
                }
              }
            }
          }

	  tmpaz = Azimuth_In[i] + (.25*AntennaBeamwidth);
          if (tmpaz > 359.999999)
            tmpaz = tmpaz - 360;
          if (((tmpcurve + tmprise) > minalt) && ((tmpcurve + tmprise) < maxalt)) {
            vincenty(radar_llb.lat, radar_llb.lon, tmpdx, tmpaz, &dest_llb);
            if ((dest_llb.lat < nlat) && (dest_llb.lat > slat) && (dest_llb.lon > wlon) && (dest_llb.lon < elon)) {
              dlat = nlat - dest_llb.lat;
              dlon = dest_llb.lon - wlon;
              tmlatvox = (int) (dlat / latspac);
              tmlonvox = (int) (dlon / lonspac);

              if ((tmlatvox != tmplatvox) || (tmlonvox != tmplonvox)) {
                if ((Data_In[i][j] > minval) && (Data_In[i][j] < maxval)) {
                  if (routine == 0) {
                    grid[tmlatvox][tmlonvox] = grid[tmlatvox][tmlonvox] + Data_In[i][j];
                    sample[tmlatvox][tmlonvox]++;
                  }
                  else if (routine == 1) {
                    if (Data_In[i][j] > grid[tmlatvox][tmlonvox])
                      grid[tmlatvox][tmlonvox] = Data_In[i][j];
                  }
                }
              }
            }
          }
	}
      }
    }
    
    
    status = nc_close(ncid);
    handle_error(status);
    
    free(Data_In);
    free(Azimuth_In);
    free(Elevation_In);
    free(GateWidth_In);
    free(Time_In); 
    free(StartRange_In);
    if ((radstr=strstr(i_filenames[k], "denton.tx"))!=NULL) {
      free(Flag_In);
      free(REF_In);
    }
    else if (((radstr=strstr(i_filenames[k], "arlington.tx"))!=NULL) || ((radstr=strstr(i_filenames[k], "midlothian.tx"))!=NULL) || ((radstr=strstr(i_filenames[k], "ftworth.tx"))!=NULL) || ((radstr=strstr(i_filenames[k], "mesquite.tx"))!=NULL)) {
      free(SNR_In);
      free(RHO_In);
      free(NCP_In);
      free(REF_In);
    }
  }

  if (routine == 0) {
    for (x = 0; x<numlats; x++) {
      for (y = 0; y<numlons; y++) {
	if (sample[x][y] > 0) {
	  grid[x][y] = grid[x][y] / (float)sample[x][y];
	  //printf("x %d y %d grid %f sample %d\n", x, y, grid[x][y], sample[x][y]);
	}
      }
    }
  }

  if (smoothing == 1) {
    float tmpgrid[numlats][numlons];
    for (x = 0; x<numlats; x++) 
      for (y = 0; y<numlons; y++) 
	tmpgrid[x][y] = grid[x][y];
    
    int z = 1;
    int e;
    for (e = 0; e < 3; e++) {
      for (x = smoothvoxels; x<(numlats-smoothvoxels); x++) {
	for (y = smoothvoxels; y<(numlons-smoothvoxels); y++) {	
	  int shepard_k = smoothvoxels * 16;
	  float shepard = 0;
	  int shepard_denom_k = smoothvoxels * 16;
	  float shepard_denom_tot = 0;
	  float shepard_tot = 0;
	  float crossdx = sqrt((latspac * latspac) + (lonspac * lonspac));
	  for (z = 1; z<(smoothvoxels+1); z++) 
	    shepard_denom_tot = shepard_denom_tot + (2*(1/((z*latspac)*(z*latspac)))) + (2*(1/((z*lonspac)*(z*lonspac)))) + (4*(1/((z*crossdx)*(z*crossdx))));
	  float shepard_denom = shepard_denom_tot/(float)shepard_denom_k;
	  for (z=1; z<smoothvoxels+1; z++)
	    shepard_tot = shepard_tot + ((((1/((z*latspac)*(z*latspac))) * grid[x-z][y]) + ((1/((z*latspac)*(z*latspac))) * grid[x+z][y]) + ((1/((z*lonspac)*(z*lonspac))) * grid[x][y-z]) + ((1/((z*lonspac)*(z*lonspac))) * grid[x][y+z])) + ((1/((z*crossdx)*(z*crossdx))) * grid[x-z][y-z]) + ((1/((z*crossdx)*(z*crossdx))) * grid[x+z][y-z]) + ((1/((z*crossdx)*(z*crossdx))) * grid[x-z][y+z]) + ((1/((z*crossdx)*(z*crossdx))) * grid[x+z][y+z]));
	  
	  shepard_tot = shepard_tot/shepard_denom;
	  if (shepard_k > 0) {
	    shepard = shepard_tot/shepard_k;
	    if (e == 0) {
	      if (tmpgrid[x][y] > 0) {
		grid[x][y] = ((3*tmpgrid[x][y]) + shepard)/4;
	      }
	      else { 
		grid[x][y] = shepard;
	      }
	    }
	    else if (e == 1) {
	      if (tmpgrid[x][y] > 0) {
                grid[x][y] = ((tmpgrid[x][y]) + shepard)/2;
              }
              else {
		grid[x][y] = shepard;
              }
            }
	    else {
	      grid[x][y] = shepard;
	    }
	  }
	  else {
	    grid[x][y] = 0;
	  }
	}
      }
    }
  }

  //tm_t_0 = (time_t)maxtime;
  //strftime(starttime,16,"%Y%m%d-%H%M%S",gmtime_r(&tm_t_0,&broketime));
  sprintf(starttime, "%s-%s", ymd, hms);
  //now we write out a grid file  
  char o_filename[256];
  char varname_out[64];
  if (strlen(varname) < 52)
    sprintf(varname_out, "Max%s", varname);
  else
    sprintf(varname_out, "MaxData");
  
  if (strlen(varname) < 52) 
    sprintf(o_filename, "%s/%s_%s.netcdf", ncdir, varname_out, starttime);
  else 
    sprintf(o_filename, "%s/maxVel_%s.netcdf", ncdir, starttime);
  
  
  static char datatype[] = "LatLonGrid";
  static int fractionaltime[] = {0};
  float tmplat[] = {(float)nlat};
  float tmplon[] = {(float)wlon};
  float tmpalt[] = {(float)minalt};
  int tmptm[] = {maxtime};
  float tmplatspac[] = {(float)latspac};
  float tmplonspac[] = {(float)lonspac};

  int ncid_out, grid_id, lat_id, lon_id, grid_dimids[2];

  status = nc_create(o_filename, NC_CLOBBER, &ncid_out);
  handle_error(status);

  status = nc_def_dim(ncid_out, "Lat", numlats, &lat_id);
  handle_error(status);
  
  status = nc_def_dim(ncid_out, "Lon", numlons, &lon_id);
  handle_error(status);

  grid_dimids[0] = lat_id;
  grid_dimids[1] = lon_id;

  status = nc_def_var(ncid_out, varname_out, NC_FLOAT, 2, grid_dimids, &grid_id);
  handle_error(status);

  status = nc_put_att_text(ncid_out, grid_id, "Units", strlen(units), units);
  handle_error(status);

  status = nc_put_att_text(ncid_out, NC_GLOBAL, "TypeName", strlen(varname_out), varname_out);
  handle_error(status);

  status = nc_put_att_text(ncid_out, NC_GLOBAL, "DataType", strlen(datatype), datatype);
  handle_error(status);

  status = nc_put_att_int(ncid_out, NC_GLOBAL, "FractionalTime", NC_INT, 1, fractionaltime);
  handle_error(status);

  status = nc_put_att_float(ncid_out, NC_GLOBAL, "Latitude", NC_FLOAT, 1, tmplat);
  handle_error(status);

  status = nc_put_att_float(ncid_out, NC_GLOBAL, "Longitude", NC_FLOAT, 1, tmplon);
  handle_error(status);
  
  status = nc_put_att_float(ncid_out, NC_GLOBAL, "Height", NC_FLOAT, 1, tmpalt);
  handle_error(status);

  status = nc_put_att_int(ncid_out, NC_GLOBAL, "Time", NC_INT, 1, tmptm);
  handle_error(status);

  status = nc_put_att_float(ncid_out, NC_GLOBAL, "LatGridSpacing", NC_FLOAT, 1, tmplatspac);
  handle_error(status);
  
  status = nc_put_att_float(ncid_out, NC_GLOBAL, "LonGridSpacing", NC_FLOAT, 1, tmplonspac);
  handle_error(status);

  status = nc_enddef(ncid_out);
  handle_error(status);

  status = nc_put_var_float(ncid_out, grid_id, &grid[0][0]);
  handle_error(status);

  status = nc_close(ncid_out);
  handle_error(status);
   
  config_destroy(&UM_config);
  
  return 0;
}

