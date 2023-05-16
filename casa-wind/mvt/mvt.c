/* $Id: mvt.c,v 0.1 2018-05-18 14:08:00 elyons Exp $ */
/* Copyright 2018 University of Massachusetts Amherst all rights reserved*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <time.h>
#include <netcdf.h>
#include <math.h>
#include <jansson.h>
#include <libconfig.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "threshold_functions.h"
	
int mvt(char i_filename[]) {
  // something to keep track of our contours later on
  int num_valid_contours = 0; 
  /*locate and verify the config file*/
  char *mvt_home = getenv("MVTHOME");
  config_t mvt_config;
  config_setting_t *acs_params;
  config_setting_t *alg_params;
  char config_fn[128];
  sprintf(config_fn, "%s/mvt_config.txt", mvt_home);
  printf("config_fn is %s\n", config_fn);

  config_init(&mvt_config);
  if (!config_read_file(&mvt_config, config_fn)) {
    printf("config file not found or contains errors.  Should be in $MVTHOME directory and called mvt_config.txt. exiting...\n");
    config_destroy(&mvt_config);
    return -1;
  }
  
  /* Check whether this will report to the MCC via sockets */
  int mcc_output = 0;
  struct sockaddr_in mcc_addr;
  int mcc_sock;
  long int mcc_port;
  const char *mcc_host;
  int connected = 0;
  
  if (config_lookup_bool(&mvt_config, "mcc_output", &mcc_output))
    printf("MCC output: %d\n", mcc_output);
  else
    printf("No mcc_output definition in mvt_config.  MCC output disabled.\n");

  if (mcc_output == 1) {
    if (!config_lookup_int(&mvt_config, "mcc_port", &mcc_port)) {
      printf("No mcc_port definition in mvt_config. Exiting...\n");
      exit(1);
    }
    
    if (!config_lookup_string(&mvt_config, "mcc_host", &mcc_host)) {
      printf("No mcc_host definition in mvt_config. Exiting...\n");
      exit(1);
    }
  }
  
  /* Check whether this will report to the ACS */
  int acs_output = 0;
  int num_acs_hosts = 0;
  int num_acs_ports = 0;
  const config_setting_t *acs_host_arr;
  const config_setting_t *acs_port_arr;
  
  if (config_lookup_bool(&mvt_config, "acs_output", &acs_output))
    printf("ACS output: %d\n", acs_output);
  else
    printf("No acs_output definition in mvt_config.  ACS output disabled.\n");
  
  if (acs_output == 1) {
    acs_params = config_lookup(&mvt_config, "acs_params");

    if (acs_params != NULL) {
      acs_host_arr = config_setting_get_member(acs_params, "acs_host");
      if (acs_host_arr != NULL) {
	num_acs_hosts = config_setting_length(acs_host_arr);
	printf("num hosts %d\n", num_acs_hosts);
      }
      else {
	printf("ACS_output is true, but no ACS hosts defined.  Exiting.\n");
	exit(-1);
      }
      acs_port_arr = config_setting_get_member(acs_params, "acs_port");
      if (acs_port_arr != NULL) {
	num_acs_ports = config_setting_length(acs_port_arr);
	
      }
      else {
	printf("ACS_output is true, but no ACS ports defined.  Exiting.\n");
	exit(-1);
      }
      if (num_acs_ports != num_acs_hosts) {
	printf("Number of ACS ports does not match number of ACS hosts.  Please ensure you have one port entry for each host in the configuration file.  Exiting.\n");
	exit(-1);
      }  
    }
    else {
      printf("No acs_params in config file.  Exiting.\n");
      exit(-1);
    }
  }
  struct sockaddr_in acs_addr[num_acs_hosts];
  int acs_sock[num_acs_hosts];
  long int acs_port[num_acs_ports];
  const char *acs_host[num_acs_hosts];
  int sock_connected[num_acs_hosts];
  int gh;
  for (gh = 0; gh < num_acs_hosts; gh++) {
    sock_connected[gh] = 0;
    acs_port[gh] = config_setting_get_int_elem(acs_port_arr, gh);
    acs_host[gh] = config_setting_get_string_elem(acs_host_arr, gh);
    if ((!acs_port[gh]) || (!acs_host[gh])) {
      printf("ACS port or host not properly defined.  Please check config file.  Exiting...\n");
      exit(-1);
    }
    else 
      printf("acs connect %d port: %d host: %s\n", gh, acs_port[gh], acs_host[gh]);

    if ((acs_sock[gh]=socket(AF_INET, SOCK_STREAM, 0))==-1){
      perror("socket");
      exit(-1);
    }

    struct timeval timeout = {2, 0};

    if (setsockopt(acs_sock[gh], SOL_SOCKET, SO_SNDTIMEO, (struct timeval *)&timeout,sizeof(struct timeval)) < 0) {
      perror("setsockopt failed\n");
    }

    memset((char *) &acs_addr[gh],'\0',sizeof(acs_addr[gh]));
    acs_addr[gh].sin_family = AF_INET;
    acs_addr[gh].sin_port = htons(acs_port[gh]);
    acs_addr[gh].sin_addr.s_addr = inet_addr(acs_host[gh]);

    if (connect(acs_sock[gh],(struct sockaddr *) &acs_addr[gh],
		sizeof (struct sockaddr_in))<0){
      perror("connect");
      //exit(-1);                                                                                                                                                                      
    }
    else {
      sock_connected[gh] = 1;
    }
  }
  /* Check whether this will record output in a file */
  /*
  int file_output = 0;
  FILE *f1;
  const char *output_fn;

  if (config_lookup_bool(&mvt_config, "file_output", &file_output))
    printf("file output: %d\n", file_output);
  else
    printf("No file_output definition in mvt_config.  File output disabled.\n");

  if (file_output == 1) {
    if(!config_lookup_string(&mvt_config, "file_name", &output_fn)) {
      printf("No file_name definition in mvt_config. exiting...\n");
      exit(1);
    }
    f1 = fopen(output_fn, "a");
    }
  */

  /* Check whether netcdf file output will be generated */
  int netcdf_output = 0;

  if (config_lookup_bool(&mvt_config, "netcdf_output", &netcdf_output))
    printf("netcdf output: %d\n", netcdf_output);
  else
    printf("No netcdf_output definition in mvt_config.  Netcdf output disabled.\n");

  /* Check whether an xml data file will be generated */
  int xml_output = 0;
  int xmlsize;
  int y;
  char *xmlbuf;
  FILE *xmlfile;
  char xmlname[128];
  struct latLon *ll_array;
  ll_array = (struct latLon*)malloc(sizeof(struct latLon));
  int edgecount = 0;
  const char *xmldir;
  xmlDocPtr doc;
  xmlNodePtr root,mvt_edges, contour_number, contour_points;

  if (config_lookup_bool(&mvt_config, "xml_output", &xml_output)){
    printf("XML output: %d\n", xml_output);
  }
  else
    printf("No xml_output definition in mvt_config.  XML output disabled.\n");

  if (config_lookup_string(&mvt_config, "xml_directory", &xmldir)){
    printf("XML directory: %s\n", xmldir);
  }
  else {
    printf("No xml_directory definition in mvt_config.  xml output disabled.\n");
    xml_output = 0;
  }
  
  /* Check whether a json data file will be generated */
  int json_output = 0;
  char *jsonbuf;
  //FILE *jsonfile;
  char jsonname[128];
  const char *jsondir;
  json_t *json_coords_arr;//  = json_array();
  json_t *jsonroot = json_object();
  json_t *features_json_arr = json_array();
  json_object_set_new( jsonroot, "type", json_string( "FeatureCollection" ) );
  json_object_set_new( jsonroot, "features", features_json_arr );
    
  if (config_lookup_bool(&mvt_config, "json_output", &json_output)){
    printf("JSON output: %d\n", json_output);
  }
  else 
    printf("No json_output definition in mvt_config.  JSON output disabled.\n");

  if (config_lookup_string(&mvt_config, "json_directory", &jsondir)){
    printf("JSON directory: %s\n", jsondir);
  }
  else {
    printf("No json_directory definition in mvt_config.  JSON output disabled.\n");
    json_output = 0;
  }

  /* Check whether a kml data file will be generated */
  /*
  int kml_output = 0;
  int kmlsize;
  char *kmlbuf;
  FILE *kmlfile;
  char kmlname[128];
  const char *kmldir;
  
  if (config_lookup_bool(&mvt_config, "kml_output", &kml_output)){
    printf("KML output: %d\n", kml_output);
  }
  else
    printf("No kml_output definition in mvt_config.  KML output disabled.\n");

  if (config_lookup_string(&mvt_config, "kml_directory", &kmldir)){
    printf("KML directory: %s\n", kmldir);
  }
  else {
    printf("No kml_directory definition in mvt_config.  kml output disabled.\n");
    kml_output = 0;
  }
  */
  
  /*read the alg_params */
  int num_contour_levels = 0;
  alg_params = config_lookup(&mvt_config, "alg_params");
  int min_contour_points;
  int min_contour_points_present;
  const config_setting_t *contour_levels;
  if (alg_params != NULL) {
    contour_levels = config_setting_get_member(alg_params, "contour_levels");
    if (contour_levels != NULL) {
      num_contour_levels = config_setting_length(contour_levels);
    }
    else {
      printf("No contour_levels defined.  Exiting.\n");
      exit(-1);
    }
    min_contour_points_present = config_setting_lookup_int(alg_params, "min_contour_points", &min_contour_points);
    if (min_contour_points_present == CONFIG_FALSE) {
      printf("min_contour_points not defined. No size threshold will be used.");
      min_contour_points = 0;
    }
    else {
      printf("Minimum contour points: %d\n", min_contour_points);
    }
  }
  else {
    printf("No alg_params in config file.  Exiting.\n");
    exit(-1);
  }
  float contours[num_contour_levels];
  int cl;
  for (cl = 0; cl < num_contour_levels; cl++) {
    contours[cl] = config_setting_get_float_elem(contour_levels, cl);
    printf("Contour level %d: %f\n", cl, contours[cl]);
  }

  /* Variables for the netcdf read */
  float **data_in;

  struct index_pair ***indexe;
  struct index_pair ***indexes;
  struct index_pair ***indexe1;
  struct index_pair ***indexes1;

  float *lat_in;
  float *lon_in;
  
  float n_lat, w_lon, lat_spac, lon_spac;
  int *lonaxis, *lataxis;
  int status,ncid,latid,lonid,data_id,h,i,j;
  char varname[128];
  size_t num_lats, num_lons;
  
  status = nc_open(i_filename, NC_NOWRITE, &ncid);
  handle_error(status);
  if (status != NC_NOERR)
    exit(-1);
  
  status = nc_inq_dimid(ncid, "Lat", &latid);    
  handle_error(status);

  status = nc_inq_dimlen(ncid, latid, &num_lats);
  handle_error(status);
  
  status = nc_inq_dimid(ncid, "Lon", &lonid);    
  handle_error(status);

  status = nc_inq_dimlen(ncid, lonid, &num_lons);
  handle_error(status);
  
  lat_in = malloc(num_lats*sizeof(float));
  lon_in = malloc(num_lons*sizeof(float));

  data_in = malloc(num_lons*sizeof(float *));
  data_in[0] = malloc(num_lons*num_lats*sizeof(float));
  for(i=0; i < num_lons; ++i)
    data_in[i] = data_in[0]+i*num_lats;

  indexe = (struct index_pair ***)malloc(sizeof(struct index_pair **) * num_contour_levels);
  for (i=0; i < num_contour_levels; i++) {
    indexe[i] = (struct index_pair **)malloc(sizeof(struct index_pair *) * num_lons);
    for (j=0; j < num_lons; ++j)
      indexe[i][j] = (struct index_pair *)malloc(sizeof(struct index_pair) * num_lats);
  }
  
  indexe1 = (struct index_pair ***)malloc(sizeof(struct index_pair **) *num_contour_levels);
  for (i=0; i <num_contour_levels; i++) {
    indexe1[i] =(struct index_pair **)malloc(sizeof(struct index_pair *) * num_lons);
    for(j=0; j< num_lons; ++j)
      indexe1[i][j] = (struct index_pair*)malloc(sizeof(struct index_pair) * num_lats);
  }

  indexes = (struct index_pair ***)malloc(sizeof(struct index_pair **) *num_contour_levels);
  for (i=0; i <num_contour_levels; i++) {
    indexes[i] =(struct index_pair **)malloc(sizeof(struct index_pair *) * num_lons);
    for(j=0; j< num_lons; ++j)
      indexes[i][j] = (struct index_pair*)malloc(sizeof(struct index_pair) * num_lats);
  }

  indexes1 = (struct index_pair ***)malloc(sizeof(struct index_pair **) *num_contour_levels);
  for (i=0; i <num_contour_levels; i++) {
    indexes1[i] =(struct index_pair **)malloc(sizeof(struct index_pair *) * num_lons);
    for(j=0; j< num_lons; ++j)
      indexes1[i][j] = (struct index_pair*)malloc(sizeof(struct index_pair) * num_lats);
  }

  status = nc_inq_varname(ncid, 0, varname);
  handle_error(status);

  status = nc_inq_varid(ncid,varname,&data_id);
  handle_error(status);

  status = nc_get_var_float(ncid, data_id, &data_in[0][0]);
  handle_error(status);
    
  status = nc_get_att_float(ncid, NC_GLOBAL, "Latitude", &n_lat);
  handle_error(status);

  status = nc_get_att_float(ncid, NC_GLOBAL, "Longitude", &w_lon);
  handle_error(status);

  status = nc_get_att_float(ncid, NC_GLOBAL, "LatGridSpacing", &lat_spac);
  handle_error(status);

  status = nc_get_att_float(ncid, NC_GLOBAL, "LonGridSpacing", &lon_spac);
  handle_error(status);

  int v;
  for (v=1; v < num_lons+1; v++) 
    lon_in[v-1] = (w_lon + (lon_spac*(v-1)));
  
  for (v=1; v < num_lats+1; v++)
    lat_in[v-1] = (n_lat - (lat_spac*(v-1))); 
  //for (v=0; v < num_lats; v++)
  //  lat_in[v] = (n_lat - (lat_spac*((num_lats - v)-1)));
  
  
  size_t filename_len;
  filename_len = strlen(i_filename);
  char starttime[16];
  char yyyy[4];
  char mon[3];
  char dd[3];
  char hh[3];
  char mm[3];
  char ss[3];
  int ymd_start = (int)filename_len - 22;
  memcpy(yyyy, &i_filename[ymd_start], 4);
  yyyy[4] = '\0';
  memcpy(mon, &i_filename[ymd_start+4], 2);
  mon[2] = '\0';
  memcpy(dd, &i_filename[ymd_start+6], 2);
  dd[2] = '\0';
  memcpy(hh, &i_filename[ymd_start+9], 2);
  hh[2] = '\0';
  memcpy(mm, &i_filename[ymd_start+11], 2);
  mm[2] = '\0';
  memcpy(ss, &i_filename[ymd_start+13], 2);
  ss[2] = '\0';
  starttime[15] = '\0';
  sprintf(starttime, "%s%s%s-%s%s%s", yyyy, mon, dd, hh, mm, ss);
  
  char alt_starttime[32];
  char mcctime[16];
  
  sprintf(mcctime,"%s%s%s%s%s00", yyyy, mon, dd, hh, mm);
  sprintf(alt_starttime, "%s-%s-%sT%s:%s:00Z", yyyy, mon, dd, hh, mm);
  
  char name_string[128];
  sprintf(name_string, "%s", varname);
  
 
  if (mcc_output == 1) {
    if ((mcc_sock=socket(AF_INET, SOCK_STREAM, 0))==-1){
      perror("socket");
    }
    else {
      memset((char *) &mcc_addr,'\0',sizeof(mcc_addr));
      mcc_addr.sin_family = AF_INET;
      mcc_addr.sin_port = htons(mcc_port);
      mcc_addr.sin_addr.s_addr = inet_addr(mcc_host);
      if (connect(mcc_sock,(struct sockaddr *) &mcc_addr, sizeof (struct sockaddr_in))<0) {
	perror("connect");
      }
      else
	connected = 1;
    }
  }
  
  if (json_output == 1) {
    sprintf(jsonname, "%s/mvt_%s_%s.geojson", jsondir, name_string, starttime);
  }

  if (xml_output == 1) {
    sprintf(xmlname, "%s/mvt_%s_%s.xml", xmldir, name_string, starttime);
    doc=xmlNewDoc(BAD_CAST "1.0");
    root=xmlNewNode(NULL,BAD_CAST "mvt_edges");
    xmlDocSetRootElement(doc,root);
    xmlNewChild(root,NULL,BAD_CAST "filename",BAD_CAST i_filename);
  }
  
  /*
  if (kml_output == 1)
    sprintf(kmlname, "%s/mvt_%s_%s.kml", kmldir, name_string, starttime);
  */
  
  lataxis = malloc(num_lats*sizeof(int));
  lonaxis = malloc(num_lons*sizeof(int));
  
  for(i=0; i<num_lats; ++i)
    lataxis[i]=i;

  for(i=0; i<num_lons; ++i) 
    lonaxis[i]=i;
  
  //initialize arrays to -1... 
  for (h=0; h<num_contour_levels; ++h) {
    for (i=0; i<num_lons; ++i) {
      for (j=0; j<num_lats; ++j) {
      indexe[h][i][j].xind = -1;
      indexe[h][i][j].yind = -1;
      indexes[h][i][j].xind = -1;
      indexes[h][i][j].yind = -1;
      indexe1[h][i][j].xind = -1;
      indexe1[h][i][j].yind = -1;
      indexes1[h][i][j].xind = -1;
      indexes1[h][i][j].yind = -1;
      }
    }
  }
  //printf("numlons %d numlats %d\n", num_lons, num_lats);
   
  //ncontour set to 1, but could be set to do multiple contours eventually... candidate for config file
  //contours array is also going to have a single entry for now, but again, could have more
  
  //function call to generate conrec contours.
  Contour(data_in, 0, (num_lats - 1), 0, (num_lons - 1), lataxis, lonaxis, num_contour_levels, contours, indexe, indexes, indexe1, indexes1);

  for (h=0; h<num_contour_levels; h++) {
    printf("processing contour level %d: %f\n", h, contours[h]);

    char hazard_string[128];
    sprintf(hazard_string, "WINDS_CASA");
        
    int k;
    int numpoints=0; //number of points in this particular contour
    int num_used_up_ips=0; //number of points total in all contours
    int numcontours=0; //number of discrete closed contours

    struct index_pair reference; //the i&j index reference of the cell that pointed to the current cell
    int ref_index; //the array that pointed the reference to the current_ip 0=indexe 1=indexe1 2=indexes 3=indexes1
    reference.xind = -1; 
    reference.yind = -1;
    
    struct index_pair *used_up_ips; //an array of all the index points we've used for all contours
    struct index_pair *contour_ips; //an array of the index points in the current contour
  
    used_up_ips = (struct index_pair*)malloc(sizeof(struct index_pair));
    contour_ips = (struct index_pair*)malloc(sizeof(struct index_pair));
    
    for (i=0; i<num_lats;i++) {
      for (j=0; j<num_lons;j=j+0) {
	
	//Print out the arrays... where each cell points and who points to it... sometimes multiple cells point to it/it points to multiple cells
	//printf("%d %d %d %d %d %d %d %d %d %d\n", i, j, indexe[h][i][j].xind, indexe[h][i][j].yind, indexes[h][i][j].xind, indexes[h][i][j].yind, indexe1[h][i][j].xind, indexe1[h][i][j].yind, indexes1[h][i][j].xind, indexes1[h][i][j].yind);
	
	//current_ip index_pair representing this i,j pair                                                                                                             
	struct index_pair current_ip;
	current_ip.xind = i;
	current_ip.yind = j;
	
	//each point in the grid has a three possibilities, post conrec.  
	
	//firstly, we may have already dealt with it... if so we move on
	if (contained(num_used_up_ips, current_ip, used_up_ips)) {
	  //printf("We've already dealt with this point\n");
	  ++j;
	  continue;
	}
	
	//secondly, it may not be part of a contour... if so, we add it to the list of points we've already dealt with and move on
	if ((indexe[h][i][j].xind == -1) && (indexes[h][i][j].xind == -1)) {
	  //printf("this point is not part of a contour\n");
	  used_up_ips[num_used_up_ips++] = current_ip;
	  num_used_up_ips++;
	  used_up_ips = (struct index_pair*)realloc(used_up_ips, (num_used_up_ips+1)*sizeof(struct index_pair));
	  ++j;
	  continue;
	}
	
	//thirdly, it may be both new AND part of a contour... 
	//printf("this is a new point and part of a contour\n");
	
	//view the arrays
	//printf("%d %d %d %d %d %d %d %d %d %d \n", i, j, indexe[h][i][j].xind, indexe[h][i][j].yind, indexes[h][i][j].xind, indexes[h][i][j].yind, indexe1[h][i][j].xind, indexe1[h][i][j].yind, indexes1[h][i][j].xind, indexes1[h][i][j].yind);
	
	//some points have multiple paths... these are harder to deal with... we want to eliminate them asap
	//lets check
	int multipath = 0;
	int kntr = 0; 
	if (indexe[h][i][j].xind != -1) 
	  ++kntr;
	if (indexes[h][i][j].xind != -1)
	  ++kntr;
	if (indexe1[h][i][j].xind != -1)
	  ++kntr;
	if (indexes1[h][i][j].xind != -1)
	  ++kntr;
	if (kntr > 2) {
	  multipath = 1;
	  //printf("this point is a multipath\n");
	}
	
	//now that that's out of the way, let's see if we have a point that referred us here...
	//if not this must be the first point in a new contour
	if(reference.xind == -1) {
	  //we don't want the first point to be a multipath... if it is, toss it back and move on
	  if (multipath == 1) {
	    ++j;
	    continue;
	  }
	  
	  json_coords_arr = json_array();
	  //otherwise start a new contour
	  //printf("this is the first point in a new contour\n");
	  //increment the number of points in this contour
	  //numpoints++;
	  //set the first point in this contour to this index point
	  
	  contour_ips[0] = current_ip;
	  
	  //now we start down the path
	  if(indexes[h][i][j].xind != -1) {
	    i = indexes[h][i][j].xind;
	    j = indexes[h][current_ip.xind][j].yind;
	    ref_index = 2;
	  }
	  else if (indexe[h][i][j].xind != -1){
	    i = indexe[h][i][j].xind;
	    j = indexe[h][current_ip.xind][j].yind;
	    ref_index = 0;
	  }
	  else if (indexes1[h][i][j].xind != -1){
	    i = indexes1[h][i][j].xind;
	    j = indexes1[h][current_ip.xind][j].yind;
	    ref_index = 3;
	  }
	  else if (indexe1[h][i][j].xind != -1){
	    i = indexe1[h][i][j].xind;
	    j = indexe1[h][current_ip.xind][j].yind;
	    ref_index = 1;
	  }
	  
	  //set the reference point and move on
	  reference = current_ip;
	  continue;
	}
	
	//ok so if we made it here we have a contour with at least one point already, and a reference
	//if this point is not already contained in the ip array of this contour, we are not closed yet
	//or it could be a multipath and we intentionally didn't add it to the array
	if (!contained(numpoints+1, current_ip, contour_ips)) {
	  //first let's figure out which index is the reference... 
	  int rfr;
	  if ((indexe[h][i][j].xind == reference.xind) && (indexe[h][i][j].yind == reference.yind))
	    rfr = 0;
	  else if ((indexe1[h][i][j].xind == reference.xind) && (indexe1[h][i][j].yind == reference.yind))
	    rfr = 1;
	  else if ((indexes[h][i][j].xind == reference.xind) && (indexes[h][i][j].yind == reference.yind))
	    rfr = 2;
	  else if ((indexes1[h][i][j].xind == reference.xind) && (indexes1[h][i][j].yind == reference.yind))
	    rfr = 3;
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
	    if ((indexe[h][i][j].xind != -1) && (rfr != 0)) {
	      i = indexe[h][i][j].xind;
	      j = indexe[h][current_ip.xind][j].yind;
	      ref_index = 0;
	    }
	    else if ((indexe1[h][i][j].xind != -1) && (rfr != 1)) {
	      i = indexe1[h][i][j].xind;
	      j = indexe1[h][current_ip.xind][j].yind;
	      ref_index = 1;
	    }
	    else if ((indexes[h][i][j].xind != -1) && (rfr != 2)) {
	      i = indexes[h][i][j].xind;
	      j = indexes[h][current_ip.xind][j].yind;
	      ref_index = 2;
	    }
	    else if ((indexes1[h][i][j].xind != -1) && (rfr != 3)) {
	      i = indexes1[h][i][j].xind;
	      j = indexes1[h][current_ip.xind][j].yind;
	      ref_index = 3;
	    }
	    else 
	      printf("uh oh\n");
	    
	    reference = current_ip;
	    continue;
	  }
      
	  //if we're here it's a multipath point in a contour
	  //printf("this is a multipath point in a contour\n");
	  //we have to toss this point and replace references to this point
	  //with references to another point... probably the closest point by proximity 
	  double dx2p;
	  double mndx = 999999999999.0; //Huge val
	  struct index_pair destination;
	  int which_index;  //0=indexe 1=indexe1 2=indexes 3=indexes1
	  
	  if ((indexe[h][i][j].xind != -1) && (rfr != 0)) {
	    inverse_vincenty(lat_in[reference.xind], lon_in[reference.yind], lat_in[indexe[h][i][j].xind], lon_in[indexe[h][i][j].yind], &dx2p);
	    if (dx2p < mndx) {
	      mndx = dx2p;
	      destination.xind = indexe[h][i][j].xind;
	      destination.yind = indexe[h][i][j].yind;
	      which_index = 0;
	    }
	  }
	  
	  if ((indexe1[h][i][j].xind != -1) && (rfr != 1)) {
	    inverse_vincenty(lat_in[reference.xind], lon_in[reference.yind], lat_in[indexe1[h][i][j].xind], lon_in[indexe1[h][i][j].yind], &dx2p);
	    if (dx2p < mndx) {
	      mndx = dx2p;
	      destination.xind = indexe1[h][i][j].xind;
	      destination.yind = indexe1[h][i][j].yind;
	      which_index = 1;
	    }
	  }
	  
	  if ((indexes[h][i][j].xind != -1) && (rfr != 2)) {
	    inverse_vincenty(lat_in[reference.xind], lon_in[reference.yind], lat_in[indexes[h][i][j].xind], lon_in[indexes[h][i][j].yind], &dx2p);
	    if (dx2p < mndx) {
	      mndx = dx2p;
	      destination.xind = indexes[h][i][j].xind;
	      destination.yind = indexes[h][i][j].yind;
	      which_index = 2;
	    }
	  }
	  
	  if ((indexes1[h][i][j].xind != -1) && (rfr != 3)) {
	    inverse_vincenty(lat_in[reference.xind], lon_in[reference.yind], lat_in[indexes1[h][i][j].xind], lon_in[indexes1[h][i][j].yind], &dx2p);
	    if (dx2p < mndx) {
	      mndx = dx2p;
	      destination.xind = indexes1[h][i][j].xind;
	      destination.yind = indexes1[h][i][j].yind;
	      which_index = 3;
	    }
	  }
	  
	  //set the reference to the destination
	  switch (ref_index) {
	  case 0:
	    indexe[h][reference.xind][reference.yind] = destination;
	    break;
	  case 1:
	    indexe1[h][reference.xind][reference.yind] = destination;
	    break;
	  case 2:
	    indexes[h][reference.xind][reference.yind] = destination;
	    break;
	  case 3:
	    indexes1[h][reference.xind][reference.yind] = destination;
	    break;
	  }
	  
	  //change the destination's reference to the current_ip over to the current ip's reference
	  if ((indexe[h][destination.xind][destination.yind].xind == current_ip.xind) && (indexe[h][destination.xind][destination.yind].yind == current_ip.yind))
	    indexe[h][destination.xind][destination.yind] = reference;
	  else if ((indexe1[h][destination.xind][destination.yind].xind == current_ip.xind) && (indexe1[h][destination.xind][destination.yind].yind == current_ip.yind))
	    indexe1[h][destination.xind][destination.yind] = reference;
	  else if((indexes[h][destination.xind][destination.yind].xind == current_ip.xind) && (indexes[h][destination.xind][destination.yind].yind == current_ip.yind))
	    indexes[h][destination.xind][destination.yind] = reference;
	  else if((indexes1[h][destination.xind][destination.yind].xind == current_ip.xind) && (indexes1[h][destination.xind][destination.yind].yind == current_ip.yind))
	    indexes1[h][destination.xind][destination.yind] = reference;
	  else 
	    printf("SOS!\n");
	  
	  //scrub the current_ip's path to the reference and the destination
	  switch (rfr) {
	  case 0:
	    indexe[h][i][j].xind = -1;
	    indexe[h][i][j].yind = -1;
	    break;
	  case 1:
	    indexe1[h][i][j].xind = -1;
	    indexe1[h][i][j].yind = -1;
	    break;
	  case 2:
	    indexes[h][i][j].xind = -1;
	    indexes[h][i][j].yind = -1;
	    break;
	  case 3:
	    indexes1[h][i][j].xind = -1;
	    indexes1[h][i][j].yind = -1;
	    break;
	  }
	  switch (which_index) {
	  case 0:
	    indexe[h][i][j].xind = -1;
	    indexe[h][i][j].yind = -1;
	    break;
	  case 1:
	    indexe1[h][i][j].xind = -1;
	    indexe1[h][i][j].yind = -1;
	    break;
	  case 2:
	    indexes[h][i][j].xind = -1;
	    indexes[h][i][j].yind = -1;
	    break;
	  case 3:
	    indexes1[h][i][j].xind = -1;
	    indexes1[h][i][j].yind = -1;
	    break;
	  }
	  
	  //we may need to further modify the path's of the other points the current multipath ip points to
	  
	  //move to the destination
	  i = destination.xind;
	  j = destination.yind;
	  continue;
	}
	
	if (contained(numpoints+1, current_ip, contour_ips)) {
	  //if the point is contained we are at the end.... maybe?
	  //increment the number of contours
	  ++numcontours;
	  printf("numcontours: %d\n", numcontours);
	  for (k=0; k<numpoints+1; ++k) {
	    //print out the contour
	    printf("contour %d pair %d i %d j %d\n", numcontours, k, contour_ips[k].xind, contour_ips[k].yind);
	    used_up_ips[num_used_up_ips] = contour_ips[k];
	    num_used_up_ips++;
	    used_up_ips = (struct index_pair*)realloc(used_up_ips, (num_used_up_ips+1)*sizeof(struct index_pair));
	  }
	  
	  //dump contour to output 
	  //double dxbetween;
	  int a;
	  if (numpoints > min_contour_points) {
	    num_valid_contours++;
	    for (a=0; a < numpoints+1; a++) {
	      //  if (a == numpoints-1) {
	      //  inverse_vincenty(lat_in[contour_ips[a].xind], lon_in[contour_ips[a].yind], lat_in[contour_ips[0].xind], lon_in[contour_ips[0].yind], &dxbetween);
	      //if (dxbetween < 1000) {
	      if (json_output == 1) {
		json_t * json_coords = json_pack("[f,f]", lon_in[contour_ips[a].yind], lat_in[contour_ips[a].xind]);
		json_array_append_new(json_coords_arr, json_coords);
	      }
	      if (xml_output == 1) {
		ll_array[edgecount].lat = lat_in[contour_ips[a].xind];
		ll_array[edgecount].lon = lon_in[contour_ips[a].yind];
		++edgecount;
		ll_array = (struct latLon*)realloc(ll_array, (edgecount+1)*sizeof(struct latLon));
	      }
	    }
	  
	    /* now append the first point back to the end to conform with geoJSON polygon standard */
	    if (json_output == 1) {
	      json_t * json_coords = json_pack("[f,f]", lon_in[contour_ips[0].yind], lat_in[contour_ips[0].xind]);
	      json_array_append_new(json_coords_arr, json_coords);
	    }
	    if (xml_output == 1) {
	      ll_array[edgecount].lat = lat_in[contour_ips[0].xind];
	      ll_array[edgecount].lon = lon_in[contour_ips[0].yind];
	      ++edgecount;
	      ll_array = (struct latLon*)realloc(ll_array, (edgecount+1)*sizeof(struct latLon));
	    }
	  
	    float windMPH = contours[h] * 2.23694;
            char magstring[8];
            sprintf(magstring, "%.0fMG", windMPH);

	    char contour_level_string[16];
	    sprintf(contour_level_string, "%.0fMPH", windMPH);

	    char id_string[128];
	    
	    sprintf(id_string, "WINDS_CASA_%s_%s_%d", starttime, contour_level_string, numcontours);

	    if (json_output == 1) {
	    
	      json_t *polygon = json_pack("{s:s,s:s,s:s,s:{s:s,s:s},s:{s:s,s:[o]}}", "type", "Feature", "id", id_string, "hazardType", hazard_string, "properties", "timestamp", alt_starttime, "mag", magstring, "geometry", "type", "Polygon", "coordinates", json_coords_arr);
	      json_array_append_new(features_json_arr, polygon);
	    }
	    if (xml_output == 1) {
	      char val[32];
	      mvt_edges = xmlNewChild(root,NULL,BAD_CAST "contours", NULL);
	      contour_number = xmlNewChild(mvt_edges,NULL, BAD_CAST "contour", NULL);
	      sprintf(val, "%d", numcontours);
	      xmlNewChild(contour_number,NULL,BAD_CAST "contour_number", BAD_CAST val);
	      contour_points = xmlNewChild(contour_number,NULL, BAD_CAST "contour_points", NULL);
	      xmlNodePtr point;
	      for (y=0; y < numpoints; ++y) {
		point = xmlNewChild(contour_points, NULL, BAD_CAST "point", NULL);
		sprintf(val, "%d", y);
		xmlNewChild(point, NULL, BAD_CAST "point_number", BAD_CAST val);
		sprintf(val, "%10.4f", ll_array[y].lat);
		xmlNewChild(point,NULL, BAD_CAST "latitude", BAD_CAST val);
		sprintf(val, "%10.4f", ll_array[y].lon);
		xmlNewChild(point,NULL, BAD_CAST "longitude", BAD_CAST val);
	      }
	    }
	  }
	  if ((mcc_output == 1) && (numpoints > 8)) {
	    for (y=0; y < numpoints; ++y) {
	      char detect[128];
	      sprintf(detect, "winds :lat %.4f :long %.4f :height 1000.0 :time %s\n", ll_array[y].lat, ll_array[y].lon, mcctime);
	      if (connected) {
		write(mcc_sock,detect,strlen(detect));
	      }
	    }
	  }
	  
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

  if (json_output == 1) {
    size_t flags = 0;
    flags |= JSON_INDENT(1);
    flags |= JSON_PRESERVE_ORDER;
    flags |= JSON_REAL_PRECISION(6);
    //json_error_t error;


    //ACS test
    printf("num_valid_contours: %d\n", num_valid_contours);
    if (acs_output == 1){

      jsonbuf = json_dumps( jsonroot, flags );
      char outstr[102800];
      sprintf(outstr, "&&&%d&&&%s", (int)strlen(jsonbuf), jsonbuf);
      size_t strn = strnlen(outstr, 102800);
      printf("strn: %d outstr: %s\n", (int)strn, outstr);
      char tmp[strn];

      int b;

      if (num_valid_contours > 0) {
	for (b=0; b<num_acs_hosts; b++) {
	  if (sock_connected[b]) {
	    send(acs_sock[b], outstr, sizeof(tmp), 0);
	  }
	}
      }
    }
    json_dump_file(jsonroot, jsonname, flags);
    json_decref( jsonroot );
    json_decref( features_json_arr );
  }

  if (xml_output == 1) {
    xmlDocDumpFormatMemory(doc,(xmlChar**)&xmlbuf,&xmlsize,1);
    xmlFreeDoc(doc);
    xmlCleanupParser();
    if(strlen(xmlname)!=0){
      xmlfile=fopen(xmlname,"w");
      if(xmlfile){
	fprintf(xmlfile,"%s\n",xmlbuf);
	fclose(xmlfile);
      } else {
	printf("could not write xml file\n");
      }
    }
  }

  if (mcc_output == 1) {
    if (connected) 
      close(mcc_sock);
  }

  if (acs_output == 1) {
    int b;
    for (b=0; b<num_acs_hosts; b++) {
      if (sock_connected[b]) {
	close(acs_sock[b]);
      }
    }
  }
  
  free(data_in);
  free(indexe);
  free(indexes);
  free(indexe1);
  free(indexes1);
  free(lat_in);
  free(lon_in);
  free(lonaxis);
  free(lataxis);
  free(ll_array);
  
  return(0);
}
