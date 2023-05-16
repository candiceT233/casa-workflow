#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <argp.h>
#include <jansson.h>
#include <libconfig.h>

#define BADVAL         -99900

static struct argp_option options[] = {
  {"configfile", 'c', "configuration_file",    0,   "Specify name of configuration file"},
  {"output_filename", 'o', "output_geojson_filename",   0,   "Specify name of output geoJSON locations file"},
  {"geojson_locations",  'g', "input_geojson_locations_file",      0,   "Specify name of input geoJSON locations file"},
  {"property", 'p', 0, 0,   "Append the QPE value to the geojson file in the properties tag"},
  {"email", 'e', 0, 0,   "Send an email alert to users identified in the config file"},
  {"AlertSent_set",'s', "no|yes",              0,   "Set AlertSent variable to yes or no"},
  { 0 }
};

#define MAX_NAME 1024

struct arguments{
  char filename[MAX_NAME];
  char config_filename[MAX_NAME];
  char output_filename[MAX_NAME];
  char geojson_filename[MAX_NAME];
  char alert_set[MAX_NAME];
  short alert_option;
  short property;
  short email;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state){
  struct arguments *arguments = state->input;

  switch(key){
  case 'c':
    if (arg[0] != '-')
      strncpy(arguments->config_filename,arg,MAX_NAME);
    else {
      printf("arguments should not begin with a - sign.  Exiting...\n");
      exit(0);
    }
    break;
  case 'g':
    if (arg[0] != '-')
      strncpy(arguments->geojson_filename,arg,MAX_NAME);
    else {
      printf("arguments should not begin with a - sign.  Exiting...\n");
      exit(0);
    }
    break;
  case 'o':
    if (arg[0] != '-')
      strncpy(arguments->output_filename,arg,MAX_NAME);
    else {
      printf("arguments should not begin with a - sign.  Exiting...\n");
      exit(0);
    }
    break;
  case 'p':
    arguments->property = 1;
    break;
  case 'e':
    arguments->email = 1;
    break;
  case 's':
    if (arg[0] != '-')
      strncpy(arguments->alert_set,arg,MAX_NAME);
    else {
      printf("arguments should not begin with a - sign.  Exiting...\n");
      exit(0);
    }
    if ((strcmp(arg, "no") != 0) && (strcmp(arg, "yes") != 0)) {
      printf("arguments to -s should be yes or no. Exiting...\n");
      exit(0);
    }
    arguments->alert_option = 1;
    break;
  case ARGP_KEY_ARG:
    if(state->arg_num>=1)
      argp_usage(state);
    strncpy(arguments->filename,arg,MAX_NAME);
    break;
  case ARGP_KEY_END:
    if(state->arg_num<1)
      argp_usage(state);
    break;
  default:
    return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

static struct argp argp = {options,parse_opt,"<input_geojson_data_file>","$Id: pointAlert.c,v 1.0 2018-10-24 11:34:23 elyons Exp $"};

int email_someone(const char *to_line, const char *from_line, const char *subject_line, const char *message_line) {
  FILE *email_connection = popen("/usr/lib/sendmail -t", "w");
  if (email_connection != NULL) {
    fprintf(email_connection, "To: %s\n", to_line);
    fprintf(email_connection, "From: %s\n", from_line);
    fprintf(email_connection, "Subject: %s\n\n", subject_line);
    fwrite(message_line, 1, strlen(message_line), email_connection);
    fwrite(".\n", 1, 2, email_connection);
    pclose(email_connection);
    return 0;
  }
  else {
    perror("Error sending mail!");
    return -1;
  }
}

int pnpoly(int nvert, float *vertx, float *verty, float testx, float testy) {
  /*
  Copyright (c) 1970-2003, Wm. Randolph Franklin                          
  Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: 
  Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.
  Redistributions in binary form must reproduce the above copyright notice in the documentation and/or other materials provided with the distribution.
  The name of W. Randolph Franklin may not be used to endorse or promote products derived from this Software without specific prior written permission.
  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                                                      

  */
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
      c = !c;
  }
  return c;
}

int main(int argc, char* argv[])
{
  struct arguments arguments;
  arguments.config_filename[0] = '\0';
  arguments.output_filename[0] = '\0';
  arguments.alert_option = 0;
  arguments.property = 0;
  arguments.email = 0;
  arguments.geojson_filename[0] = '\0';
  int geojson_specified = 0;
  json_t *empty = json_object();
  
  argp_parse(&argp,argc,argv,0,0,&arguments);

  if (arguments.geojson_filename[0] == '\0') {
    printf("No geojson location input specified. There is nothing to do!\n");
    exit(0);
  }
  else {
    geojson_specified = 1;
  }

  if (arguments.output_filename[0] == '\0') {
    printf("No output filename specified.  Writing back to the same geojson input filename\n");
    sprintf(arguments.output_filename, "%s", arguments.geojson_filename);
  }
  else {
    int dumpstatus = json_dump_file(empty, arguments.output_filename, 0);
    if (dumpstatus == -1) {
      printf("GeoJSON file dump failed!\n");
    }
  }
  
  char *pointAlert_home;
  if (arguments.config_filename[0] == '\0') {
    if ((pointAlert_home = getenv("POINT_ALERT_HOME")) != NULL) {
      sprintf(arguments.config_filename, "%s/pointAlert_config.txt", pointAlert_home);
    }
    else {
      printf("Please specify config file with -c option or declare POINT_ALERT_HOME environment variable.  exiting...\n");
      exit(0);
    }
  }
  //printf("config_filename: %s\n", arguments.config_filename);

  config_t pointAlert_config;
  
  config_init(&pointAlert_config);

  if (!config_read_file(&pointAlert_config, arguments.config_filename)) {
    printf("config file not found or contains errors.  exiting...\n");
    config_destroy(&pointAlert_config);
    exit(0);
  }

  const char *email_from;
  const config_setting_t *email_to_arr;
  int ema;
  int num_email_recipients = 1;
  //double email_threshold;
  config_setting_t *email_to_params;
  const char *data_type;

  if(!config_lookup_string(&pointAlert_config, "data_type", &data_type)) {
    printf("No data_type definition. exiting...\n");
    exit(0);
  }
  else {
    //printf("Data type: %s\n", data_type);
  }
  
  if (arguments.email) {
    
    if(!config_lookup_string(&pointAlert_config, "email_from", &email_from)) {
      printf("No email_from definition.  exiting...\n");
      exit(0);
    }
    else {
      //printf("Email from: %s\n", email_from);
    }
        
    email_to_params = config_lookup(&pointAlert_config, "email_to_params");
    
    if (email_to_params != NULL) {
      email_to_arr = config_setting_lookup(email_to_params, "email_to");
      if (email_to_arr != NULL) {
	num_email_recipients = config_setting_length(email_to_arr);
      }
      else {
	printf("No Email Recipients defined.  Exiting.\n");
	exit(0);
      }
    }
  }
  
  const char *email_addresses[num_email_recipients];

  if (arguments.email) {
    for (ema = 0; ema < num_email_recipients; ema++) {
      email_addresses[ema] = config_setting_get_string_elem(email_to_arr, ema);
      if (!email_addresses[ema]) {
	printf("Email Address %d not properly defined.  Exiting...\n", ema);
	exit(0);
      }
      else {
	//printf("Emailing: %s\n", email_addresses[ema]);
      }
    }
  }
  
  //geojson location input
  if (geojson_specified) {
    json_t *json, *pjson;
    json_error_t error;
    json_t *features, *polyfeatures;
    size_t num_features, num_polyfeatures;
    int k,l;
    json_t *yes = json_string("yes");
    json_t *no = json_string("no");

    //first load location input data
    json = json_load_file(arguments.geojson_filename, 0, &error);
    if(!json) {
      printf("JSON ERROR: %s\n", error.text);
      exit(-1);
    }

    if (!json_is_object(json)){
      printf("json type: %d\n", json_typeof(json));
      printf("geoJSON point file does not contain an object.  Please try another geoJSON file. \n");
      exit(0);
    }

    features = json_object_get(json, "features");
    if (features == NULL) {
      printf("geoJSON point file does not contain the necessary \"features\" field.  Exiting...\n");
      exit(0);
    }
    if (!json_is_array(features)) {
      printf("\"features\" field in point file is not an array.  Exiting...");
      exit(0);
    }

    num_features = json_array_size(features);
    if (num_features < 1) {
      printf("no features present in this point geoJSON.  Exiting...");
      exit(0);
    }
    
    //point file is ok, now check polygon file

    pjson = json_load_file(arguments.filename, 0, &error);
    if(!pjson) {
      printf("JSON ERROR: %s\n", error.text);
      exit(-1);
    }
    
    if (!json_is_object(pjson)){
      printf("json type: %d\n", json_typeof(pjson));
      printf("geoJSON polygon file does not contain an object.  Please try another geoJSON file. \n");
      exit(0);
    }
    
    polyfeatures = json_object_get(pjson, "features");
    if (features == NULL) {
      printf("geoJSON polygon file does not contain the necessary \"features\" field.  Exiting...\n");
      exit(0);
    }
    if (!json_is_array(polyfeatures)) {
      printf("\"features\" field in polygon file is not an array.  Exiting...");
      exit(0);
    }
    
    num_polyfeatures = json_array_size(polyfeatures);
    if (num_polyfeatures < 1) {
      printf("no features present in this polygon geoJSON.  Exiting...");
      exit(0);
    }

    for (k=0; k<num_features; k++) {
      json_t *this_feature = json_array_get(features, k);
      //float lat = BADVAL;
      //float lon = BADVAL;
      if (json_is_object(this_feature)) {
	json_t *this_feature_type = json_object_get(this_feature, "type");
	if (json_is_string(this_feature_type)) { 
	  if (strcmp(json_string_value(this_feature_type), "Feature") == 0) {
	    //const char *key;
	    //json_t *value;
	    //json_object_foreach(this_feature, key, value) { printf("key %s val %d\n", key, json_typeof(value));}
	    json_t *this_feature_properties = json_object_get(this_feature, "properties");
	    
	    if(arguments.alert_option == 1) {
	      if (strcmp(arguments.alert_set, "yes") == 0) {
		int setAlert = json_object_set(this_feature_properties, "AlertSent", yes);
		if (setAlert == -1) {
		  printf("Failed to set AlertSent property\n");
		}
	      }
	      else if (strcmp(arguments.alert_set, "no") == 0) {
		int setAlert = json_object_set(this_feature_properties, "AlertSent", no);
		if (setAlert == -1) {
		  printf("Failed to set AlertSent property\n");
		}
	      }
	      else {
		printf("argument not recognized. exiting...");
		exit(0);
	      }
	    }
	    else {

	      int zero = 0;
	      int resetAlert = json_object_set(this_feature_properties, "alert", json_integer(zero));
	      json_object_del(this_feature_properties, "hazardType");
	      json_object_del(this_feature_properties, "timestamp");
	      json_object_del(this_feature_properties, data_type);
	      if (resetAlert == -1) {
		printf("failed to set the alert property\n");
	      }
	      
	      json_t *this_feature_alert = json_object_get(this_feature_properties, "AlertSent");
	      json_t *this_feature_id = json_object_get(this_feature_properties, "name");
	      //json_object_foreach(this_feature_properties, key, value) { printf("key %s val %s\n", key, json_string_value(value));}

	      json_t *this_feature_geometry = json_object_get(this_feature, "geometry");
	      if (json_is_object(this_feature_geometry)) {
		json_t *this_feature_geometry_type = json_object_get(this_feature_geometry, "type");
		if (json_is_string(this_feature_geometry_type)) {
		  if (strcmp(json_string_value(this_feature_geometry_type), "Point") == 0) {
		    json_t *this_coord_pair = json_object_get(this_feature_geometry, "coordinates");
		    if (json_is_array(this_coord_pair)) {
		      if (json_array_size(this_coord_pair) == 2) {
			json_t *this_lon = json_array_get(this_coord_pair, 0);
			json_t *this_lat = json_array_get(this_coord_pair, 1);
			if ((json_is_real(this_lon)) && (json_is_real(this_lat))) {
			  float lat = json_real_value(this_lat);
			  float lon = json_real_value(this_lon);
			  //printf("lat: %f lon: %f\n", json_real_value(this_lat), json_real_value(this_lon));
			  
			  //ok now parse the data file and look for objects in there
			  for (l=0; l<num_polyfeatures; l++) {
			    json_t *this_polyfeature = json_array_get(polyfeatures, l);
			    if (json_is_object(this_polyfeature)) {
			      json_t *this_polyfeature_type = json_object_get(this_polyfeature, "type");
			      //json_t *this_polyfeature_id = json_object_get(this_polyfeature, "id");
			      if (json_is_string(this_polyfeature_type)) {
				if (strcmp(json_string_value(this_polyfeature_type), "Feature") == 0) {
				  json_t *this_polyfeature_hazard = json_object_get(this_polyfeature, "hazardType");
				  json_t *this_polyfeature_properties = json_object_get(this_polyfeature, "properties");
				  json_t *this_polyfeature_time = json_object_get(this_polyfeature_properties, "timestamp");
				  json_t *this_polyfeature_val = json_object_get(this_polyfeature_properties, data_type);
				  json_t *this_polyfeature_geometry = json_object_get(this_polyfeature, "geometry");
				  if (json_is_object(this_polyfeature_geometry)) {
				    json_t *this_polyfeature_geometry_type = json_object_get(this_polyfeature_geometry, "type");
				    if (json_is_string(this_polyfeature_geometry_type)) {
				      if (strcmp(json_string_value(this_polyfeature_geometry_type), "Polygon") == 0) {
					json_t *this_polyfeature_geometry_coordinates = json_object_get(this_polyfeature_geometry, "coordinates");
					if (json_is_array(this_polyfeature_geometry_coordinates)) {
					  size_t num_polygons = json_array_size(this_polyfeature_geometry_coordinates);
					  int p;

					  for (p=0; p<num_polygons; p++) {
					    json_t *this_polygon = json_array_get(this_polyfeature_geometry_coordinates, p);
					    if (json_is_array(this_polygon)) {
					      size_t num_coords = json_array_size(this_polygon);
					      int c;
					      float *latpoints, *lonpoints;
					      size_t goodvals = 0;
					      latpoints = malloc(num_coords * sizeof(float));
					      lonpoints = malloc(num_coords * sizeof(float));
					      for (c=0; c<num_coords; c++) {
						json_t *llpair = json_array_get(this_polygon, c);
						if (json_is_array(llpair)) {
						  if (json_array_size(llpair) == 2) {
						    if ((json_is_real(json_array_get(llpair, 0))) && (json_is_real(json_array_get(llpair, 1)))) {
						      lonpoints[c] = json_real_value(json_array_get(llpair, 0));
						      latpoints[c] = json_real_value(json_array_get(llpair, 1));
						      goodvals++;
						    }
						  }
						}
					      }
					      if (goodvals == num_coords) {
						int inside = pnpoly(num_coords, lonpoints, latpoints, lon, lat);
						if (inside == 1) {
						  if(arguments.property) {
						    int alert = 1;
						    int setAlert = json_object_set(this_feature_properties, "alert", json_integer(alert));
						    int settimestatus = json_object_set(this_feature_properties, "timestamp", this_polyfeature_time);
						    int setHazardType = json_object_set(this_feature_properties, "hazardType", this_polyfeature_hazard);
						    int setDataVal = json_object_set(this_feature_properties, data_type, this_polyfeature_val);
						    if ((setAlert == -1) || (settimestatus == -1) || (setDataVal == -1) || (setHazardType == -1)) {
						      printf("GeoJSON property append failed!\n");
						    }
						  }
						  if(arguments.email) {
						    if(strcmp(json_string_value(this_feature_alert), "no") == 0) {
						      char email_subject[256];
						      sprintf(email_subject, "CASA automated notification for %s\n", json_string_value(this_feature_id));
						      char email_message[1028];
						      sprintf(email_message, "siteLocation: %s\nalert type: %s %s: %s\ntimestamp: %s\n\n\nTo be removed from this mailing list, please send email to elyons@engin.umass.edu", json_string_value(this_feature_id), json_string_value(this_polyfeature_hazard), data_type, json_string_value(this_polyfeature_val), json_string_value(this_polyfeature_time));
						      for (ema = 0; ema < num_email_recipients; ema++) {
							email_someone(email_addresses[ema], email_from, email_subject, email_message);
						      }
						      int setAlert = json_object_set(this_feature_properties, "AlertSent", yes);
						      if (setAlert == -1) {
							printf("Failed to set AlertSent property\n");
						      }
						    }
						  }
						}
					      }
					      free(latpoints);
					      free(lonpoints);
					    }
					  }
					}
				      }
				    }
				  }
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    int dumpstatus = json_dump_file(json, arguments.output_filename, 0);
    if (dumpstatus == -1) {
      printf("GeoJSON file dump failed!\n");
    }
  }

  return 0;
}
