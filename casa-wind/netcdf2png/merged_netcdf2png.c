/* $Id: netcdf2png.c,v 1.2 2007-11-27 01:47:09 luko Exp $
   modified 3/13/08 elyons SNR Filtering Exp $
   Copyright 2005-2006 Dynamic Sensing Technologies (All rights reserved) */
#define _GNU_SOURCE
#include <png.h>
#include <netcdf.h>
#include <math.h>
#include <stdlib.h>
#include <argp.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <time.h>

#define ERROR -1
#define TRUE 1
#define FALSE 0

struct picdimensions{
  int width;
  int height;
  double llx,lly,llz;
  double ulx,uly,ulz;
  double lrx,lry,lrz;
  double plotmin;
  double plotmax;
  int transparency;
  int palletsize;
  int *pallet;
};

struct mergedimensions{
  size_t num_lats;
  size_t num_lons;
  char timeformat[16];
};


/* Write XML file describing what we're doing - netcdf filename, etc... */
int write_xml(char *filename,struct picdimensions *picdims,
	      struct mergedimensions *datadims,int ncid,
	      char **xmlbuf, int *buffersize){
  int i;
  char attname[1024];
  char value[1024];
  xmlDocPtr doc;
  xmlNodePtr root,ncatts,ddims,pdims;
    
  doc=xmlNewDoc(BAD_CAST "1.0");
  root=xmlNewNode(NULL,BAD_CAST "root");
  xmlDocSetRootElement(doc,root);

  /* Add filename */
  xmlNewChild(root,NULL,BAD_CAST "filename",BAD_CAST filename);
  
  /* Add timestamp */
  sprintf(value,"%s",datadims->timeformat);
  xmlNewChild(root,NULL,BAD_CAST "time", BAD_CAST value);
  
  /* Add netcdf global attributes */
  ncatts = xmlNewChild(root,NULL,BAD_CAST "netcdf_attributes",NULL);
  for(i=0;;i++){
    if(nc_inq_attname(ncid,NC_GLOBAL,i,attname)==NC_NOERR){
      nc_type atttype;
      size_t length;
      /* Get attribute type */
      nc_inq_att(ncid,NC_GLOBAL,attname,&atttype,&length);
      switch(atttype){
	int ival; float fval; double dval;
      case NC_CHAR:
	nc_get_att_text(ncid,NC_GLOBAL,attname,value);
	value[length]='\0';
	break;
      case NC_INT:
	nc_get_att_int(ncid,NC_GLOBAL,attname,&ival);	
	sprintf(value,"%d",ival);
	break;
      case NC_FLOAT:
	nc_get_att_float(ncid,NC_GLOBAL,attname,&fval);	
	sprintf(value,"%12.8f",fval);
	break;
      case NC_DOUBLE:
	nc_get_att_double(ncid,NC_GLOBAL,attname,&dval);
	sprintf(value,"%12.8f",dval);
	break;
      default:
	printf("File contains unsupported attribute\n");
	continue;
      }
      xmlNewChild(ncatts,NULL,BAD_CAST attname,BAD_CAST value);
    } else {
      break;
    }
  }
  
  /* Add mergedimensions */
  ddims = xmlNewChild(root,NULL,BAD_CAST "mergedimensions",NULL);
  sprintf(value,"%d",(int)(datadims->num_lats));
  xmlNewChild(ddims,NULL,BAD_CAST "num_lats",BAD_CAST value);
  sprintf(value,"%d",(int)(datadims->num_lons));
  xmlNewChild(ddims,NULL,BAD_CAST "num_lons",BAD_CAST value);
 
  /* Add picdimensions */
  pdims = xmlNewChild(root,NULL,BAD_CAST "picdimensions",NULL);
  sprintf(value,"%d",picdims->width);
  xmlNewChild(pdims,NULL,BAD_CAST "width",BAD_CAST value);
  sprintf(value,"%d",picdims->height);
  xmlNewChild(pdims,NULL,BAD_CAST "height",BAD_CAST value);
  sprintf(value,"%f,%f,%f",picdims->llx,picdims->lly,picdims->llz);
  xmlNewChild(pdims,NULL,BAD_CAST "lowerleft",BAD_CAST value);
  sprintf(value,"%f,%f,%f",picdims->ulx,picdims->uly,picdims->ulz);
  xmlNewChild(pdims,NULL,BAD_CAST "upperleft",BAD_CAST value);
  sprintf(value,"%f,%f,%f",picdims->lrx,picdims->lry,picdims->lrz);
  xmlNewChild(pdims,NULL,BAD_CAST "lowerright",BAD_CAST value);
  sprintf(value,"%f",picdims->plotmin);
  xmlNewChild(pdims,NULL,BAD_CAST "plotmin",BAD_CAST value);
  sprintf(value,"%f",picdims->plotmax);
  xmlNewChild(pdims,NULL,BAD_CAST "plotmax",BAD_CAST value);
  sprintf(value,"%d",picdims->palletsize);
  xmlNewChild(pdims,NULL,BAD_CAST "palletsize",BAD_CAST value);
 
  /* Put XML file in specified buffer */
  xmlDocDumpFormatMemory(doc,(xmlChar**)xmlbuf,buffersize,1);
  xmlFreeDoc(doc);
  xmlCleanupParser();

  return(0);
}


int make_pallet_from_png(char *filename, int *palletsize, int **pallet){
  FILE *fp;
  int png_transforms=0;
  png_uint_32 height, width;
  int bit_depth,color_type,interlace_type,compression_type, filter_method;
  int i;
  png_byte **row_pointers;

  fp=fopen(filename,"rb");
  if (!fp){
    return (ERROR);
  }

  /* Open PNG file for reading */
  png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
					       NULL,NULL,NULL);
  if (!png_ptr)
    return (ERROR);
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr){
    png_destroy_read_struct(&png_ptr,(png_infopp)NULL, (png_infopp)NULL);
    return (ERROR);
  }
  png_infop end_info = png_create_info_struct(png_ptr);
  if (!end_info){
    png_destroy_read_struct(&png_ptr, &info_ptr,(png_infopp)NULL);
    return (ERROR);
  }
  if (setjmp(png_jmpbuf(png_ptr))){
    png_destroy_read_struct(&png_ptr, &info_ptr,&end_info);
    fclose(fp);
    return (ERROR);
  }
  png_init_io(png_ptr, fp);

  png_read_png(png_ptr, info_ptr, png_transforms, NULL);

  /* Get color scale data */
  row_pointers = png_get_rows(png_ptr, info_ptr);

  png_get_IHDR(png_ptr, info_ptr,&width,&height,
	       &bit_depth, &color_type, &interlace_type,
	       &compression_type, &filter_method);

  /* Check to see that this is an 8-bit, RGB[a] Image */
  if(bit_depth!=8){
    printf("Color depth not 8 bits!\n");
    png_destroy_read_struct(&png_ptr,&info_ptr,&end_info);
    fclose(fp);
    return (ERROR);
  }

  /* Make the pallet */
  *palletsize=width;
  (*pallet)=(int*)malloc((width)*sizeof(int));

  switch(color_type){
  case PNG_COLOR_TYPE_RGB:
    for(i=0;i<width;i++){
      (*pallet)[i] =((int)(row_pointers[0])[i*3+2]);
      (*pallet)[i]+=((int)(row_pointers[0])[i*3+1])<<8;
      (*pallet)[i]+=((int)(row_pointers[0])[i*3+0])<<16;
      (*pallet)[i]+=0xFF000000;
    }
    break;
  case PNG_COLOR_TYPE_RGBA:
    for(i=0;i<width;i++){
      (*pallet)[i] =((int)((row_pointers[0])[i*4+0]));
      (*pallet)[i]+=((int)((row_pointers[0])[i*4+1]))<<8;
      (*pallet)[i]+=((int)((row_pointers[0])[i*4+2]))<<16;
      (*pallet)[i]+=0xFF000000;
    }
    break;
  default:
    printf("color type unsupported!!\n");
    break;
  }
  
  png_destroy_read_struct(&png_ptr,&info_ptr,&end_info);
  fclose(fp);

  return 0;
}


int number2color(int palletsize,int *pallet,int transparency,
		 double min,double max,double val){
  union {
    char c[4];
    int  i;
  } pix_num;
  double temp;
  int itemp;

  temp=(val-min)/(max-min);
  temp*=(palletsize-1);
  itemp=palletsize-(int)rint(temp);
  if(itemp<0){
    pix_num.i=pallet[0];
    pix_num.c[3]=transparency;
  }else if(itemp>(palletsize-1)){
    /* Make off-range in the minus direction fully transparent */
    pix_num.i=0;
  }else{
    pix_num.i=pallet[itemp];
    pix_num.c[3]=transparency;
  }


  return(pix_num.i);
}

void data2picture(float **databuf, struct mergedimensions *datadims,
		  char **picbuf, struct picdimensions *picdims){
  int i,j;

  double pix_val;
  union {
    char c[4];
    int  i;
  } pix_num;
      
  /* Perform resampling to plotting coordinates */
  for(j=0;j<373;j++){
    for(i=0;i<373;i++){
      if (databuf[j][i] == 0)
	continue;
      pix_val=databuf[j][i];
      pix_num.i=number2color(picdims->palletsize,picdims->pallet,
			     picdims->transparency,
			     picdims->plotmin,picdims->plotmax,pix_val);
      /* Place pixel in image buffer */
      int dj = 0; int di = 0;
      int endj = 3; int endi = 3;
      for(dj=0; dj<endj; dj++){
	for(di=0; di<endi; di++){
	  picbuf[j*endj+dj][i*(endi*4)+(di*4)]=pix_num.c[0];
	  picbuf[j*endj+dj][i*(endi*4)+((di*4)+1)]=pix_num.c[1];
	  picbuf[j*endj+dj][i*(endi*4)+((di*4)+2)]=pix_num.c[2];
	  picbuf[j*endj+dj][i*(endi*4)+((di*4)+3)]=pix_num.c[3];
	}
      }
    }
  }
}

static struct argp_option options[] = {
  {"opacity",   'q' , "val", 0, "Set opacity 0-transparent to 255-opaque"},
  {"output",    'o' , "name", 0, "name of the output png (default: plot.png)"},
  {"zrange",    'z' , "min:max", 0, "Set the min/max z values for plotting"},
  {"pallet",    'c' , "filename", 0, "Use the specified color pallet"},
  {"xml",       'x' , "filename", 0, "Write XML log to filename"},
  {"size",      'g' , "W,H", 0, "Size of plot (in pixels)"},
  {0}
};

#define MAX_NAME 1024

struct arguments{
  char filename[MAX_NAME];
  char output[MAX_NAME];
  char palletname[MAX_NAME];
  char xmlname[MAX_NAME];
  int width, height;
  double llx,lly,llz;
  double ulx,uly,ulz;
  double lrx,lry,lrz;
  int transparency;
  double plotmin,plotmax;
};
 
static error_t parse_opt(int key, char *arg, struct argp_state *state){
  struct arguments *arguments = state->input;

  switch(key){
  case 'q':
    sscanf(arg,"%d",&(arguments->transparency));
    break;
  case 'o':
    strncpy(arguments->output,arg,MAX_NAME);
    break;
  case 'c':
    strncpy(arguments->palletname,arg,MAX_NAME);
    break;
  case 'x':
    strncpy(arguments->xmlname,arg,MAX_NAME);
    break;
  case 'z':
    sscanf(arg,"%lf,%lf",&(arguments->plotmin),&(arguments->plotmax));
    break;
  case 'g':
    sscanf(arg,"%d,%d",&(arguments->width),&(arguments->height));
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

static struct argp argp = {options,parse_opt,"FILENAME","$Id: netcdf2png.c,v 1.12 2007-11-27 01:47:09 luko Exp $"};

int main(int argc, char *argv[]){
  int i, j;
  /* PNG structures */
  png_structp png_ptr;
  png_infop info_ptr;
  png_text text_ptr[1];
  FILE *pngfile;  
  char varname[NC_MAX_NAME];

  /* XML description buffer */
  int xmlsize;
  char *xmlbuf;
  FILE *xmlfile;
  
  /* NetCDF variables */
  int ncid,varid;
  /* Data plotting information/state */
  struct picdimensions picdims;
  struct mergedimensions datadims;
  /* Options */
  struct arguments arguments;
  
  /* Set default options */
  strncpy(arguments.filename,"-",MAX_NAME);
  strncpy(arguments.palletname,"",MAX_NAME);
  strncpy(arguments.xmlname,"",MAX_NAME);
  strncpy(arguments.output,"plot.png",MAX_NAME);
  arguments.plotmin=NAN; arguments.plotmax=NAN;
  arguments.transparency=210;
  arguments.width=1120;
  arguments.height=1120;
  arguments.llx=-25.0;arguments.lly=-25.0;arguments.llz=0.0;
  arguments.ulx=-25.0;arguments.uly=+25.0;arguments.ulz=0.0;
  arguments.lrx=+25.0;arguments.lry=-25.0;arguments.lrz=0.0;

  /* Parse command line */
  argp_parse(&argp,argc,argv,0,0,&arguments);
  
  /* Clean up */
  if(isnan(arguments.plotmin) || isnan(arguments.plotmax)){
    arguments.plotmin=0;
    arguments.plotmax=60.0;
  }
  
  if(strlen(arguments.palletname)==0)
    strncpy(arguments.palletname, "/home/ldm/netcdf2png/colorscales/max_wind.png",MAX_NAME);
  
  /* Make picture description */
  picdims.width=arguments.width;
  picdims.height=arguments.height;
  picdims.transparency=arguments.transparency;
  picdims.plotmin=arguments.plotmin;
  picdims.plotmax=arguments.plotmax;
  
  if(make_pallet_from_png(arguments.palletname,
			  &(picdims.palletsize), &(picdims.pallet))){
    printf("Error reading pallet, giving up\n");
    exit(-1);
  }
  
  /* Open netcdf file to see if we even can plot this at all */
  char junkname[NC_MAX_NAME+1];
  size_t array_start[2];
  size_t array_size[2];

  float **databuf;
  
  /*for obtaining the start time*/
  static size_t start[]={0};
  static size_t count[]={1};
  int starttime;
  //size_t timelen;
  //int status;
  time_t starttime_t;
  struct tm brokendown_time;

  nc_open(arguments.filename,NC_NOWRITE,&ncid);

  //nc_inq_attlen(ncid, NC_GLOBAL, "Time", &timelen);
  //starttime = (int *) malloc(timelen * sizeof(int));
  nc_get_att_int(ncid, NC_GLOBAL, "Time", &starttime);
  starttime_t = (time_t)starttime;
  strftime(datadims.timeformat,16,"%Y%m%d%H%M%S",gmtime_r(&starttime_t,&brokendown_time));

  /* Read position parameters, etc... */
  /* FIXME - This need some error checking */

  nc_inq_dimid(ncid,"Lat",&varid);
  nc_inq_dim(ncid,varid,junkname,&(datadims.num_lats));

  nc_inq_dimid(ncid,"Lon",&varid);
  nc_inq_dim(ncid,varid,junkname,&(datadims.num_lons));
  
  //ok lets make a picture
  /* Set up the PNG stuff */
  pngfile = fopen(arguments.output, "wb");
  if (!pngfile){
    return (ERROR);
  }
  
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png_ptr)
    return (ERROR);
  info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr){
    png_destroy_write_struct(&png_ptr,(png_infopp)NULL);
    return (ERROR);
  }
  if (setjmp(png_jmpbuf(png_ptr))){
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(pngfile);
    return (ERROR);
  }
  png_init_io(png_ptr, pngfile);
  png_set_IHDR(png_ptr, info_ptr, arguments.width, arguments.height,
	       8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  png_write_info(png_ptr, info_ptr);
  
  /* Read data */
  nc_inq_varname(ncid, 0, varname);
  nc_inq_varid(ncid,varname,&varid);
    
  databuf=(float**)malloc(datadims.num_lats*sizeof(float*));
  array_size[0]=1;
  array_size[1]=datadims.num_lons;
  for(i=0;i<datadims.num_lats;i++){
    array_start[0]=i;
    array_start[1]=0;    
    databuf[i]=(float*)malloc(datadims.num_lons*sizeof(float));
    nc_get_vara_float(ncid,varid,array_start,array_size,databuf[i]);
  }
  
  /* Write xml description of what kind of picture we're making */
  write_xml(arguments.filename,&picdims,&datadims,ncid,&xmlbuf,&xmlsize);
  
  /* Optionally write xml comments to a file */
  if(strlen(arguments.xmlname)!=0){
    xmlfile=fopen(arguments.xmlname,"w");
    if(xmlfile){
      fprintf(xmlfile,"%s\n",xmlbuf);
      fclose(xmlfile);
    } else {
      printf("could not write xml file\n");
    }
  }
  
  /* Fill in picture data */
  int bytes_per_pixel=4;
  png_byte** row_pointers;

  /* Allocate picture space */
  row_pointers=(png_byte**)malloc(arguments.height*sizeof(png_byte*));
  for (i = 0; i < arguments.height; i++) {
    row_pointers[i] = (png_byte*)malloc(arguments.width*bytes_per_pixel);
  }

  /* Fill in comments */
  text_ptr[0].compression=PNG_TEXT_COMPRESSION_NONE;
  text_ptr[0].key=strdup("Netcdf2png Parameters");
  text_ptr[0].text=xmlbuf;
  
  png_set_text(png_ptr,info_ptr,text_ptr,1);
  data2picture(databuf,&datadims,(char **)row_pointers,&picdims);
  
  /* Write to the picture */
  png_write_image(png_ptr, row_pointers);
  png_write_info(png_ptr, info_ptr);
  png_write_end(png_ptr, NULL);
  
  return 0;
}
