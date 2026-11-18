#ifdef _NETCDF

#include "netcdfhlp.h"
#include <ncFile.h>
#include <ncDim.h>
#include <ncException.h>
#include <iomanip>
#include <netcdf>
#include <map>
#include <time.h>
#include <sstream>
#include <algorithm>


using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

/////////////////////
// CLASS NETCDFHLP //
/////////////////////

void netcdfhlp::findVarNameFromGridNameExt(char* varname, char* gridname) {

	std::string tmpname(gridname);
	while (tmpname.find_first_of("/\\") != std::string::npos)
	{
		tmpname = tmpname.substr(tmpname.find_first_of("/\\") + 1);
	}

	while (tmpname.find_first_of(".") != std::string::npos)
	{
		tmpname = tmpname.substr(0, tmpname.find_first_of(".")) + tmpname.substr(tmpname.find_first_of(".") + 1);
	}

	strcpy(varname, tmpname.c_str());
	return;

}


void netcdfhlp::findVarNameFromGridName(char* varname, char* gridname) {

	std::string tmpname(gridname);
	while (tmpname.find_first_of("/\\") != std::string::npos)
	{
		tmpname = tmpname.substr(tmpname.find_first_of("/\\")+1);
	}

	if (tmpname.find(".") != std::string::npos)
	{
		tmpname = tmpname.substr(0, tmpname.find("."));
	}

	if (tmpname.find("_timeseries") != std::string::npos)
	{
		tmpname = tmpname.substr(0, tmpname.find("_timeseries"));
	}

	while (tmpname.find_first_of("_-") != std::string::npos)
	{
		tmpname = tmpname.substr(0, tmpname.find_first_of("_-")) + tmpname.substr(tmpname.find_first_of("_-") + 1);
	}

	strcpy(varname, tmpname.c_str());
	return;

}


void netcdfhlp::findFileNameFromGridName(char* filename, char* gridname, bool isSumFile) {

	char *tmpname = (char *)calloc(_stringsize_+1,sizeof(char));
	strcpy(tmpname, gridname);

	// instead of deleting all after the first _ or .
	strcpy(filename, gridname);


	// Append netcdf-file ending, but only if not already present
	int index = std::string(filename).find(".nc");
	if (index<0)
	{
		if (isSumFile) {
			sprintf(filename + strlen(filename), "%s", ".nc");
			//		sprintf(filename, "%s%s", filename, "_sum.nc");
		}
		else {
			sprintf(filename + strlen(filename), "%s", ".nc");
		}
	}
	free(tmpname);

	return;

}

//////////////////////////
// CLASS NETCDFGEOMETRY //
//////////////////////////

netcdfGeometry::netcdfGeometry() {
	xll_ = 0.;
	yll_ = 0.;
	cellsize_ = 1.;
}

netcdfGeometry::~netcdfGeometry() {

}

void netcdfGeometry::setGeometry(double xll, double yll, double cellsize) {
	xll_ = xll;
	yll_ = yll;
	cellsize_ = cellsize;
}

double netcdfGeometry::lat(double x, double y) {
	// TODO
   return 18.;
}

double netcdfGeometry::lon(double x, double y) {
	// TODO
   return 34.;
}

////////////////////////
// CLASS NETCDFWRITER //
////////////////////////

#ifndef SG
netcdfWriter::netcdfWriter(wasim::mystream *Stream) {
	WSstream = Stream;
	has_latgrid_ = false;
	has_longrid_ = false;
}
#else
netcdfWriter::netcdfWriter() {
	has_latgrid_ = false;
	has_longrid_ = false;
}
#endif

netcdfWriter::~netcdfWriter() {
	for(int ii=0; ii<openFiles_.size(); ii++) {
		free(openFiles_[ii]);
	}
 	for(int ii=0; ii<listOfFilesForVars_.size(); ii++) {
 		free(listOfFilesForVars_[ii]);
 	}
 	for(int ii=0; ii<listOfVarsInFile_.size(); ii++) {
 		free(listOfVarsInFile_[ii]);
 	}
}

bool netcdfWriter::fileIsOpen(char* filename) {
	for(int ii=0; ii<openFiles_.size(); ii++) {
		if(strcmp(openFiles_[ii],filename)==0) return true;
	}
	return false;
}

void netcdfWriter::setFileOpen(char* filename) {
	char *newfile = (char *)calloc(_stringsize_+1,sizeof(char));
	strcpy(newfile, filename);
	openFiles_.push_back(newfile);
}

bool netcdfWriter::varExistsInFile(char* filename, char* varname) {

	for(int ii=0; ii<listOfVarsInFile_.size(); ii++) {
		if(strcmp(listOfFilesForVars_[ii],filename)==0 &&
		   strcmp(listOfVarsInFile_[ii],varname)==0) return true;
	}
	return false;
}

void netcdfWriter::setVarInFile(char* filename, char* varname) {
	char *newfile = (char *)calloc(_stringsize_+1,sizeof(char));
	char *newvar = (char *)calloc(_stringsize_+1,sizeof(char));
	strcpy(newfile, filename);
	strcpy(newvar, varname);
	listOfFilesForVars_.push_back(newfile);
	listOfVarsInFile_.push_back(newvar);
}

void netcdfWriter::setLatGrid(double* latgrid) {
	has_latgrid_ = true;
	latgrid_ = latgrid;
}

void netcdfWriter::setLonGrid(double* longrid) {
	has_longrid_ = true;
	longrid_ = longrid;
}

void netcdfWriter::writeTimestep(bool doublemode, double* dataD,
								 float* dataF, char* filename,
								 char* varname, int nrows, int ncols,
								 double nodata,double cellsize,
								 int year, int  month, int  day, int  hour, int  minute)
{

	try
	{
		bool isTimeVariable = false;
		std::string dimTName;
		if (year > 0 && month > 0 && day > 0 && hour >= 0 && minute >= 0)
		{
			isTimeVariable = true;
			dimTName = std::string("time");
		}
		else
		{
			dimTName = std::string("t");
		}

		double actTimeValue = 0;
		if (isTimeVariable)
		{
			// calculate minutes since "1970-01-01 01:00:00.0 +0100"
			time_t acttime;
			tm* tm1 = new(tm);
			bool bApplyRule = false;
			tm1->tm_hour = hour;
//			tm1->tm_min = useMinutes ? minute : 0; // 0..59
			tm1->tm_min = minute; // 0..59
			tm1->tm_sec = 0; // 0..59
			tm1->tm_mon = month - 1;  // 0..11
			tm1->tm_year = year - 1900; // e.g. 2010 = 110
			tm1->tm_mday = day;
			tm1->tm_wday = 0;
			tm1->tm_yday = 0;
			tm1->tm_isdst = 0;
#if defined WIN64
			acttime = _mkgmtime(tm1); // Windows 64
#elif defined WIN32
			acttime = _mkgmtime(tm1); // Windows 32
#elif defined SG
			acttime = mktime(tm1); // disadvantag here: Daylight saving time is used
#else
			acttime = timegm(tm1); // Linux
#endif
			//acttime = mktime(tm1);
			actTimeValue = (double)(acttime) / 60.0 - 60.0; // since the reference time starts at 1970-01-01 01:00:00 and not at 00:00:00, we need to subtract 60 minutes from the counter
		}

		if (!fileIsOpen(filename)) {

			// Create a new netCDF file and write header, attributes and
			// lat/lon-grids if necessary

			setFileOpen(filename);
			setVarInFile(filename, varname);
			// new C++ API uses exceptions, so we do not need to check the file after declaring
			NcFile dataFile(filename, NcFile::replace, NcFile::nc4);

			NcDim dimX = dataFile.addDim("x", ncols);
			NcDim dimY = dataFile.addDim("y", nrows);
			NcDim dimT = dataFile.addDim(dimTName);  //adds an unlimited dimension

			double *xx = new double[ncols];
			for (int i = 0; i < ncols; i++) { xx[i] = xll_ + (double)i*cellsize_; }
			double *yy = new double[nrows];
			double yol = yll_ + (nrows - 1)*cellsize_;
			for (int j = 0; j < nrows; j++) { yy[j] = yol - (double)j*cellsize_; }
			double tt = 0.;

			NcVar xVar = dataFile.addVar("x", ncDouble, dimX);
			xVar.setCompression(false, true, 3);
//			xVar.putAtt("units", "m East of reference point");
			xVar.putAtt("standard_name", "projection_x_coordinate");
//			xVar.putAtt("long_name", "x coordinate according to CH1903");
			xVar.putAtt("units", "m");
			xVar.putAtt("axis", "X");
			xVar.putAtt("missing_value", ncDouble, nodata);
			xVar.putAtt("_FillValue", ncDouble, nodata);
			xVar.putVar(xx);

			NcVar yVar = dataFile.addVar("y", ncDouble, dimY);
			yVar.setCompression(false, true, 3);
//			yVar.putAtt("units", "m North of reference point");
			yVar.putAtt("standard_name", "projection_y_coordinate");
//			yVar.putAtt("long_name", "y coordinate according to CH1903");
			yVar.putAtt("units", "m");
			yVar.putAtt("axis", "Y");
			yVar.putAtt("missing_value", ncDouble, nodata);
			yVar.putAtt("_FillValue", ncDouble, nodata);
			yVar.putVar(yy);

			NcVar timeVar = dataFile.addVar(dimTName, ncDouble, dimT);
//			vector<size_t> timesize;
//			timesize.push_back(0);
//			timeVar.putVar(timesize,&tt);
			if (isTimeVariable)
			{
				timeVar.putAtt("standard_name", "time");
				timeVar.putAtt("long_name", "time");
				timeVar.putAtt("axis", "T");
				timeVar.putAtt("units", "minutes since 1970-01-01 01:00:00.0 +0100");
			}

			if (has_latgrid_) {
				std::vector<NcDim> vdims;
				vdims.push_back(dimX);
				vdims.push_back(dimY);
				NcVar latVar = dataFile.addVar("lat", ncDouble, vdims);
				latVar.setCompression(false, true, 3);
				latVar.putAtt("units", "degrees north");
				vector<size_t> vstart, vcount;
				vstart.push_back(0);
				vstart.push_back(0);
				vcount.push_back(ncols);
				vcount.push_back(nrows);
				latVar.putVar(vstart, vcount, latgrid_);
			}

			if (has_longrid_) {
				std::vector<NcDim> vdims;
				vdims.push_back(dimX);
				vdims.push_back(dimY);
				NcVar lonVar = dataFile.addVar("lon", ncDouble, vdims);
				lonVar.setCompression(false, true, 3);
				lonVar.putAtt("units", "degrees east");
				vector<size_t> vstart, vcount;
				vstart.push_back(0);
				vstart.push_back(0);
				vcount.push_back(ncols);
				vcount.push_back(nrows);
				lonVar.putVar(vstart, vcount, longrid_);
			}

			NcVar myVar;
			std::vector<NcDim> vdims;
			vdims.push_back(dimT);
			vdims.push_back(dimX);
			vdims.push_back(dimY);

			if (doublemode)
			{
				myVar = dataFile.addVar(varname, ncDouble, vdims);
				myVar.putAtt("missing_value", ncDouble, nodata);
				myVar.putAtt("_FillValue", ncDouble, nodata);
// test how to rename a variable in an existing file
// 			myVar.rename("xyz");
			}
			else
			{
				myVar = dataFile.addVar(varname, ncFloat, vdims);
				myVar.putAtt("missing_value", ncFloat, nodata);
				myVar.putAtt("_FillValue", ncFloat, nodata);
// test how to rename a variable in an existing file
//				myVar.rename("xyz");
			}

			myVar.setCompression(false, true, 3);
			myVar.putAtt("default_cellsize", ncDouble, cellsize);

			// Write the first time step
			vector<size_t> vstart, vcount;
			vstart.push_back(0);
			vcount.push_back(1);
			int newRecordId = 0;
			if (isTimeVariable)
			{
				timeVar.putVar(vstart, vcount, &actTimeValue);
			}
			else
			{
				timeVar.putVar(vstart, vcount, &newRecordId);
			}

			vstart.push_back(0);
			vstart.push_back(0);
			vcount.push_back(ncols);
			vcount.push_back(nrows);
			if (doublemode)
			{
				myVar.putVar(vstart, vcount, dataD);
			}
			else
			{
				myVar.putVar(vstart, vcount, dataF);
			}

			delete[] xx;
			delete[] yy;

			dataFile.close();

		} else if (!varExistsInFile(filename, varname)) {

			setVarInFile(filename, varname);
			// Add a variable to an existing file

			NcFile dataFile(filename, NcFile::write, NcFile::nc4);


// test how to rename variable in an existing file
//			NcVar myVar1 = dataFile.getVar("xyz");
//			myVar1.rename("zyx");

			NcDim dimT = dataFile.getDim(dimTName);
			NcDim dimX = dataFile.getDim("x");
			NcDim dimY = dataFile.getDim("y");

			NcVar myVar;
			std::vector<NcDim> vdims;
			vdims.push_back(dimT);
			vdims.push_back(dimX);
			vdims.push_back(dimY);

			if (doublemode)
			{
				myVar = dataFile.addVar(varname, ncDouble, vdims);
				myVar.putAtt("missing_value", ncDouble, nodata);
				myVar.putAtt("_FillValue", ncDouble, nodata);
			}
			else
			{
				myVar = dataFile.addVar(varname, ncFloat, vdims);
				myVar.putAtt("missing_value", ncFloat, nodata);
				myVar.putAtt("_FillValue", ncFloat, nodata);
			}

			myVar.setCompression(false, true, 3);
			myVar.putAtt("default_cellsize", ncDouble, cellsize);

			// Write the first time step
			vector<size_t> vstart, vcount;
			vstart.push_back(0);
			vcount.push_back(1);
			NcVar timeVar = dataFile.getVar(dimTName);
			int newRecordId = 0;
			if (isTimeVariable)
			{
				timeVar.putVar(vstart, vcount, &actTimeValue);
			}
			else
			{
				timeVar.putVar(vstart, vcount, &newRecordId);
			}

			vstart.push_back(0);
			vstart.push_back(0);
			vcount.push_back(ncols);
			vcount.push_back(nrows);
			if (doublemode)
			{
				myVar.putVar(vstart, vcount, dataD);
			}
			else
			{
				myVar.putVar(vstart, vcount, dataF);
			}

			dataFile.close();

		} else {

			// Append a new time step to an existing netCDF file

			NcFile dataFile(filename, NcFile::write, NcFile::nc4);

			int newRecordId;
			// how many grids are already contained in the nc file? This is the unlimited dimension (or the time step dimension, the one that changes slowest)
			newRecordId = dataFile.getDim(dimTName).getSize();
			NcVar timeVar = dataFile.getVar(dimTName);

			// creating empty vectors holding start and count elements for each dimension
			vector<size_t> vstart;
			vector<size_t> vcount;

			// first we use the vectors for the time dimension which is one, so only one element has to be created in each vector
			vstart.push_back(newRecordId);
			vcount.push_back(1); // only one new record is added to the time variable at the end (at location newRecordId)

			if (isTimeVariable)
			{
				timeVar.putVar(vstart, vcount, &actTimeValue);
			}
			else
			{
				timeVar.putVar(vstart, vcount, &newRecordId);
			}

			// now writing the grid into the three dimensional data variable
			NcVar myVar = dataFile.getVar(varname);
			// so we need to add two more elements for dimensions x and y to the vectors
			vstart.push_back(0); // columns start at element 0
			vstart.push_back(0); // rows start at element 0
			vcount.push_back(ncols); // all columns should be written
			vcount.push_back(nrows); // all rows as well
			// write depending on mode either double or float values
			if (doublemode)
			{
				myVar.putVar(vstart, vcount, dataD);
			}
			else
			{
				myVar.putVar(vstart, vcount, dataF);
			}


			dataFile.close();

		}
	}
	catch (NcException &e)
	{
#ifndef SG
		(*WSstream) << "FATAL: error creating "
			<< filename
			<< " for writing in netCDF: "
			<< e.what()
			<< "[end-entry]";
#endif
		throw(26);
	}

}


////////////////////////
// CLASS NETCDFREADER //
////////////////////////

#ifndef SG
netcdfReader::netcdfReader(wasim::mystream *Stream):
isTimeVariable(false)
{
	WSstream = Stream;
	//dataFile = NULL;
	mapDataFiles.clear();
}
#else
netcdfReader::netcdfReader():
isTimeVariable(false)
{
	//dataFile = NULL;
	mapDataFiles.clear();
}
#endif

netcdfReader::~netcdfReader() {
	//dataFile->close();
	for (std::map<std::string, netCDF::NcFile*>::iterator iter = mapDataFiles.begin(); iter != mapDataFiles.end(); iter++)
	{
		iter->second->close();
	}
	mapDataFiles.clear();
}

void netcdfReader::readTimestep(bool doublemode, double* &stackbufD, char* filename,char *pvarname,
					int &nc_rows,int &nc_cols,double &nc_xll,double &nc_yll,double &nc_cellsize,
					double &nc_nodata,int &nc_levels,int actLevel,int year, int month , int day , int hour , int minute)
{

	// TODO find out if the netcdf contains double or single precision
	// and make Wasim act correctly in either case

	char *varname = (char *)calloc(_stringsize_+1,sizeof(char));
	char *layername = (char *)calloc(_stringsize_+1,sizeof(char));
	int layer = -1;

	double actTimeValue=0;
	if (year > 0 && month > 0 && day > 0 && hour >= 0 && minute >= 0)
	{
		// calculate minutes since "1970-01-01 01:00:00.0 +0100"
		time_t acttime;
		tm *tm1 = new(tm);
		bool bApplyRule = false;

//		tm1->tm_hour = (int)((60.0 * hour + minutes_to_interval_end - timestep) / 60.0); // 0...23 --> if hour = 24, then (60*24 + 0 - 1440)/60 = 0, should usually be 0 for daily time steps.

		tm1->tm_hour = hour;
		//			tm1.tm_min = useMinutes ? minute : 0; // 0..59
		tm1->tm_min = minute; // 0..59
		tm1->tm_sec = 0; // 0..59
		tm1->tm_mon = month - 1;  // 0..11
		tm1->tm_year = year - 1900; // e.g. 2010 = 110
		tm1->tm_mday = day;
		tm1->tm_wday = 0;
		tm1->tm_yday = 0;
		tm1->tm_isdst = 0;


#if defined WIN64
			acttime = _mkgmtime(tm1); // WINDOWS 64
#elif defined WIN32
			acttime = _mkgmtime(tm1); // WINDOWS 32
#elif defined (SG)
			acttime = mktime(tm1); // disadvantag here: Daylight saving time is used
#else
			acttime = timegm(tm1);  // Linux
#endif

		//acttime = mktime(&tm1);

//		actTimeValue = (double)(acttime) / 60.0 + 60.0; // since mktime returns seconds since 1070-01-01 00:00:00 and we need CET, we add 60 minutes
		// A date/time provided in UTC wouldn't have to be changed here. But the date and time provided in WaSiM runs at the moment,
		// are provided in CET (i.e. UTC +1.0). That means, that the acttime contains a value that is actually 1 hour or 3600 seconds or 60 minutes too large, bacause the netCDF we use
		// will interpret the minutes since 1970-01-01 01:00:00 (and not 00:00:00), i.e. the counter must always be set 60 minutes smaller for those specific netCDFs
		actTimeValue = (double)(acttime) / 60.0 - 60.0;

		delete tm1;
	}

	NcVar data;
	NcVar timeVar;

	strcpy(varname, filename);
	strcpy(layername, filename);

	int varnameindex = index(filename, "<");
	int layerindex = index(filename, ">");

	if (layerindex != 0 && varnameindex != 0) {
		// We have a < and a > in the filename, giving varname and timestep/layer.
		// We separate the filename in front of the <> from the
		// variable name and the timestep/layer.
		substr(filename, 0, varnameindex - 1);
		substr(varname, varnameindex + 1, layerindex - varnameindex - 2);
		substr(layername, layerindex + 1, 9);
		if (strlen(layername) > 0)
		{
			layer = atoi(layername);
		}
		else
		{
			layer = -1;
		}
	}
	else
	{
		// either only a varname or nether varname nor layername is present, so layer is set to -1
		layer = -1;
	}
	if ((actLevel >= 0)&&(layer==-1))
	{
		layer = actLevel; // return the layer that is provided as optional parameter, otherwise the first layer (with index 0, the default of the parameter)
	}

	try
	{
		// new C++ API
		// look for the file if it is already open:
		std::map<std::string, netCDF::NcFile*>::iterator mapiter;
		mapiter = mapDataFiles.find(filename);
		if (mapiter == mapDataFiles.end())
		{ // new file, so a new map entry has to be added
			if (mapDataFiles.size() >= _MAX_OPEN_NETCDF_FILES_)
			{
#ifndef SG
				(*WSstream) << "FATAL: more than "<<_MAX_OPEN_NETCDF_FILES_<<" different netCDF files can not be hold open at the same time. If you read time series of input grids from netCDF, please put all the single nc files into one large file or into some groups of files" << "[end-entry]";
#endif
				throw(1100);
			}
			mapDataFiles.insert(pair<std::string, netCDF::NcFile*>(filename, new NcFile(filename, NcFile::read)));
			mapiter = mapDataFiles.find(filename);
		}

		dataFile = mapiter->second;

		if (dataFile == NULL)
		{
			dataFile = new NcFile(filename, NcFile::read);
		}
		//NcFile dataFile(filename, NcFile::read);


		// new C++ API
		isTimeVariable = false;
		int num_dims = dataFile->getDimCount();
		std::multimap<std::string, netCDF::NcDim> mm_dims = dataFile->getDims();
		for (std::multimap<std::string, netCDF::NcDim>::iterator iter = mm_dims.begin(); iter != mm_dims.end(); iter++)
		{
			if (iter->first.find("x") != std::string::npos)		{
				nc_cols = iter->second.getSize();

			}
			else if (iter->first.find("y") != std::string::npos) {
				nc_rows = iter->second.getSize();
			}
			else if (iter->first.find("time") != std::string::npos) {
				nc_levels = iter->second.getSize();
				isTimeVariable = true;
			}
			else if (iter->first.find("t") != std::string::npos) {
				nc_levels = iter->second.getSize();
			}

		}

		if (isTimeVariable && layer==-1) // get layer number from actTimeValue instead of from external layer
		{
			timeVar = dataFile->getVar("time");
			long timeSize = nc_levels;
			std::vector<double>* timeVector = new std::vector<double>(nc_levels,0);
			int s = timeVector->size();
			timeVar.getVar(timeVector->data());

			std::vector<double>::iterator iterlayer = std::find(timeVector->begin(), timeVector->end(), actTimeValue);
			if (iterlayer != timeVector->end())
			{
				layer = iterlayer - timeVector->begin();
			}
			else
			{
				layer = -1;
			}
			int ys = 0;
		}

		if (layer >=0 && varnameindex != 0) {
			// layer = atoi(layername);
			data = dataFile->getVar(varname);
		}
		else {
			// if no varname is given, it can still be a netCDF-grid from stateini or a standard grid that has a specific name
			// so we look for a varname that equals either
			// a) the gridname without extension
			// b) the gridname with extension
			// c) any name (will read then layer 0 from first variable that is neither "x" nor "y" nor "t")
			bool varFound = false;

			netcdfhlp::findVarNameFromGridName(varname, filename);
			data = dataFile->getVar(varname);

			if (data.isNull()) {
				netcdfhlp::findVarNameFromGridNameExt(varname, filename);
				data = dataFile->getVar(varname);
				if (data.isNull()) {

					int varcnt = dataFile->getVarCount();
					std::multimap<std::string, NcVar> vars = dataFile->getVars();

					for (std::multimap<std::string, NcVar>::iterator iter = vars.begin(); iter != vars.end(); iter++)
					{
						if (
							(iter->first.compare("x") != 0) &&
							(iter->first.compare("y") != 0) &&
							(iter->first.compare("t") != 0) &&
							(iter->first.compare("lat") != 0) &&
							(iter->first.compare("lon") != 0) &&
							(iter->first.compare("crs") != 0) &&
							(iter->first.compare("time") != 0)
						   )
						{
							sprintf(varname, "%s", iter->first.c_str());
							varFound = true;
							iter = --vars.end(); // break the search here, just in case there are other variables. We take the first one.
#ifndef SG
							(*WSstream) << "INFO: netcdf input file names should be provided as filename<varname>timestep (example: \"temperature.nc<tair>12\")." << "[end-entry]";
#endif
							data = dataFile->getVar(varname);
						}
					}

				}
				else
				{
					varFound = true;
				}
			}
			else
			{
				varFound = true;
			}
		}
		if (data.isNull()) {
#ifndef SG
			(*WSstream) << "FATAL: cannot find variable "
				<< "\"" << varname << "\""
				<< " in file "
				<< "\"" << filename << "\""
				<< " for reading from netCDF."
				<< "[end-entry]";
#endif
			throw(26);
		}

		double scale_factor = 1.;
		double add_offset = 0;

		// new C++ API
		nc_cellsize = -1.0; // if found as attribute, it will be taken from there, otherwise it must be calculated from x and y coordinates, which is not working for one cell models. So it is recommended to put an attribute default_cellsize into the data layer
		std::map<std::string, netCDF::NcVarAtt> mm_atts = data.getAtts();
		for (map<string, NcVarAtt>::iterator iter = mm_atts.begin(); iter != mm_atts.end(); iter++)
		{
			if (iter->first.find("scale_factor") != std::string::npos)
			{
				iter->second.getValues(&scale_factor);
			}
			else if (iter->first.find("add_offset") != std::string::npos)
			{
				iter->second.getValues(&add_offset);
			}
			else if (iter->first.find("missing_value") != std::string::npos)
			{
				iter->second.getValues(&nc_nodata);
			}
			else if (iter->first.find("_FillValue") != std::string::npos) // as an alternative to missing_value
			{
				iter->second.getValues(&nc_nodata);
			}
			else if (iter->first.find("default_cellsize") != std::string::npos)
			{
				iter->second.getValues(&nc_cellsize);
			}
		}

		netCDF::NcDim dim1 = data.getDim(1);
		netCDF::NcDim dim2 = data.getDim(2);

		// getting yll and info on how rows and cols are stored
		bool bottomup = true;
		double *ybuf = new double[nc_rows];
		std::vector<size_t> vstart, vcount;
		vstart.push_back(0);
		vcount.push_back(nc_rows);
		NcVar varYDim;
		varYDim = dataFile->getVar("y");
		if (!varYDim.isNull()) {
			varYDim.getVar(vstart, vcount, ybuf);
			if ((nc_rows > 1) && (ybuf[1] < ybuf[0]))
			{
				bottomup = false;
				nc_yll = ybuf[nc_rows - 1];
				if (nc_cellsize < 0)
				{
					if (nc_rows > 1)
					{
						nc_cellsize = ybuf[0] - ybuf[1];
					}
					else
					{
						nc_cellsize = 10; // default value, 10m
					}
				}
			}
			else
			{
				nc_yll = ybuf[0];
				if (nc_cellsize < 0)
				{
					if (nc_rows > 1)
					{
						nc_cellsize = ybuf[1] - ybuf[0];
					}
					else
					{
						nc_cellsize = 10; // default value, 10m for a one cell model (if only netCDF is read in, some other method for storing the cellsize must be implemented
					}
				}
			}
		}
		delete[] ybuf;
		vstart.clear();
		vcount.clear();

		// getting xll
		double *xbuf = new double[nc_cols];
		vstart.push_back(0);
		vcount.push_back(nc_cols);
		NcVar varXDim;
		varXDim = dataFile->getVar("x");
		if (!varXDim.isNull()) {
			varXDim.getVar(vstart, vcount, xbuf);
			nc_xll = xbuf[0];
		}
		delete[] xbuf;
		vstart.clear();
		vcount.clear();


		// Go to the correct timestep/layer and read the data
		// new C++ API
		stackbufD = new double[(size_t)nc_cols * nc_rows];
		vstart.push_back(layer);
		vstart.push_back(0);
		vstart.push_back(0);
		vcount.push_back(1);
		if (dim1.getName()==std::string("x"))
		{	// first dimension is "x"
			vcount.push_back(nc_cols);
			vcount.push_back(nc_rows);
		}
		else
		{	// first dimension is "y"
			vcount.push_back(nc_rows);
			vcount.push_back(nc_cols);
		}
		if (data.getType() == netCDF::ncDouble)
		{
			data.getVar(vstart, vcount, stackbufD);
			for (size_t ii = 0; ii < (size_t)nc_cols * nc_rows; ++ii)
			{
				stackbufD[ii] = stackbufD[ii] * scale_factor + add_offset;
			}
		}
		else
		{
			float *stackbufF = new float[(size_t)nc_cols * nc_rows];
			data.getVar(vstart, vcount, stackbufF);
			for (size_t ii = 0; ii < (size_t)nc_cols * nc_rows; ++ii)
			{
				stackbufD[ii] = double(stackbufF[ii]) * scale_factor + add_offset;
			}
			delete[] stackbufF;
		}

		if (dim2.getName() == "x")
		{ // transposing rows and columns because x and y coordinates where not in x-y-order but in y-x-order
			double *temp = new double[(size_t)nc_cols * nc_rows];
			memmove(temp, stackbufD, (size_t)nc_cols * nc_rows * sizeof(double));
			size_t counter = 0;
			if (bottomup == true)
			{
				for (long i = (long)nc_cols * (nc_rows - 1); i < (long)nc_cols * nc_rows; i++)
				{
					for (long j = i; j >= 0; j -= nc_cols)
					{
						stackbufD[counter] = temp[j];
						counter++;
					}
				}
			}
			else
			{
				for (int i = 0; i < nc_cols; i++)
				{
					for (long j = i; j < (long)nc_rows * nc_cols; j += nc_cols)
					{
						stackbufD[counter] = temp[j];
						counter++;
					}
				}
			}
			delete[] temp;
		}
		else
		{
			if (dim1.getName() == "x")
			{
				if (bottomup == true)
				{
					double *temp = new double[(size_t)nc_cols * nc_rows];
					memmove(temp, stackbufD, (size_t)nc_cols * nc_rows * sizeof(double));
					size_t counter = 0;
					for (int i = 0; i < nc_cols; i++)
					{
						for (int j = (i+1)*nc_rows-1; j >= i*nc_rows; j--)
						{
							stackbufD[counter] = temp[j];
							counter++;
						}
					}
					delete[] temp;
				}
				// no topdown required since this is the default anyway
			}
		}


		//dataFile->close();

		free(layername);
		if (pvarname) {
			strcpy(pvarname,varname);
		}
		free(varname);
	}
	catch (NcException &e)
	{
#ifndef SG
		(*WSstream) << "FATAL: exception occured when reading from netCDF file "
			<< filename << ". Internal message is: " << e.what()
			<< "[end-entry]";
#endif
		throw(26);
	}
}


std::string netcdfReader::readTimeValue(char* filename,int actLevel)
{
	NcVar timeVar;
	NcVar data;

	if (this->isTimeVariable==true)
	{

		std::map<std::string, netCDF::NcFile*>::iterator mapiter;
		mapiter = mapDataFiles.find(filename);
		if (mapiter == mapDataFiles.end())
		{ // new file, so a new map entry has to be added
			if (mapDataFiles.size() >= _MAX_OPEN_NETCDF_FILES_)
			{
	#ifndef SG
				(*WSstream) << "FATAL: more than "<<_MAX_OPEN_NETCDF_FILES_<<" different netCDF files can not be hold open at the same time. If you read time series of input grids from netCDF, please put all the single nc files into one large file or into some groups of files" << "[end-entry]";
	#endif
				throw(1100);
			}
			mapDataFiles.insert(pair<std::string, netCDF::NcFile*>(filename, new NcFile(filename, NcFile::read)));
			mapiter = mapDataFiles.find(filename);
		}
		dataFile = mapiter->second;

		timeVar = dataFile->getVar("time");
		long timeSize = dataFile->getDim("time").getSize();

		std::vector<double>* timeVector = new std::vector<double>(timeSize,0);
		int s = timeVector->size();
		timeVar.getVar(timeVector->data());

//		time_t dActTime = (*timeVector)[actLevel]*60.0+3600;
		time_t dActTime = (*timeVector)[actLevel]*60.0;     // noc corrections will be made on time. If used in SG, Daylight Saving Time might be used instead for display

		tm *tmstart;
		tmstart = gmtime(&dActTime); // gmtime returns a pointer to a structure, no need to allocate that memory
		if (tmstart == NULL)
		{
#ifndef SG
			(*WSstream) << "FATAL: actual time value in netCDF is invalid (either before 1970-01-01 or after 3000-12-31)" << "[end-entry]";
#endif
			throw(1101);
		}
		int iYear = tmstart->tm_year + 1900;
		int iMonth = tmstart->tm_mon + 1;
		int iDay = tmstart->tm_mday;
		int iHour = tmstart->tm_hour;
		int iMinute = tmstart->tm_min;


		std:stringstream ss;
		ss << setw(4) << setfill('0') << iYear << "-"
		   << setw(2) << setfill('0') << iMonth  << "-"
		   << setw(2) << setfill('0') << iDay << endl;
		std::string sout1;
		ss >> sout1;

		ss.clear();
		ss << setw(2) << setfill('0') << iHour  << ":"
		   << setw(2) << setfill('0') << iMinute  << ":00" << endl;
		std::string sout2;
		ss>>sout2;

		sout1.append(" ");
		sout1.append(sout2);
		return sout1;
	}
	else
	{
		return std::string("");
	}

}

#endif
