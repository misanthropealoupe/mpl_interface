#include "mpl_interface.hpp"
#include "dedisp-contain.hpp"
#include "numpy/arrayobject.h"
#include <string>
#include <tuple>
#include <iostream>
#include <random>

#define NOARGS PyTuple_Pack(0)
#define DISPCONST 4.149e3

double disp_delay(double f, double dm)
{
	return dm * DISPCONST * (1.0f/(f * f));
}

double disp_delta(double f, double dm, search_data sd){
	return disp_delay(f,dm) - disp_delay(sd.f0,dm);
}

double *test_data(int r, int c)
{
	double *ret = new double[r*c];
	for(int i = 0; i < r; i++){
		for(int j = 0; j < c; j++){
			ret[i*c + j] = (double) (i*c + j);
		}
	}
	return ret;
}

PyObject *getPython2DArray(double* ar, int nrow, int ncol)
{
	PyObject *ret = PyList_New(nrow);
	for(int i = 0; i < nrow; i++){
		PyObject *thisrow = PyList_New(ncol);
		for(int j = 0; j < ncol; j++){
			PyList_SetItem(thisrow,j,PyFloat_FromDouble(ar[i*ncol + j]));
		}
		PyList_SetItem(ret,i,thisrow);
	}
	return ret;
}

PyObject *getPythonArray(double* ar, int n)
{
	PyObject *ret = PyList_New(n);
	for(int i = 0; i < n; i++){
		PyList_SetItem(ret,i,PyFloat_FromDouble(ar[i]));
	}
	return ret;
}

PyObject *getPythonLinspace(double a, double b, int n)
{
	PyObject *ret = PyList_New(n);
	double delta = (a - b)/((double) (n + 1));
	for(int i = 0; i < n; i++){
		PyList_SetItem(ret,i,PyFloat_FromDouble(delta*(1 + i)));
	}
	return ret;
}

void imshow_show(double *ar, int nrow, int ncol, PyObject* mpl)
{
	PyObject* imshowFun = PyObject_GetAttrString(mpl,(char*)"imshow");
	PyObject* a = getPython2DArray(ar,nrow,ncol);
	PyObject* args = PyTuple_Pack(1,a);
	PyObject* noargs = PyTuple_Pack(0); 
	PyObject* result = PyObject_CallObject(imshowFun,args);
	PyObject* showFun = PyObject_GetAttrString(mpl,(char*)"show");
	PyObject_CallObject(showFun,noargs);
}

void imshow_save_simple(double *ar, int nrow, int ncol, std::string fname)
{
	if(Py_IsInitialized() == 0){
		Py_Initialize();
	}
	PyObject* mplBaseString = PyString_FromString((char*)"matplotlib");
	PyObject* mplBase = PyImport_Import(mplBaseString);
	call(mplBase, "use", PyTuple_Pack(1,PyString_FromString("Agg")));
	PyObject* mplString = PyString_FromString((char*)"matplotlib.pyplot");
	PyObject* mpl = PyImport_Import(mplString);

	imshow_save(ar,nrow,ncol,mpl,fname,true);
	//Py_Finalize();
}

void imshow_save(double *ar, int nrow, int ncol, PyObject* mpl, std::string fname, bool colorbar)
{
	PyObject* cbarfun = PyObject_GetAttrString(mpl, (char*)"colorbar");
	PyObject* axFun = PyObject_GetAttrString(mpl,(char*)"gca");
	PyObject* ax = PyObject_CallObject(axFun,NOARGS);
	PyObject* imshowFun = PyObject_GetAttrString(ax,(char*)"imshow");
	PyObject* a = getPython2DArray(ar,nrow,ncol);
	PyObject* args = PyTuple_Pack(1,a);
	PyObject* saveargs = PyTuple_Pack(1,PyString_FromString(fname.c_str()));
	PyObject* result = PyObject_CallObject(imshowFun,args);

	if(colorbar){	
		PyObject_CallObject(cbarfun,PyTuple_Pack(1,result));
	}

	PyObject* saveFun = PyObject_GetAttrString(mpl,(char*)"savefig");
	std::cout << "saving plot\n";
	PyObject_CallObject(saveFun,saveargs);
	PyObject* closeFun = PyObject_GetAttrString(mpl,(char*)"close");
	PyObject_CallObject(closeFun,NOARGS);
}

double* zeros_d(int n)
{
	double* ret = new double[n];
	for(int i = 0; i < n; i++){
		ret[i] = 0.0f;
	}
	return ret;
}

double* fill_d(int n, double fillval)
{
	double* ret = new double[n];
	for(int i = 0; i < n; i++){
		ret[i] = fillval;
	}
	return ret;
}

double *randArray(double mean, double std, int n)
{
	std::random_device rd;
	std::mt19937 e2(rd());
	std::normal_distribution<> gaussian(mean, std);
	double* ret = new double[n];
	for(int i = 0; i < n; i++){
		ret[i] = gaussian(e2);
	}
	return ret;
}

PyObject *callKey(PyObject* target, std::string fname, PyObject *args, PyObject *keywords)
{
	PyObject* fun = PyObject_GetAttrString(target,fname.c_str());
	return PyObject_Call(fun,args,keywords);
}

PyObject *call(PyObject* target, std::string fname, PyObject *args)
{
	PyObject* fun = PyObject_GetAttrString(target,fname.c_str());
	return PyObject_CallObject(fun,args);
}

std::tuple<int,double*> dedisperse(double *ft, double dm, double fillval, search_data sd)
{
	double dt = sd.t_frame/((double) sd.nt);
	int max_delay = (int) (disp_delta(sd.f1,sd.maxdm,sd)/dt);
	int nt_ext = sd.nt + max_delay;
	double* ft_remap = fill_d(sd.nf * nt_ext,fillval);
	int delay;
	double df = (sd.f1 - sd.f0)/((double) sd.nf);
	double freq;
	for(int i = 0; i < sd.nf; i++){
		freq = sd.f0 + i * df;
		delay = (int) (disp_delta(freq,dm,sd)/dt);
		// printf("delay %i\n",delay);
		// printf("copy %i\n",sd.nt - delay);
		for(int j = 0; j < sd.nt; j++){
			ft_remap[i * nt_ext + j + (max_delay - delay)] = ft[i * sd.nt + j];
		}
	}
	return std::make_tuple(max_delay,ft_remap);
}

#define ARLEN 50
PyObject *getNumpyTickLabels(PyArrayObject* locs, int indLen, int decimalLen, double offset, double delta)
{
	int nlabels = PyArray_DIMS(locs)[0];
	PyObject* labels = PyList_New(nlabels);
	char* buf = new char[ARLEN];
	char* smallbuf;
	char* formatString1 = new char[ARLEN];
	int integerLen;
	int strlen;
	int fmlen = sprintf(formatString1,"%%.%if", decimalLen);
	//printf("fmlen %i\n",fmlen);
	char* formatString = new char[fmlen];
	strncpy(formatString,formatString1,fmlen);

	//std::cout << formatString << "\n";
	for(int i = 0; i < nlabels; i++){
		double tickval = (double) ((double*) PyArray_GETPTR1(locs,i))[0];
		//std::cout << tickval << "\n";
		tickval = tickval * delta + offset;
		strlen = sprintf(buf, "%.2f", tickval);
		smallbuf = new char[strlen];
		strncpy(smallbuf,buf,strlen);
		//printf("%s\n", smallbuf);
		//std::cout << smallbuf << "\n";
		PyList_SetItem(labels,i,PyString_FromString(smallbuf));
		delete[] smallbuf;
	}
	delete[] buf;
	delete[] formatString;
	delete[] formatString1;
	return labels;
}

void plot_save(double* data, int len, std::string fname, PyObject *mpl)
{
	//PyObject *keywords = Py_BuildValue("{s:O}", "dpi", Py_BuildValue("i", 300));
	//callKey(mpl,"figure",NOARGS,keywords);
	call(mpl,"figure",NOARGS);
	call(mpl,"plot", PyTuple_Pack(1,getPythonArray(data,len)));
	call(mpl, "savefig", PyTuple_Pack(1,PyString_FromString(fname.c_str())));
	call(mpl, "close", NOARGS);
}

void plot_save_simple(double* data, int len, std::string fname)
{
	if(Py_IsInitialized() == 0){
		Py_Initialize();
	}	
	PyObject* mplBaseString = PyString_FromString((char*)"matplotlib");
	PyObject* mplBase = PyImport_Import(mplBaseString);
	call(mplBase, "use", PyTuple_Pack(1,PyString_FromString("Agg")));
	PyObject* mplString = PyString_FromString((char*)"matplotlib.pyplot");
	PyObject* mpl = PyImport_Import(mplString);
	plot_save(data,len,fname,mpl);
	//Py_Finalize();
}

void multipanel_transient_save(double *ft, double *dmt, int ntdm, PyObject *mpl, trigger_data td, search_data sd)
{
	double ddm = sd.maxdm/((double) sd.ndm);
	std::tuple<int,double*>tp = dedisperse(ft,td.dm_ind * ddm,-1.0f,sd);
	double* ft_rm = std::get<1>(tp);
	int max_delay  = std::get<0>(tp);
	int ntf = sd.nt + max_delay;
	// double* ft_rm = ft;
	// int ntf = sd.nt;
	PyObject *ft_list = getPython2DArray(ft_rm, sd.nf, ntf);
	PyObject *dmt_list = getPython2DArray(dmt, sd.ndm, ntdm);

	double* int_intensity = zeros_d(sd.nf);
	int t0 = td.t_ind;
	for(int i = 0; i < sd.nf; i++){
		//printf("f %i\n",i);
		for(int j = 0; j < td.t_width; j++){
			//printf("\tt %i\n",j);
			int_intensity[i] += ft_rm[i * ntf + t0 + max_delay + j];
		}
	}
	PyObject *intensityList = getPythonArray(int_intensity,sd.nf);
	//SEGFAULT WHATTF?
	PyObject *fig_axes = call(mpl,"subplots", PyTuple_Pack(2,Py_BuildValue("i",3),Py_BuildValue("i",1)));

	PyArrayObject* axes = (PyArrayObject*) PyTuple_GetItem(fig_axes,1);
	PyObject* axft = PyArray_GETITEM(axes,PyArray_GETPTR1(axes,0));
	PyObject* axdmt = PyArray_GETITEM(axes,PyArray_GETPTR1(axes,1));
	PyObject* axint = PyArray_GETITEM(axes,PyArray_GETPTR1(axes,2));
	//plot ft
	PyObject *keywords = Py_BuildValue("{s:O}", "aspect",PyString_FromString("auto"));

	callKey(axft,"imshow",PyTuple_Pack(1,ft_list),keywords);
	callKey(axdmt,"imshow",PyTuple_Pack(1,dmt_list),keywords);
	call(axint,"plot",PyTuple_Pack(1,intensityList));

	//make axis tick labels
	double df = (sd.f1 - sd.f0)/((double) sd.nf);
	double dt = (sd.t_frame)/((double) sd.nt);
	PyArrayObject* ft_freqlocs = (PyArrayObject*) call(axft,"get_yticks",NOARGS);
	PyArrayObject* ft_timelocs = (PyArrayObject*) call(axft,"get_xticks",NOARGS);
	PyObject* ft_freqlabels = getNumpyTickLabels(ft_freqlocs,sd.nf,2,sd.f0,df);
	PyObject* ft_timelabels = getNumpyTickLabels(ft_timelocs,ntf,2,0.0f,dt);

	PyArrayObject* dmt_dmlocs = (PyArrayObject*) call(axdmt,"get_yticks",NOARGS);
	PyArrayObject* dmt_timelocs = (PyArrayObject*) call(axdmt,"get_xticks",NOARGS);
	PyObject* dmt_dmlabels = getNumpyTickLabels(dmt_dmlocs,sd.ndm,2,0.0f,ddm);
	PyObject* dmt_timelabels = getNumpyTickLabels(dmt_timelocs,sd.nt,2,0.0f,dt);

	PyArrayObject* int_freqlocs = (PyArrayObject*) call(axint,"get_xticks",NOARGS);
	PyObject* int_freqlabels = getNumpyTickLabels(int_freqlocs,sd.nf,2,sd.f0,df);

	call(axft,"set_yticklabels",PyTuple_Pack(1,ft_freqlabels));
	call(axft,"set_xticklabels",PyTuple_Pack(1,ft_timelabels));
	call(axft,"set_xlabel",PyTuple_Pack(1,PyString_FromString("Time (s)")));
	call(axft,"set_ylabel",PyTuple_Pack(1,PyString_FromString("Freq (MHz)")));

	call(axdmt,"set_yticklabels",PyTuple_Pack(1,dmt_dmlabels));
	call(axdmt,"set_xticklabels",PyTuple_Pack(1,dmt_timelabels));
	call(axdmt,"set_xlabel",PyTuple_Pack(1,PyString_FromString("Time (s)")));
	call(axdmt,"set_ylabel",PyTuple_Pack(1,PyString_FromString("DM")));

	call(axint,"set_xticklabels",PyTuple_Pack(1,int_freqlabels));
	call(axint,"set_xlabel",PyTuple_Pack(1,PyString_FromString("Freq (MHz)")));
	//call(mpl,"title",PyTuple_Pack(1,PyString_FromString("FT/DMT/INT-INTENSITY")));
	call(mpl, "savefig", PyTuple_Pack(1,PyString_FromString(td.fname.c_str())));
	call(mpl, "close", NOARGS);
}

void multipanel_transient_save_simple(double *ft, double *dmt, double f0, double f1, int nf, double t_chunk, int nt, double maxdm, int ndm, trigger_data td)
{
	search_data sd;
	sd.nt = nt;
	sd.t_frame = t_chunk;
	sd.dt = (double) (t_chunk/((double) nt));
	sd.ndm = ndm;
	sd.maxdm = maxdm;
	sd.ddm = (double) (maxdm/((double) ndm));
	sd.nf = nf;
	sd.f0 = f0;
	sd.f1 = f1;
	sd.df = (double) ((f1 - f0)/((double) nf));

	if(Py_IsInitialized() == 0){
		Py_Initialize();
	}
	PyObject* mplString = PyString_FromString((char*)"matplotlib.pyplot");
	PyObject* mpl = PyImport_Import(mplString);
	multipanel_transient_save(ft,dmt,nt,mpl,td,sd);
	//Py_Finalize();
}
// void test_imshow(PyObject *mpl)
// {
// 	double *dat = test_data(100,100);
// 	imshow_save(dat,100,100,mpl,"mytest.png");
// }

float* d_to_f(double* in, int len){
	float* ret = new float[len];
	for(int i = 0; i < len; i++){
		ret[i] = (float) in[i];
	}
	return ret;
}

void test_multipanel(PyObject *mpl)
{
	int width = 10;
	search_data sd = {800.0,400.0,1024,10.0,1024,200.0,1024};
	trigger_data td;
	td.fname = std::string("test_trigger.png");
	td.t_width = width;
	td.t_ind = 512;
	td.dm_ind = 1024;
	td.snr = 8.0f;
	td.nt_trig = 1024;
	td.ndm_trig = 1024;
	td.nf_trig = 1024;
	double *ft = randArray(0.0,1.0,1024 * 1024);
	double *dmt = zeros_d(1024 * 1024);
	double dt = sd.t_frame/((double) sd.nt);
	double df = (sd.f1 - sd.f0)/((double) sd.nf);
	double freq;
	int delay;
	for(int i = 0; i < 1024; i++){
		freq = sd.f0 + i * df;
		delay = (int) (disp_delta(freq,sd.maxdm,sd)/dt);
		//printf("delay %i\n",delay);
		for(int j = 0; j < width; j++){
			ft[i * 1024 + 512 + delay + j] += 1.0;
		}
	}

	td.dm_data = d_to_f(dmt,1024 * 1024);
	td.ft_data = d_to_f(ft,1024 * 1024);
	multipanel_transient_save(ft,dmt,1024,mpl,td,sd);
}

int main()
{
	Py_Initialize();
	import_array();
	PyObject* mplString = PyString_FromString((char*)"matplotlib.pyplot");
	PyObject* mplModule = PyImport_Import(mplString);
	//test_imshow(mplModule);
	test_multipanel(mplModule);
	Py_Exit(0);
}

// int main()
// {
// 	Py_Initialize();
// 	PyObject* mplString = PyString_FromString((char*)"matplotlib.pyplot");
// 	PyObject* mplModule = PyImport_Import(mplString);
// 	PyObject* plotFun = PyObject_GetAttrString(mplModule,(char*)"plot");
// 	PyObject* a = PyList_New(5);
// 	for(int i = 0; i < 5; i++){
// 		PyList_SetItem(a,i,PyFloat_FromDouble((double) i));
// 	}
// 	PyObject* args = PyTuple_Pack(1,a);
// 	PyObject* noargs = PyTuple_Pack(0); 
// 	PyObject* result = PyObject_CallObject(plotFun,args);
// 	PyObject* showFun = PyObject_GetAttrString(mplModule,(char*)"show");
// 	result = PyObject_CallObject(showFun,noargs);
// 	Py_Exit(0);
// }