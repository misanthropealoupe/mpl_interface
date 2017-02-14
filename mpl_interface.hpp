#include "Python.h"
#include <string>

typedef struct search_data{
	double f0;
	double f1;
	int nf;
	double t_frame;
	int nt;
	double maxdm;
	int ndm;
} search_data;

typedef struct trigger_data{
	std::string fname;
	int t_width;
	//by convention take the start of pulse; pulse
	//lies in [t_ind, t_ind + t_width)
	int t_ind;
	int dm_ind;
	double snr;
} trigger_data;

void imshow_save(double *ar, int nrow, int ncol, PyObject* mpl, char* fname);
void imshow_show(double *ar, int nrow, int ncol, PyObject* mpl);
void multipanel_transient_save(double *ft, double *dmt, int ntdm, PyObject *mpl, trigger_data td, search_data sd);