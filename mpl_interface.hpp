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
	double dt;
	double df;
	double ddm;
} search_data;

class trigger_data;

// struct trigger_data{
// 	std::string fname;
// 	int t_width;
// 	//by convention take the start of pulse; pulse
// 	//lies in [t_ind, t_ind + t_width)
// 	int t_ind;
// 	int dm_ind;
// 	double snr;
// 	int nt_trig;
// 	int ndm_trig;
// 	int nf_trig;
// 	float* dm_data;
// 	float* ft_data;
// 	int* dm_dim;
// 	int* ft_dim;    

// 	bool operator < (const trigger_data& b) const
//     {
//         return (snr < b.snr);
//     }
// };
// typedef struct trigger_data trigger_data;

// class trigger_data{
// 	public:
// 		std::string fname;
// 		int t_width;
// 		//by convention take the start of pulse; pulse
// 		//lies in [t_ind, t_ind + t_width)
// 		int t_ind;
// 		int dm_ind;
// 		double snr;
// 		int nt_trig;
// 		int ndm_trig;
// 		int nf_trig;
// 		float* dm_data;
// 		float* ft_data;
// 		int* dm_dim;
// 		int* ft_dim;    

// 		bool operator < (const trigger_data& b) const
// 	    {
// 	        return (snr < b.snr);
// 	    }
// };
// //typedef struct trigger_data trigger_data;


PyObject *callKey(PyObject* target, std::string fname, PyObject *args, PyObject *keywords);
PyObject *call(PyObject* target, std::string fname, PyObject *args);

void imshow_save(double *ar, int nrow, int ncol, PyObject* mpl, std::string fname,bool cbar);
void imshow_save_simple(double *ar, int nrow, int ncol, std::string fname);
void imshow_show(double *ar, int nrow, int ncol, PyObject* mpl);
void multipanel_transient_save(double *ft, double *dmt, int ntdm, PyObject *mpl, trigger_data td, search_data sd);
void multipanel_transient_save_simple(double *ft, double *dmt, double f0, double f1, int nf, double t_chunk, int nt, double maxdm, int ndm, trigger_data td);
void plot_save(double* data, int len, std::string fname, PyObject *mpl);
void plot_save_simple(double* data, int len, std::string fname);