/* 
   Copyright (c) 2008 - Chris Buckley. 

   Permission is granted for use and modification of this file for
   research, non-commercial purposes. 
*/
#include "common.h"
#include "sysfunc.h"
#include "trec_eval.h"
#include "functions.h"
#include "trec_format.h"

static int 
te_calc_err (const EPI *epi, const REL_INFO *rel_info,
		  const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval);

/* See trec_eval.h for definition of TREC_MEAS */
TREC_MEAS te_meas_err =
    {"err",
     "  Expected Reciprocal Rank at cutoffs.\n\
        Olivier Chapelle, Donald Metzler, Ya Zhang, and Pierre Grinspan. \n\
        Expected reciprocal rank for graded relevance. \n\
        CIKM 2009, 621--630, New York, NY, USA, 2009. \n\
        \n\
        default p(stop|rel) set to 0.90\n",
     te_init_meas_s_float,
     te_calc_err,
     te_acc_meas_s,
     te_calc_avg_meas_s,
     te_print_single_meas_s_float,
     te_print_final_meas_s_float,
     NULL, -1};

static int 
te_calc_err (const EPI *epi, const REL_INFO *rel_info,
		  const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval)
{
    RES_RELS res_rels;
    double sum,p;
    long i;
   
    if (UNDEF == te_form_res_rels (epi, rel_info, results, &res_rels))
	return (UNDEF);

    sum = 0.0;
    p = 1.0;
    for (i = 0; i < res_rels.num_ret; i++) {
    	if (res_rels.results_rel_list[i] >= epi->relevance_level) {
            sum += (epi->err_R / (i+1.0) ) * p;
            p = p * (1.0 - epi->err_R);
    	}
    }
	eval->values[tm->eval_index].value = sum;

    return (1);
}
