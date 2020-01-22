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
double log2(double x);

static int 
te_calc_gerr_cut (const EPI *epi, const REL_INFO *rel_info,
		  const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval);
static long long_cutoff_array[] = {1, 5, 10, 15, 20, 30, 100, 200, 500, 1000};
static PARAMS default_gerr_cutoffs = {
    NULL, sizeof (long_cutoff_array) / sizeof (long_cutoff_array[0]),
    &long_cutoff_array[0]};

/* See trec_eval.h for definition of TREC_MEAS */
TREC_MEAS te_meas_gerr_cut =
    {"gerr_cut",
     "  Generalized Expected Reciprocal Rank at cutoffs .\n\
        see section 7.2, \n\
        Olivier Chapelle, Donald Metzler, Ya Zhang, and Pierre Grinspan. \n\
        Expected reciprocal rank for graded relevance. \n\
        CIKM 2009, 621--630, New York, NY, USA, 2009. \n\
        \n\
        default p(stop|rel) set to 0.90\n\
        default gamma set to 0.80\n\
        Cutoffs must be positive without duplicates\n\
        Default params: -m gerr_cut.1,5,10,15,20,30,100,200,500,1000\n",
     te_init_meas_a_float_cut_long,
     te_calc_gerr_cut,
     te_acc_meas_a_cut,
     te_calc_avg_meas_a_cut,
     te_print_single_meas_a_cut,
     te_print_final_meas_a_cut,
     (void *) &default_gerr_cutoffs, -1};

static int 
te_calc_gerr_cut (const EPI *epi, const REL_INFO *rel_info,
		  const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval)
{
    long  *cutoffs = (long *) tm->meas_params->param_values;
    long cutoff_index = 0;
    RES_RELS res_rels;
    double sum,p;
    long i;
   
    if (UNDEF == te_form_res_rels (epi, rel_info, results, &res_rels))
	return (UNDEF);

    sum = 0.0;
    p = 1.0;
    for (i = 0; i < res_rels.num_ret; i++) {
        if (i == cutoffs[cutoff_index]){
            int effectiveMaxRel = res_rels.num_rel;
            if (i < res_rels.num_rel){
                effectiveMaxRel = i;
            }
            double norm = epi->err_R * ( (1.0 - pow((epi->rbp_p * (1.0 - epi->err_R)), effectiveMaxRel)) / (1.0 - (epi->rbp_p * (1.0 - epi->err_R))) );
            eval->values[tm->eval_index + cutoff_index].value = sum / norm;
            if (++cutoff_index == tm->meas_params->num_params)
                break;
        }
    	if (res_rels.results_rel_list[i] >= epi->relevance_level) {
            sum += epi->err_R * p;
            p = p * epi->rbp_p * (1.0 - epi->err_R);
    	}else{
    	    p = p * epi->rbp_p;
    	}
    }
    /* calculate values for those cutoffs not achieved */
    while (cutoff_index < tm->meas_params->num_params) {
        int i = cutoffs[cutoff_index];
        int effectiveMaxRel = res_rels.num_rel;
        if (i < res_rels.num_rel){
            effectiveMaxRel = i;
        }
        double norm = epi->err_R * ( (1.0 - pow((epi->rbp_p * (1.0 - epi->err_R)), effectiveMaxRel)) / (1.0 - (epi->rbp_p * (1.0 - epi->err_R))) );
        
        eval->values[tm->eval_index + cutoff_index].value = sum / norm;
        cutoff_index++;
    }

    return (1);
}
