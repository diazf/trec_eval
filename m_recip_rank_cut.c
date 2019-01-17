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
te_calc_recip_rank_cut (const EPI *epi, const REL_INFO *rel_info,
		    const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval);
static long long_cutoff_array[] = {5, 10, 25, 50, 100};
static PARAMS default_recip_rank_cutoffs = {
NULL, sizeof (long_cutoff_array) / sizeof (long_cutoff_array[0]),
&long_cutoff_array[0]};

/* See trec_eval.h for definition of TREC_MEAS */
TREC_MEAS te_meas_recip_rank_cut =
    {"recip_rank_cut",
    "    Reciprocal Rank of the first relevant retrieved doc.\n\
    Measure is most useful for tasks in which there is only one relevant\n\
    doc, or the user only wants one relevant doc.\n\
    Cutoffs must be positive without duplicates\n\
    Default params: -m recip_rank_cut.5,10,25,50,100\n\
    Based on the original recip_rank implementation. \n",
     te_init_meas_a_float_cut_long, 
     te_calc_recip_rank_cut,
     te_acc_meas_a_cut,
     te_calc_avg_meas_a_cut,
     te_print_single_meas_a_cut,
     te_print_final_meas_a_cut,
     (void *) &default_recip_rank_cutoffs, -1};

static int 
te_calc_recip_rank_cut (const EPI *epi, const REL_INFO *rel_info,
                        const RESULTS *results, const TREC_MEAS *tm,
                        TREC_EVAL *eval)
{
    long  *cutoffs = (long *) tm->meas_params->param_values;
    RES_RELS res_rels;
    long i;

    if (UNDEF == te_form_res_rels (epi, rel_info, results, &res_rels))
	return (UNDEF);
    
    for (i = 0; i < res_rels.num_ret; i++) {
    	if (res_rels.results_rel_list[i] >= epi->relevance_level)
    	    break;
    }
    if (i<res_rels.num_ret){
        for (int cutoff_index = 0; cutoff_index < tm->meas_params->num_params ; cutoff_index++){
            if (i < cutoffs[cutoff_index]){
                eval->values[tm->eval_index + cutoff_index].value = (double) 1.0 / (double) (i+1);
            }
        }
    }

    return (1);
}
