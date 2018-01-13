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
te_calc_num_nonjudged_ret(const EPI *epi, const REL_INFO *rel_info,
			       const RESULTS *results,  const TREC_MEAS *tm,
			       TREC_EVAL *eval);

/* See trec_eval.h for definition of TREC_MEAS */
TREC_MEAS te_meas_num_nonjudged_ret =
    {"num_nonjudged_ret",
     "    Number of non-judged documents retrieved for topic. \n\
    Not an evaluation number per se, but gives details of retrieval results.\n\
    Summary figure is sum of individual topics, not average.\n",
     te_init_meas_s_long,
     te_calc_num_nonjudged_ret,
     te_acc_meas_s, 
     te_calc_avg_meas_empty,
     te_print_single_meas_s_long,
     te_print_final_meas_s_long,
     NULL, -1};

static int 
te_calc_num_nonjudged_ret (const EPI *epi, const REL_INFO *rel_info,
			       const RESULTS *results,  const TREC_MEAS *tm,
			       TREC_EVAL *eval)
{
    RES_RELS res_rels;

    if (UNDEF == te_form_res_rels (epi, rel_info, results, &res_rels))
	return (UNDEF);

    eval->values[tm->eval_index].value = (double)
	res_rels.num_nonpool +
	res_rels.num_unjudged_in_pool;

    return (1);
}
