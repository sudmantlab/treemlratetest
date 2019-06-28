
all: rates.html

rates.html: divergences.txt dup_del_rates_by_lineage_sum_table.txt 
	python assess_dup_del_rate_significance.py --fn_rates dup_del_rates_by_lineage_sum_table.txt --fn_divergence divergences.txt >$@
		
