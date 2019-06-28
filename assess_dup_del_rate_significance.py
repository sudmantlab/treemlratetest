import csv
import numpy as np
import scipy.stats as stat
from sys import stderr
import argparse


def html_table(lol):
    ret=""
    #ret.append('data:text/html,'
    ret+= '<table>'
    for sublist in lol:
        ret+= '  <tr><td>'
        ret+= '    </td><td>'.join(sublist)
        ret+= '  </td></tr>'
    ret+='</table>'
    return ret

def test_1_v_2_rates(events_by_lin,branch_lens_by_lin,all_lins,test_lins,div_by=1e6,units=["Mb","My"],description='',lt=False):
    
    sum_events_all=np.sum(np.array([events_by_lin[lin] for lin in all_lins ]))/div_by
    sum_lens_all=np.sum(np.array([branch_lens_by_lin[lin] for lin in all_lins ]))
    lambda_all=sum_events_all/sum_lens_all
    print >>stderr, "\tMODEL 1 %f%s %f%s %f %s/%s"%(sum_events_all,units[0], sum_lens_all,units[1],sum_events_all/sum_lens_all,units[0],units[1])
    bls=branch_lens_by_lin
    bes={k:v/div_by for k,v in events_by_lin.iteritems()}

    L1_terms=[(-lambda_all*bls[lin])+(np.log(lambda_all*bls[lin])*bes[lin])-np.log(bes[lin]) for lin in all_lins]

    L1=np.sum(L1_terms)
    
    df=0
    rest_lins=all_lins
    L2_test_terms=[]
    L2_test_lambdas=[]
    for test_lin in test_lins:
        df+=1
        
        sum_events_test=np.sum(np.array([events_by_lin[lin] for lin in test_lin]))/div_by
        sum_lens_test=np.sum(np.array([branch_lens_by_lin[lin] for lin in test_lin ]))
        lambda_test=sum_events_test/sum_lens_test
        rest_lins=list(set(rest_lins)-set(test_lin))
        L2_test_lambdas.append(lambda_test)
        L2_test_terms+=[(-lambda_test*bls[lin])+(np.log(lambda_test*bls[lin])*bes[lin]) - np.log(bes[lin]) for lin in test_lin]

        #print >>stderr, L2_test_terms,np.product(L2_test_terms), L2_rest_terms,np.product(L2_rest_terms)
        
    sum_events_rest=np.sum(np.array([events_by_lin[lin] for lin in rest_lins  ]))/div_by
    sum_lens_rest=np.sum(np.array([branch_lens_by_lin[lin] for lin in rest_lins  ]))
    if len(rest_lins)==0:
        print >>stderr, "IF YOU ARE TESTING ALL, you need to leave one out!"
        exit(1)
    lambda_rest=sum_events_rest/sum_lens_rest
    L2_rest_terms=[(-lambda_rest*bls[lin])+(np.log(lambda_rest*bls[lin])*bes[lin]) - np.log(bes[lin]) for lin in rest_lins]
    L2=np.sum(L2_test_terms)+np.sum(L2_rest_terms)
    
    print >>stderr, "\tMODEL 2 TEST lambdas %s %s/%s rest %f %s/%s on %d degrees of freedom"%(",".join(["%.2f"%l for l in L2_test_lambdas]),units[0],units[1],lambda_rest,units[0],units[1],df)
    
    v=2*(L2-L1)
    print >>stderr, "\t",v,stat.chisqprob(v, df)
    p=stat.chisqprob(v, df)
    table=[description,'%s/%s'%(units[0],units[1]),'&#955;=%.3f'%lambda_all,'&#955;1 = %.3f<br>%s'%(lambda_rest,"<br>".join(["&#955;%d = %.3f"%(i+2,l) for i,l in enumerate(L2_test_lambdas)])),"%d"%df,"%s%.4g"%(lt and "<" or "",p)]
    return table
if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_rates")
    parser.add_argument("--fn_divergence")
    o = parser.parse_args()

    Fin=open(o.fn_rates)
    Fin_div=open(o.fn_divergence)
    
    branch_lMY_by_lin={}
    branch_l_by_lin={}
    for ldict in csv.DictReader(Fin_div,delimiter='\t'):
        branch_lMY_by_lin[ldict['lineage']]=float(ldict['branch_t'])
        branch_l_by_lin[ldict['lineage']]=float(ldict['branch_d'])*2.867e3 #3 mb genome
    dup_bp_by_lin={}
    del_bp_by_lin={}
    
    del_sites_by_lin={}
    dup_sites_by_lin={}

    del_genes_by_lin={}
    dup_genes_by_lin={}
    
    for ldict in csv.DictReader(Fin,delimiter='\t'):
        #print >>stderr, ldict
        lineage=ldict['lineage']
        dup_bp_by_lin[ldict['lineage']]=float(ldict['copy corrected duplicated basepairs'])
        #del_bp_by_lin[ldict['lineage']]=float(ldict['deleted basepairs'])
        del_bp_by_lin[ldict['lineage']]=float(ldict['deleted bp >5kb'])
        
        dup_sites_by_lin[ldict['lineage']]=float(ldict['duplicated sites'])
        #del_sites_by_lin[ldict['lineage']]=float(ldict['deleted sites'])
        del_sites_by_lin[ldict['lineage']]=float(ldict['deleted sites >5kb'])
        #print >>stderr, ldict['duplicated bp / million years (Mbp)'], "%.2f"%(dup_bp_by_lin[lineage]/branch_lMY_by_lin[lineage]/1e6)
    
    t=[]
    print >>stderr, "TEST THE AFRICAN GREAT APE BRANCH (DUP)"    
    t.append(test_1_v_2_rates(dup_bp_by_lin,branch_lMY_by_lin,['Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Pab', 'Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ptre-Ptrs-Ptrt-Ptrv', 'Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ppa', 'Hde-Hsa', 'Gbeg-Ggod-Ggog','Ggod-Ggog', 'Pab-Ppy', 'Ppy','Gbeg',"Hsa"],[['Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv']],description="acceleration in African great ape against all other branches  "))
    t.append(test_1_v_2_rates(dup_bp_by_lin,branch_l_by_lin,['Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Pab', 'Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ptre-Ptrs-Ptrt-Ptrv', 'Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ppa', 'Hde-Hsa', 'Gbeg-Ggod-Ggog','Ggod-Ggog', 'Pab-Ppy', 'Ppy','Gbeg',"Hsa"],[['Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv']],div_by=1e6,units=['bp','sub']))
    t.append(test_1_v_2_rates(dup_sites_by_lin,branch_lMY_by_lin,['Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Pab', 'Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ptre-Ptrs-Ptrt-Ptrv', 'Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ppa', 'Hde-Hsa', 'Gbeg-Ggod-Ggog','Ggod-Ggog', 'Pab-Ppy', 'Ppy','Gbeg',"Hsa"],[['Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv']],div_by=10,units=['sites','MY'],lt=True))
    
    #print >>stderr, "TEST THE AFRICAN GREAT APE BRANCH AND GORILLA AND the HSA-CHIMP-ANC (DUP)"    
    #t.append(test_1_v_2_rates(dup_bp_by_lin,branch_lMY_by_lin,['Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Pab', 'Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ptre-Ptrs-Ptrt-Ptrv', 'Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ppa', 'Hde-Hsa', 'Gbeg-Ggod-Ggog','Ggod-Ggog', 'Pab-Ppy', 'Ppy','Gbeg',"Hsa"],[['Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv','Gbeg-Ggod-Ggog','Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv']],description="acceleration in African great ape and gorilla and human-chimp linages against all branches "))
    #t.append(test_1_v_2_rates(dup_bp_by_lin,branch_l_by_lin,['Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Pab', 'Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ptre-Ptrs-Ptrt-Ptrv', 'Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ppa', 'Hde-Hsa', 'Gbeg-Ggod-Ggog','Ggod-Ggog', 'Pab-Ppy', 'Ppy','Gbeg',"Hsa"],[['Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv','Gbeg-Ggod-Ggog','Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv']],units=['bp','sub']))
    #t.append(test_1_v_2_rates(dup_sites_by_lin,branch_lMY_by_lin,['Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Pab', 'Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ptre-Ptrs-Ptrt-Ptrv', 'Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ppa', 'Hde-Hsa', 'Gbeg-Ggod-Ggog','Ggod-Ggog', 'Pab-Ppy', 'Ppy','Gbeg',"Hsa"],[['Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv','Gbeg-Ggod-Ggog','Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv']],div_by=10,units=['sites',"MY"],lt=True))
    
    #print >>stderr, "TEST THE AFRICAN GREAT APE BRANCH AND GORILLA AND the HSA-CHIMP-ANC (DUP) ALL DIFFERENT RATES"    
    #t.append(test_1_v_2_rates(dup_bp_by_lin,branch_lMY_by_lin,['Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ptre-Ptrs-Ptrt-Ptrv','Ppa', 'Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Hde-Hsa', 'Gbeg-Ggod-Ggog', 'Pab-Ppy'],[['Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv'],['Gbeg-Ggod-Ggog'],['Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv']]))
    #t.append(test_1_v_2_rates(dup_bp_by_lin,branch_l_by_lin,['Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ptre-Ptrs-Ptrt-Ptrv','Ppa', 'Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Hde-Hsa', 'Gbeg-Ggod-Ggog', 'Pab-Ppy'],[['Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv'],['Gbeg-Ggod-Ggog'],['Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv']],units=['bp','sub']))

    #print >>stderr, "TEST THE AFRICAN GREAT APE BRANCH AND GORILLA AND the HSA-CHIMP-ANC (DUP) 3, same in gorilla and HSA-chimp DIFFERENT RATES"    
    #t.append(test_1_v_2_rates(dup_bp_by_lin,branch_lMY_by_lin,['Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ptre-Ptrs-Ptrt-Ptrv','Ppa', 'Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Hde-Hsa', 'Gbeg-Ggod-Ggog', 'Pab-Ppy'],[['Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv'],['Gbeg-Ggod-Ggog'],['Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv']],description="acceleration in African great ape, chimpanzee-human and gorilla branches each with different rates against the rest of the tree >1.8Mya"))
    #t.append(test_1_v_2_rates(dup_bp_by_lin,branch_l_by_lin,['Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ptre-Ptrs-Ptrt-Ptrv','Ppa', 'Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Hde-Hsa', 'Gbeg-Ggod-Ggog', 'Pab-Ppy'],[['Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv'],['Gbeg-Ggod-Ggog'],['Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv']],units=['bp','sub']))
    
    #print >>stderr, "ALL different rates"
    #t.append(test_1_v_2_rates(dup_bp_by_lin,branch_lMY_by_lin,['Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Pab', 'Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ptre-Ptrs-Ptrt-Ptrv', 'Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ppa', 'Hde-Hsa', 'Gbeg-Ggod-Ggog','Ggod-Ggog', 'Pab-Ppy', 'Ppy','Gbeg',"Hsa"],[['Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv'],['Gbeg-Ggod-Ggog'],['Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv'],['Ppa-Ptre-Ptrs-Ptrt-Ptrv'],['Ptre-Ptrs-Ptrt-Ptrv'],['Ppa'],['Hde-Hsa'],["Hsa"],["Gbeg"],["Pab-Ppy"],["Ppy"],["Ggod-Ggog"]],description="different rates among all branches "))
    #t.append(test_1_v_2_rates(dup_sites_by_lin,branch_lMY_by_lin,['Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Pab', 'Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ptre-Ptrs-Ptrt-Ptrv', 'Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ppa', 'Hde-Hsa', 'Gbeg-Ggod-Ggog','Ggod-Ggog', 'Pab-Ppy', 'Ppy','Gbeg',"Hsa"],[['Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv'],['Gbeg-Ggod-Ggog'],['Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv'],['Ppa-Ptre-Ptrs-Ptrt-Ptrv'],['Ptre-Ptrs-Ptrt-Ptrv'],['Ppa'],['Hde-Hsa'],["Hsa"],["Gbeg"],["Pab-Ppy"],["Ppy"],["Ggod-Ggog"]],units=['sites','MY'],div_by=100,lt=True))

    #t.append(test_1_v_2_rates(dup_bp_by_lin,branch_l_by_lin,['Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Pab', 'Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ptre-Ptrs-Ptrt-Ptrv', 'Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ppa', 'Hde-Hsa', 'Gbeg-Ggod-Ggog','Ggod-Ggog', 'Pab-Ppy', 'Ppy','Gbeg',"Hsa"], [['Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv'],['Gbeg-Ggod-Ggog'],['Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv'],['Ppa-Ptre-Ptrs-Ptrt-Ptrv'],['Ptre-Ptrs-Ptrt-Ptrv'],['Ppa'],['Hde-Hsa'],["Hsa"],["Gbeg"],["Pab-Ppy"],["Ppy"],["Ggod-Ggog"]],units=['bp','sub']))
    
    #print >>stderr, "TEST WITH ALL BRANCHES"
    #t.append(test_1_v_2_rates(dup_bp_by_lin,branch_lMY_by_lin,['Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Pab', 'Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ptre-Ptrs-Ptrt-Ptrv', 'Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ppa', 'Hde-Hsa', 'Gbeg-Ggod-Ggog','Ggod-Ggog', 'Pab-Ppy', 'Ppy','Gbeg',"Hsa"],[['Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv','Gbeg-Ggod-Ggog','Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv']]))
    #t.append(test_1_v_2_rates(dup_bp_by_lin,branch_l_by_lin,['Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Pab', 'Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ptre-Ptrs-Ptrt-Ptrv', 'Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ppa', 'Hde-Hsa', 'Gbeg-Ggod-Ggog','Ggod-Ggog', 'Pab-Ppy', 'Ppy','Gbeg','Hsa'],[['Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv','Gbeg-Ggod-Ggog','Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv']],units=['bp','sub']))
    #t.append(test_1_v_2_rates(dup_bp_by_lin,branch_l_by_lin,['Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Ptre-Ptrs-Ptrt-Ptrv','Ppa', 'Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv', 'Hde-Hsa', 'Gbeg-Ggod-Ggog', 'Pab-Ppy'],[['Gbeg-Ggod-Ggog-Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv'],['Gbeg-Ggod-Ggog'],['Hde-Hsa-Ppa-Ptre-Ptrs-Ptrt-Ptrv'],['Ppa-Ptre-Ptrs-Ptrt-Ptrv'],['Ptre-Ptrs-Ptrt-Ptrv'],['Ppa'],['Hde-Hsa']],units=['bp','sub']))

    t.insert(0,["description","units","Model 1","Model 2",'degrees of freedom','p-value'])
    print  "data:text/html,",html_table(t)
    print >>stderr, "-------------------"    

