source('function_hap.R')
setwd(dir='input_for_R/')


list_snp_fis=c(2,4,6,6,8,9,16,21,5)
corresponding_amp=c('GSMUA_Achr6G28840_SWA2_548195_amp1.csv','GSMUA_Achr7G04520_SWA2_458345_amp1.csv','GSMUA_Achr7G04520_SWA2_458345_amp1.csv','GSMUA_Achr7G07160_SWA2_45268_amp1.csv','GSMUA_Achr7G07160_SWA2_45268_amp1.csv','GSMUA_Achr7G07160_SWA2_45268_amp1.csv','GSMUA_Achr7G07160_SWA2_45268_amp1.csv','GSMUA_Achr7G21610_SWA2_21429_amp1.csv','GSMUA_Achr8G28550_SWA2_21616_amp1.csv')


t=fast_analysis('GSMUA_Achr6G28840_SWA2_548195_amp1_ind.txt','GSMUA_Achr6G28840_SWA2_548195_amp1.csv') 
# to view all haplotype : plot_group(t,rownames(t),dim(t)[1],rep(1,dim(t)[1]))
# another example: make_tree(t,'',F,read.table('GSMUA_Achr6G28840_SWA2_548195_amp1_ind.txt',sep='\t',header=F,colClasses="character"),F,'GSMUA_Achr6G28840_SWA2_548195_amp1.csv')

d=fast_analysis('GSMUA_Achr7G04520_SWA2_458345_amp1_ind.txt','GSMUA_Achr7G04520_SWA2_458345_amp1.csv')

e=fast_analysis('GSMUA_Achr7G04520_SWA2_458345_amp2_ind.txt','GSMUA_Achr7G04520_SWA2_458345_amp2.csv')

f=fast_analysis('GSMUA_Achr7G07160_SWA2_45268_amp1_ind.txt','GSMUA_Achr7G07160_SWA2_45268_amp1.csv')

g=fast_analysis('GSMUA_Achr7G21610_SWA2_21429_amp1_ind.txt','GSMUA_Achr7G21610_SWA2_21429_amp1.csv')
h=fast_analysis('GSMUA_Achr7G21610_SWA2_21429_amp2_ind.txt','GSMUA_Achr7G21610_SWA2_21429_amp2.csv')
j=fast_analysis('GSMUA_Achr8G28550_SWA2_21616_amp1_ind.txt','GSMUA_Achr8G28550_SWA2_21616_amp1.csv')
k=fast_analysis('GSMUA_Achr8G28550_SWA2_21616_amp2_ind.txt','GSMUA_Achr8G28550_SWA2_21616_amp2.csv')



# test=unlist(strsplit(test_hap,split=''))
# test[5]=paste(c('*',test[5],'*'),collapse='')
# paste(test,collapse='')

export_haplotype(t,'GSMUA_Achr6G28840_SWA2_548195_amp1_uniq_hap.fasta')
export_haplotype(d,'GSMUA_Achr7G04520_SWA2_458345_amp1_uniq_hap.fasta')
export_haplotype(e,'GSMUA_Achr7G04520_SWA2_458345_amp2_uniq_hap.fasta')
export_haplotype(f,'GSMUA_Achr7G07160_SWA2_45268_amp1_uniq_hap.fasta')
export_haplotype(g,'GSMUA_Achr7G21610_SWA2_21429_amp1_uniq_hap.fasta')
export_haplotype(h,'GSMUA_Achr7G21610_SWA2_21429_amp2_uniq_hap.fasta')
export_haplotype(j,'GSMUA_Achr8G28550_SWA2_21616_amp1_uniq_hap.fasta')
export_haplotype(j,'GSMUA_Achr8G28550_SWA2_21616_amp2_uniq_hap.fasta')