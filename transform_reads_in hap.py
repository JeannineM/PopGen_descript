import numpy as np

def extract_hap(seq, output,vect_ind,vect_hap,list_pos,ind):
  new_seq=[]
  for j in list_pos:
    new_seq.append(seq[j])
  hap=''.join(new_seq)

  
  if hap in vect_hap:
    ind_hap=vect_hap.index(hap)
  else:
    ind_hap=len(vect_hap)
    vect_hap.append(hap)
    output=np.vstack((output,np.zeros((1, 313))))
    
  if ind in vect_ind:
    ind_ind=vect_ind.index(ind)
  else:
    ind_ind=len(vect_ind)
    vect_ind.append(ind)
    
  output[ind_hap,ind_ind]+=1
  return output;

def count_hap(name_gene):
  v=open(name_gene+'.SNPs','r')

  w=open('amplicon_pos.txt','r')
  w=w.readlines()
  for i in range(len(w)):
    w[i]=w[i].replace("|","_")
    w[i]=w[i].rsplit('\t')

  amp=[]
  for i in range(len(w)):
    if w[i][0]==name_gene:
      amp.append(w[i][1:3])
    
  list_pos_ref=[]
  ref_SNP=[]
  for k in v:
    list_pos_ref.append(int(k.rsplit('\t')[1]))
    ref_SNP.append(k.rsplit('\t')[2])

  counter=0
  snp_per_amp=[0]
  for o in range(len(amp)):
    for j in list_pos_ref:
      if j in range(int(amp[o][0]),int(amp[o][1])):
	counter+=1
    snp_per_amp.append(counter)
  
  for l in range(len(snp_per_amp)-1):
    u=open('align_seq/'+name_gene+'_amp'+str(l+1)+'_Q10_reads.aln','r')
    bool_1=False
    seq=''
    list_pos=list_pos_ref[:]
    for line in u:
      if(line[0:6]=='>0_REF' or bool_1):
	if(line[0]!='>' and bool_1):
	  bool_2=True
	  seq=seq+line.strip('\n')
	bool_1=True
	if (line[0]=='>' and len(seq)>0):
	    bool_1=False

    mod=[]
    while ('-' in seq):
      mod.append(seq.index('-'))
      seq=seq.replace('-','_',1)

    for j in range(len(list_pos)):
      list_pos[j]+=-1
      for k in mod:
	if k<=list_pos[j]:
	  list_pos[j]+=1
    
    for j in range(len(list_pos)):
      if seq[list_pos[j]]!=ref_SNP[j]:
	print(''.join(str(list_pos_ref)))
	print(''.join(str(mod)))
	print(''.join(str(list_pos)))
	print(name_gene+' amp '+str(l-1))
	print('error '+str(j)+' '+ str(list_pos[j])+' '+seq[list_pos[j]]+' '+ref_SNP[j])

  
    u=open('align_seq/'+name_gene+'_amp'+str(l+1)+'_Q10_reads.aln','r')

    output=np.zeros((1, 313))
    vect_hap=[]
    vect_ind=[]

    i=u.readline()
    info=i.replace('_','\t',2).rsplit('\t')
    for i in u:
      if i[0]=='>':
	output=extract_hap(seq, output,vect_ind,vect_hap,list_pos[snp_per_amp[l]:snp_per_amp[l+1]],info[1])
	seq=''
	info=i.replace('_','\t',2).rsplit('\t')
      else:
	seq=seq+i.strip('\n')
    output=extract_hap(seq, output,vect_ind,vect_hap,list_pos[snp_per_amp[l]:snp_per_amp[l+1]],info[1])
    

    np.savetxt(name_gene+'_amp'+str(l+1)+'.csv',np.transpose(output[:(len(output)-1)]),delimiter='\t',header='\t'.join(vect_hap),comments='')
      
    wz=file(name_gene+'_amp'+str(l+1)+'_ind.txt','w')
    wz.write('\t'.join(vect_ind))
    wz.close()
  return;

#vect_name_gene=['GSMUA_Achr7G07160_SWA2_45268','GSMUA_Achr7G21610_SWA2_21429','GSMUA_Achr8G28550_SWA2_21616','GSMUA_Achr6G28840_SWA2_548195','GSMUA_Achr7G04520_SWA2_458345']
#name_gene=vect_name_gene[1]
#for name in vect_name_gene:
#  count_hap(name)

