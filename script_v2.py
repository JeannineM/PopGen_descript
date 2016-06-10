import numpy as np

def extract_hap(seq, output,vect_ind,vect_hap,list_pos,amp):
  new_seq=[]
  for j in list_pos:
    new_seq.append(seq[j])
  hap=''.join(new_seq)
  
  for l in range(len(amp)):
    local_hap=hap[snp_per_amp[l]:snp_per_amp[l+1]]
    if '_'*len(local_hap)!=local_hap:
    
      if local_hap in vect_hap[l]:
	ind_hap=vect_hap[l].index(local_hap)
      else:
	ind_hap=len(vect_hap[l])
	vect_hap[l].append(local_hap)
	output[l]=np.vstack((output[l],np.zeros((1, 313))))
	
      if info[1] in vect_ind:
	ind_ind=vect_ind.index(info[1])
      else:
	ind_ind=len(vect_ind)
	vect_ind.append(info[1])
      
    output[l][ind_hap,ind_ind]+=1
      
  return;

def count_hap(name_gene):
  u=open('align_seq/'+name_gene+'_Q10_reads.aln','r')

  v=open(name_gene+'.SNPs','r')

  w=open('amplicon_pos.txt','r')
  w=w.readlines()
  for i in range(len(w)):
    w[i]=w[i].replace("|","_")
    w[i]=w[i].rsplit('\t')
    
  gene_name=name_gene


  amp=[]
  for i in range(len(w)):
    if w[i][0]==gene_name:
      amp.append(w[i][1:3])
    
  list_pos=[]
  ref_SNP=[]
  for k in v:
    list_pos.append(int(k.rsplit('\t')[1]))
    ref_SNP.append(k.rsplit('\t')[2])

  output={}
  vect_hap={}
  snp_per_amp=[0]
  counter=0
  for o in range(len(amp)):
    output[o]=np.zeros((1, 313))
    vect_hap[o]=[]
    for j in list_pos:
      if j in range(int(amp[o][0]),int(amp[o][1])):
	counter+=1
    snp_per_amp.append(counter)

  bool_1=False
  seq=''
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
    for k in mod:
      if k<=list_pos[j]:
	list_pos[j]+=1
    list_pos[j]+=-1

  for j in range(len(list_pos)):
    if seq[list_pos[j]]!=ref_SNP[j]:
      print('error')

  u=open('align_seq/'+name_gene+'_Q10_reads.aln','r')


  vect_ind=[]

  i=u.readline()
  info=i.replace('_','\t',2).rsplit('\t')
  for i in u:
    if i[0]=='>':
      extract_hap(seq, output,vect_ind,vect_hap,list_pos,amp)
      seq=''
      info=i.replace('_','\t',2).rsplit('\t')
    else:
      seq=seq+i.strip('\n')
  extract_hap(seq, output,vect_ind,vect_hap,list_pos,amp)
  

  for l in range(len(amp)):
    np.savetxt(gene_name+'_amp_'+str(l+1)+'.csv',np.transpose(output[l][:(len(output[l])-1)]),delimiter='\t',header='\t'.join(vect_hap[l]),comments='')
    
  wz=file(gene_name+'_ind.txt','w')
  wz.write('\t'.join(vect_ind))
  wz.close()

name_gene=['GSMUA_Achr7G07160_SWA2_45268','GSMUA_Achr7G21610_SWA2_21429_450','GSMUA_Achr8G28550_SWA2_21616','GSMUA_Achr6G28840_SWA2_548195','GSMUA_Achr7G04520_SWA2_458345']

count_hap(name_gene[0])