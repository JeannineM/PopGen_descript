def recreate_hap(name_gene):
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
	  
    for j in range(2):
      amp[l][j]=int(amp[l][j])-1
      for k in mod:
	if k<=amp[l][j]:
	  amp[l][j]+=1
    
    for j in range(len(list_pos)):
      if seq[list_pos[j]]!=ref_SNP[j]:
	print(''.join(str(list_pos_ref)))
	print(''.join(str(mod)))
	print(''.join(str(list_pos)))
	print(name_gene+' amp '+str(l-1))
	print('error '+str(j)+' '+ str(list_pos[j])+' '+seq[list_pos[j]]+' '+ref_SNP[j])
	
    z=open('realign_uniq_hap/'+name_gene+'_amp'+str(l+1)+'_uniq_hap.fasta','r')
    wz=file('realign_uniq_hap/'+name_gene+'_amp'+str(l+1)+'_uniq_hap_amplicon.fasta','w')
    for line in z:
      if line[0]=='>':
	wz.write(line)
      else:
	hap=line.rstrip('\n')
	new_seq=seq[:]
	sub_list_snp=list_pos[snp_per_amp[l]:snp_per_amp[l+1]]
	for j in range(len(sub_list_snp)):
	  new_seq=new_seq[:sub_list_snp[j]]+hap[j]+new_seq[sub_list_snp[j]+1:]
	wz.write(new_seq[amp[l][0]:(amp[l][1]+1)]+'\n')

    z.close()
    wz.close()
  return;

vect_name_gene=['GSMUA_Achr7G07160_SWA2_45268','GSMUA_Achr7G21610_SWA2_21429','GSMUA_Achr8G28550_SWA2_21616','GSMUA_Achr6G28840_SWA2_548195','GSMUA_Achr7G04520_SWA2_458345']
#name_gene=vect_name_gene[1]
for name in vect_name_gene:
  recreate_hap(name)