name_gene=['GSMUA_Achr7G21610_SWA2_21429','GSMUA_Achr7G07160_SWA2_45268','GSMUA_Achr8G28550_SWA2_21616','GSMUA_Achr6G28840_SWA2_548195','GSMUA_Achr7G04520_SWA2_458345']


for name in name_gene:
  uu=open(name+'.reads','r')


  u=uu.readlines()
  uu.close()
  counter=1
  indice=1
  v=open(name+'_'+str(indice)+'_reads.fasta','w')
  for i in range(len(u)):
    if (counter%20000==0):
      v.close()
      indice+=1
      v=open(name+'_'+str(indice)+'_reads.fasta','w')
    u[i]='>'+str(counter)+'\t'+u[i]
    u[i]=u[i].replace('\t','_',3)
    u[i]=u[i].replace('\t','\n',1)
    v.write(u[i])
    counter+=1
  
  v.close()