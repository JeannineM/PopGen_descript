import random

for i in range(16):
  n=random.randint(1,2147483647)
  u=open('bash_script_n'+str(i+1)+'.sh','w')
  text="R CMD BATCH --no-save --no-restore '--args seed="+str(n)+ " output=\"model1_seed_"+str(n)+"_n"+str(i+1)+".txt\"' run_FD_4pop_cluster.R Routput_"+str(i+1)+".Rout\n"
  u.write(text)
  u.close()
  
