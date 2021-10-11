import math

def stat1(data):
  
  alleleA=[] # same dimensions as data but holds the first allele in a locus
  alleleB=[] # same dimensions as data but holds the second allele in a locus
  
  for i in range(len(data)): # fills up alleleA and alleleB
      temA=[]
      temB=[]
      for j in range(1,len(data[i])):
          temA.append(data[i][j][:2])
          temB.append(data[i][j][2:])
      
      alleleA.append(temA)
      alleleB.append(temB)
  
  allcnt=[] #a 1D array, each element is a dictionary of alleles with frequency counts per loci
  homoloci=[] # maintains frequency counts of homologous alleles per loci
  for j in range(len(alleleA[0])): # fills up allcnt
      newdic={}
      temp=0
      for i in range(len(alleleA)):
          if(alleleA[i][j]==alleleB[i][j]):
              temp+=1

          if(alleleA[i][j] in newdic):
              newdic[alleleA[i][j]]+=1
          else:
              newdic[alleleA[i][j]]=1
          
          if(alleleB[i][j] in newdic):
              newdic[alleleB[i][j]]+=1
          else:
              newdic[alleleB[i][j]]=1
      
      homoloci.append(temp/len(data)) 
      allcnt.append(newdic)
      
  print(allcnt[0])
  di=[] # a 1D array that holds the departures of each loci from Hardy-Weinberg equilibrium
  totspots=len(data)*2 #total number of alleles per locus. used in computing allele freq per locus

  for i in range(len(allcnt)): #fills up di 
      vals=[]
      for key,value in allcnt[i].items():
          vals.append(value/totspots)
      di.append(homoloci[i]*homoloci[i] - vals[0]*vals[0])
      
  hits=0
  
  running_sum=0
  
  for i in range(len(alleleA)):
      for j in range(len(alleleB)):
          keysi=[]
          keysj=[]

          for key,value in allcnt[i].items():
              keysi.append(key)
              
          for key,value in allcnt[j].items():
              keysj.append(key)
          
          if(len(keysj)<2):
              continue
          
          alA=keysi[0]
          alB=keysj[1]
          
          hits=0
          for k in range(len(alleleA[0])):
              if( (alleleA[i][k]==keysi or alleleB[i][k]==keysi) and (alleleA[j][k]==keysj or alleleB[j][k]==keysj) ):
                  hits+=1
          ai=allcnt[i][alA]/totspots
          bj=allcnt[j][alB]/totspots
          x=(hits/len(data) - ai*bj)/((ai*(1-ai)+di[i])*(bj*(1-bj)+di[j]))
          running_sum+=math.sqrt(abs(x))
  
  numloci=len(alleleA[0])
  print(2*running_sum/(numloci*(numloci-1)))
  #return 2*running_sum/(numloci*(numloci-1)) 
  
  