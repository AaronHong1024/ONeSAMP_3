def stat1(data):
  print(data[0][1][2:])
  
  alleleA=[]
  alleleB=[]
  
  print(len(data))
  print(len(data[0]))
  
  for i in range(len(data)):
      temA=[]
      temB=[]
      for j in range(1,len(data[i])):
          temA.append(data[i][j][:2])
          temB.append(data[i][j][2:])
      
      alleleA.append(temA)
      alleleB.append(temB)
  
  allcnt=[]
  for j in range(len(alleleA[0])):
      newdic={}
      for i in range(len(alleleA)):
          if(alleleA[i][j] in newdic):
              newdic[alleleA[i][j]]+=1
          else:
              newdic[alleleA[i][j]]=1
          
          if(alleleB[i][j] in newdic):
              newdic[alleleB[i][j]]+=1
          else:
              newdic[alleleB[i][j]]=1
      
      allcnt.append(newdic)
      
  print(allcnt[0])
  di=[]
  totspots=len(data)*2
  print(totspots)
  for i in range(len(allcnt)):
      vals=[]
      for key,value in allcnt[i].items():
          vals.append(value/totspots)
      
      if(len(vals)>1):
          di.append(vals[0]*vals[0]+vals[1]*vals[1]+2*vals[0]*vals[1]-1)
      else:
          di.append(0)
      
  hits=0
  
  for i in range(len(alleleA)):
      for j in range(len(alleleB)):
          for k in range(len(alleleA[0])):
              if(alleleA[i][j]==alleleB[i][j]):
                  hits+=1
              
  
  print(hits)
  
  