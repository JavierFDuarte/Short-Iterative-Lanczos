import pandas as pd     #read and write data
df1 = pd.read_table("CoefEvTemp.txt",sep=r'\s{2,}')
print(df1)

n = 4 

re_coef=list(df1["recoef"])
im_coef=list(df1["imcoef"])
ci=[ ]
for i in range(0,n):
    ci.append(complex(re_coef[i],im_coef[i]))
#print("re_cn=",re_coef)
#print("im_cn=",im_coef)
print("ci=",ci)
#print(df1)



kk=0
initialstate = { }
for i in range(0,n):
    fname = "temporalfile_v"+str(kk)+".csv"
    df = pd.read_csv(filepath_or_buffer=fname)
    ky = list(df["ky"])
    rpt= list(df["rpt"])
    ipt= list(df["ipt"])
    #print(ky,rpt,ipt)
    v = { }
    s=0
    for j in ky:
        v[j] = v.get(j,0.0)+complex(rpt[s],ipt[s])
        s=s+1
    print(v)
    for j in v.keys():
        v[j]=ci[i]*v[j]
    print(v)
    print("-------------------")
    for k in v.keys():
        initialstate[k]= initialstate.get(k,0.0) + v[k]
    
    kk=kk+1

print("initialstate=",initialstate)
    

kyinitialstate = [ ]
reinitialstate = [ ]
iminitialstate = [ ]
for i in initialstate.keys():
    kyinitialstate.append(i)
    reinitialstate.append(initialstate[i].real)
    iminitialstate.append(initialstate[i].imag)

dic = { "instateky" : kyinitialstate , "instatere": reinitialstate, "instateim": iminitialstate}

print(dic)

data = pd.DataFrame(dic)
print(data)

data.to_csv(path_or_buf="initialstate.csv")