import pandas as pd
vector_array = [{0: (1+0j)}, {"Nulo": (0+0j)}, {"Nulo": (0+0j)}, {"Nulo": (0+0j)}]
kk=0
for i in vector_array:
    filename = "temporalfile_v"+str(kk)+".csv"
    key = [ ] 
    realpart = [ ]
    imagpart = [ ]
#    print(filename)
    for j in i.keys():
        key.append(j)
        realpart.append(i[j].real)
        imagpart.append(i[j].imag)
#        print(j)
#        print(i[j].real)
#        print(i[j].imag)
    kk=kk+1
    dic = { "ky" : key , "rpt": realpart , "ipt": imagpart}
    df = pd.DataFrame(dic)
    df.to_csv(path_or_buf=filename)
    print(df)