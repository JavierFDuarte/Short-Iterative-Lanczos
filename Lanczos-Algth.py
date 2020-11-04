#Partimos de un Estado inicial que es autoestado del oscilador armónico y creamos base de Lanczos por aplicación de operadores de aniquilación
#y creación básicamente. Entonces los vectores de Base Lanczos VAN A SER a lo sumo combinación lineal de autoest del oscilador armónico

#En el código los Vectores Base de Lanczos (las combinac lineales) están representadas por Diccionarios siendo:
# key -----> El orden de 1 de los autoestados términos de la combinación lineal
# valor ----> El coeficiente que acompaña al autoestado
#
#   |v> = SUM[ c_k*|k> ]  -------> Dic = { k : c_k }


##########################################################################################################################################
##########################################################################################################################################

# Import Necesary Packages
import math as mth      #sqrt function
import cmath as cmth    #Complex algebra
import pandas as pd     #read and write data
import numpy as np      #read and write data

##########################################################################################################################################
##########################################################################################################################################

#Functions
#Raising operator
def raising(dc):
    dic2 = { }
    ls2 = list(dc.keys())
    ls2.sort()
    for i in ls2:
        dic2[i+1] = mth.sqrt(i+1)*dc[i]
    return dic2

#Lowering operator
def lowering(dc):
    dic3 = { }
    ls3 = list(dc.keys())
    ls3.sort()
    for i in ls3:
        if i == 0:
            pass
        else:
            dic3[i-1] = mth.sqrt(i)*dc[i]
    return dic3

#Number operator
def num(dc):
    dic4 = { }
    ls4 = list(dc.keys())
    ls4.sort()
    for i in ls4:
        dic4[i] = i*dc[i]
    return dic4

#Hamiltonian operator
def hamiltonian(dc):
    dc_a = raising(raising(dc))
    dc_b = lowering(lowering(dc))
    dc_c = num(dc)
    ham = { }
    ls5 = list(dc.keys())
    ls5.sort()
    for i in ls5:
        ham[i] = -dc[i]  -2.*dc_c[i]
    for i in dc_a.keys():
        ham[i] = dc_a[i] + ham.get(i,0.)
    for i in dc_b.keys():
        ham[i] = dc_b[i] + ham.get(i,0)
    l = list(ham.keys())
    l.sort()
    ham2 = { }
    for i in l:
        ham2[i] = -1*ham[i]
    return ham2

#Kronecker Delta function
def delta(u,v):
    if u == v:
        d=1.0
    else:
        d=0.0
    return d

#Inner Product functional
def inner_prod(dc1,dc2):
    ls6 = list(dc1.keys())
    ls6.sort()
    ls7 = list(dc2.keys())
    ls7.sort()
    suma = 0.0
    for i in ls6:
        for j in ls7:
            v = dc1[i]
            w = dc2[j]
            suma = suma + v.conjugate()*w*delta(i,j)
    return suma

#Normalization operation
def normaliz(dc):
    ls8 = list(dc.keys())
    ls8.sort()
    norm = inner_prod(dc,dc).real
    norm=mth.sqrt(norm)
    normalizado = { }
    for i in ls8:
        normalizado[i] = 1./norm * dc[i]
    return normalizado , norm

##########################################################################################################################################
##########################################################################################################################################

#Lanczos Algorithm

n=4

#Initial State read from file "initialstate.txt"
#Coefficients and recursive Krylov's Vectors stored on lists
df1 = pd.read_table("initialstate.txt",sep=r'\s{2,}')
print(df1   .shape)
ns = list(df1["n"])
re_cns=list(df1["r_cn"])
im_cns=list(df1["i_cn"])
print("n=",ns)
print("re_cn=",re_cns)
print("im_cn=",im_cns)
print(df1)


alpha_array = [ ]
array_beta = [ ]
array_vectors = [ ]

#Initial State
kk = 0
v0 = { }
for i in ns:
    v0[i] = complex(re_cns[kk],im_cns[kk])
    kk=kk+1

v0,normaignorar = normaliz(v0)
#print("v0 =", v0)
array_vectors.append(v0)


alpha = inner_prod(v0 , hamiltonian(v0))
alpha_array.append(alpha)
#print(alpha_array[0])



vk = array_vectors[0]
hv=hamiltonian(vk)
v1 = hv

for i in v0.keys():
    v1[i] = hv.get(i,0.0) - alpha_array[0]*v0[i]
    if abs(v1[i]) <= 0.00000001:
        del v1[i]
#print(v1)

v1,beta = normaliz(v1)
array_beta.append(beta)
array_vectors.append(v1)



###########################################
###########################################


for k in range(1,n-1):
    alpha = inner_prod(array_vectors[k] , hamiltonian(array_vectors[k]))
    alpha_array.append(alpha)



    hv=hamiltonian(array_vectors[k])
    vkplus = hv
    vk = array_vectors[k]
    for i in vk.keys():
        vkplus[i] = hv.get(i,0.0) - alpha_array[k]*vk[i]
        if abs(vkplus[i]) <= 0.00000001:
            del vkplus[i]

    vkles = array_vectors[k-1]
    for i in vkles.keys():
        vkplus[i] = vkplus.get(i,0.0) - array_beta[k-1]*vkles[i]
        if abs(vkplus[i]) <= 0.00000001:
            del vkplus[i]




    vkplus,beta = normaliz(vkplus)
    #print("BETA=",beta)
    array_vectors.append(vkplus)
    array_beta.append(beta)

alpha = inner_prod(vkplus , hamiltonian(vkplus))
alpha_array.append(alpha)



###########################################
###########################################
print("                              ")
print("Alpha_Array")
print(alpha_array)

print("                              ")
print("Beta_Array")
print(array_beta)

print("                              ")
print("Vector_Array")
print(array_vectors)
##########################################
##########################################

realalpha = [ ]
imagalpha = [ ]


for k in alpha_array:
    realalpha.append(k.real)
    imagalpha.append(k.imag)

print(realalpha)
print(imagalpha)

#Save Coefficients on file
coefalpha = {"realalpha":realalpha,"imagalpha":imagalpha}
dfcoefalpha = pd.DataFrame(coefalpha)
print(dfcoefalpha)
filename1= open("coeficientesAlpha.txt","w")
np.savetxt(filename1,dfcoefalpha.values)
filename1.close()

coefbeta = {"beta":array_beta}
dfcoefbeta = pd.DataFrame(coefbeta)
print(dfcoefbeta)
filename2= open("coeficientesBeta.txt","w")
np.savetxt(filename2,dfcoefbeta.values)
filename2.close()



# TENGO QUE GUARDAR TAMBIÉN EL VECTOR ARRAY DE ALGUNA FORMA PARA DEJAR LISTO PARA SIGUIENTE ITERACION!!!!!

kk=0
for i in array_vectors:
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



# Check Funcionamiento Raising and Lowering Functions
#estado =  {5:1 , 0:1 , 2:1}
#print(estado)
#estado2 = raising(raising(estado))
#print(estado2)
#estado3 = lowering(lowering(estado))
#print(estado3)
#estado4 = num(estado)
#print(estado4)

# Check Funcionamiento Hamiltonian Function
#estado = {2:1. , 3:1.}
#h = hamiltonian(estado)
#print(h)
#funciona carajo funciona!

#Check Funcionamiento Inner_Product Function
#dic = {100:1/mth.sqrt(5.) , 101:1/mth.sqrt(3.) , 149:1/mth.sqrt(3.)}
#dicc = {104:1/mth.sqrt(3.) , 100:1/mth.sqrt(3.) , 159:1/mth.sqrt(3.)}
#print(inner_prod(dic,dicc))
