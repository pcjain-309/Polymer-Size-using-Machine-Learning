f=open("def.chain","r+")
l=f.readlines()
m = int(l[8].split(None,1)[0])
print(m)
l[8]=l[8].replace(str(m),str(m+1))
f.close()
f1=open("def.chain","w")
f1.writelines(l)
f1.close()