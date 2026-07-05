seq="GAGCCACCGCGATA"
gc_dif=0 #initial skew set to zero
print(gc_dif,end=" ") 
for i in seq:
    if i=="C":
        gc_dif-=1
        print(gc_dif,end=" ")
    elif i=="G":
        gc_dif+=1
        print(gc_dif,end=" ")
    else: 
        print(gc_dif,end=" ")


