import re

seq="CGCGATACGTTACATACATGATAGACCGCGCGCGATCATATCGCGATTATC"


print (len(re.findall(f'(?={"CGCG"})', seq)))