import sys
rot=lambda s:''.join(chr((ord(c)-65+13)%26+65) if c.isalpha() else c for c in s)
PA=sys.argv[1]; PB=rot(PA); n=int(sys.argv[2]); sp=sys.argv[3]
eng="$LUB^B?#}>C$=?>$CLUP%}YH"
A=PA+sp+eng.replace("#",PA[:n]); B=PB+sp+eng.replace("#",PB[:n])
open("cross.txt","w").write(f"R A {A}\nR B {B}\n")
open("pp.txt","w").write(f"{A} {B}\n{B} {A}\n{A} {A}\n{B} {B}\n")
print(A);print(B)
