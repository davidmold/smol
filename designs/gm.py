import subprocess, sys
rot=lambda s:''.join(chr((ord(c)-65+13)%26+65) if c.isalpha() else c for c in s)
PM="WGEWLHHRJX"; PG=rot(PM); LBL=sys.argv[1] if len(sys.argv)>1 else "LUBO"; pad=sys.argv[2] if len(sys.argv)>2 else ""
M=PM+"QQ"+"^B?"+PM[:4]+"}"+"$=?>%"+"^F$"+LBL+">B^C$=?>%}"
enc=subprocess.check_output(["../target/release/smol-codon","encode",M]).decode().strip()
G=PG+"QQ"+rot(LBL)+pad+enc
open("gm.txt","w").write(f"R M {M}\nR G {G}\n")
print("M",len(M),M); print("G",len(G),G)
print("decode", subprocess.check_output(["../target/release/smol-codon","decode",enc]).decode().strip())
