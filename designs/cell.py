import subprocess
rot=lambda s:''.join(chr((ord(c)-65+13)%26+65) if c.isalpha() else c for c in s)
PM="WGEWLHHRJX"; PG=rot(PM)
# cell cycle: gate, then forever: copy G (blind) -> release, toggle translate, rewind R, translate G -> release M, toggle back, rewind, repeat
M=PM+"QQ"+"^B?"+PM[:4]+"}"+"^F>B"+"$=?>%"+"^C^F>B^F>"
enc=subprocess.check_output(["../target/release/smol-codon","encode",M]).decode().strip()
G=PG+"QQ"+enc
open("cell.txt","w").write(f"R M {M}\nR G {G}\n")
print("M",len(M),M); print("G",len(G),G)
