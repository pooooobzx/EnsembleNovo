

f = open("denovo_clip_new.txt")
casa = list(f.readlines())
f = open("denovo_prime_new.txt")
prime = list (f.readlines())
f = open("denovo_reverse_clip_new.txt")
reverse = list(f.readlines())


all_correct = 0
only_prime = 0
only_reverse= 0
only_casa = 0
casa_prime = 0
casa_reverse = 0
prime_revese=0
for idx, _ in enumerate(reverse):
    r =  not "incorrect" in reverse[idx]
    p = not "incorrect" in prime[idx]
    c = not "incorrect" in casa[idx]
    if r and p and c:
        all_correct += 1
    if r and p and not c:
        prime_revese += 1
    if r and c and not p :
        casa_reverse +=1 
    if p and c and not r:
        casa_prime += 1 
    if p and not c and not r:
        only_prime += 1 
    if c and not p and not r:
        only_casa +=1 
    if r and not p and not c :
        only_reverse += 1
print("all_correct,", all_correct)
print("only_prime,", only_prime)
print("only_casa,", only_casa)
print("only_reverse,", only_reverse)
print("casa_reverse,", casa_reverse)
print("casa_prime", casa_prime)
print("prime_reverse", prime_revese)
    