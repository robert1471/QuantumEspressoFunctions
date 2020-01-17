from itertools import permutations

elements = ["Mn", "Mn", "Fe", "Fe", "Fe", "Fe"]
Permutations = list(permutations(elements))

k = 0
while k != 15:
    k = len(Permutations)
    for i, p in enumerate(Permutations):
        for j, q in enumerate(Permutations):
            if i != j:
                if Permutations[i] == Permutations[j]:
                    Permutations.remove(Permutations[i])




print(len(Permutations))
for p in Permutations:
    print(p)
