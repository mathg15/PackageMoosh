import moosh as ms

# Définition des couches
Layer1 = ms.Layer(1, 1, 100).layer()
Layer2 = ms.Layer(2, 1, 100).layer()

# Définition de la structure de la multicouche
structure = ms.Structure(
    Layer1,
    Layer2
).structure()

# On utilise ensuite la fonction Angular pour construire le spectre
ms.Angular(structure, 0, 600).angular()
