import moosh as ms

# Exemple :

# Structure avec 3 couches de même épaisseur

# On défini chaque couche de la structure
Layer1 = ms.Layer(1, 1, 100)  # Epsilon = 1, Mu = 1, Height = 100 nm
Layer2 = ms.Layer(2, 1, 100)  # Epsilon = 2, Mu = 1, Height = 100 nm
Layer3 = ms.Layer(4, 1, 100)  # Epsilon = 4, Mu = 1, Height = 100 nm

# On assemble les couches dans la structure
structure = ms.Structure(
    Layer1,
    Layer2,
    Layer3
)
# On obtient 3 arrays : [Eps1, Eps2, Eps3], [Mu1, Mu2, Mu3]
# et [Height1, Height2, Height3]

# On peut utiliser des structures pré-définies : Mirroir de Bragg
braggStructure = ms.StructureBragg(10)

# Calcul des coefficients de réflexion de la structure
Coef = ms.Reflection(braggStructure, 1)
Coef.angular(600)  # Réflexion en fonction de l'angle
Coef.spectrum(75)  # Réflexion en fonction de la longueur d'onde
