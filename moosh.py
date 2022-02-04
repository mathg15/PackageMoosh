import numpy as np
import matplotlib as plt


class Layer:

    def __init__(self, eps, mu, height):
        self.eps = eps
        self.mu = mu
        self.height = height

    def layer(self):
        return np.array([[self.eps, self.mu, self.height]])


class Structure:

    def __init__(self, *args):
        layer = args[0]
        for i in range(len(args) - 1):
            layer = np.append(layer, args[i + 1], axis=0)
        self.st = layer

    def structure(self):
        epsStruc = self.st[:, 0]
        muStruc = self.st[:, 1]
        heightStruc = self.st[:, 2]
        return epsStruc, muStruc, heightStruc


class Angular:

    def __init__(self, structure, polarisation, lambda_):
        self.Eps, self.Mu, self.Height = structure
        self.polarisation = polarisation
        self.lambda_ = lambda_
    
    def cascade(self, A, B):
        # On combine deux matrices de diffusion A et B en une seule matrice de diffusion (2*2) S
        t = 1 / (1 - B[0, 0] * A[1, 1])
        S = np.array([[A[0, 0] + A[0, 1] * B[0, 0] * A[1, 0] * t, A[0, 1] * B[0, 1] * t],
                      [B[1, 0] * A[1, 0] * t, B[1, 1] + A[1, 1] * B[0, 1] * B[1, 0] * t]], dtype=complex)
        return S


    def coefficient(self, thetacoef, _lambda):


        # On considère que la première valeur de la hauteur est 0
        self.Height[0] = 0

        # On definie k0 à partir de lambda
        k0 = (2 * np.pi) / _lambda

        # g est la longeur de Type
        g = len(Type)

        if self.polarisation == 0:
            f = self.Eps
        elif self.pol == 1:
            f = self.Mu

        # Définition de alpha et gamma en fonction de TypeEps, TypeMu, k0 et Theta
        alpha = np.sqrt(self.Eps[0] * self.Mu[0]) * k0 * np.sin(thetacoef * (np.pi / 180))
        gamma = np.sqrt(self.Eps * self.Mu * (k0 ** 2) - np.ones(g) * (alpha ** 2))
        # print("gamma=",gamma)

        # On fait en sorte d'avoir un résultat positif en foction au cas où l'index serait négatif
        if np.real(self.Eps[0]) < 0 and np.real(self.Mu[0]):
            gamma[0] = -gamma[0]

        # On modifie la détermination de la racine carrée pour obtenir un stabilité parfaite
        if g > 2:
            gamma[1: g - 2] = gamma[1:g - 2] * (1 - 2 * (np.imag(gamma[1:g - 2]) < 0))
        # Condition de l'onde sortante pour le dernier milieu
        if np.real(self.Eps[g - 1]) < 0 and np.real(self.Mu[g - 1]) < 0 and np.real(
                np.sqrt(self.Eps[g - 1] * self.Mu * (k0 ** 2) - (alpha ** 2))):
            gamma[g - 1] = -np.sqrt(self.Eps[g - 1] * self.Mu[g - 1] * (k0 ** 2) - (alpha ** 2))
        else:
            gamma[g - 1] = np.sqrt(self.Eps[g - 1] * self.Mu[g - 1] * (k0 ** 2) - (alpha ** 2))

        # Définition de la matrice Scattering
        T = np.zeros((2 * g, 2, 2), dtype=complex)
        T[0] = [[0, 1], [1, 0]]

        # Cacul des matrices S
        for k in range(g - 1):
            # Matrice de diffusion des couches
            t = np.exp(1j * gamma[k] * hauteur[k])
            T[2 * k + 1] = np.array([[0, t], [t, 0]])

            # Matrice de diffusion d'interface
            b1 = gamma[k] / (f[k])
            b2 = gamma[k + 1] / (f[k + 1])
            T[2 * k + 2] = [[(b1 - b2) / (b1 + b2), (2 * b2) / (b1 + b2)],
                            [(2 * b1) / (b1 + b2), (b2 - b1) / (b1 + b2)]]

        # Matrice de diffusion pour la dernière couche
        t = np.exp(1j * gamma[g - 1] * hauteur[g - 1])
        T[2 * g - 1] = [[0, t], [t, 0]]

        A = np.zeros((2 * g - 1, 2, 2), dtype=complex)
        A[0] = T[0]

        for j in range(len(T) - 2):
            A[j + 1] = self.cascade(A[j], T[j + 1])

        # Coefficient de reflexion de l'ensemble de la structure
        r = A[len(A) - 1][0, 0]

        # Coefficient de transmission de l'ensemble de la structure
        tr = A[len(A) - 1][1, 0]

        # Coefficient de réflexion de l'énergie
        Re = np.abs(r) ** 2

        # Coefficient de transmission de l'énergie
        Tr = (np.abs(tr) ** 2) * gamma[g - 1] * f[0] / (gamma[0] * f[g - 1])

        return r, tr, Re, Tr

    def angular(self, lambda__):

        
        rangeAngle = np.linspace(0, 90, 200)

        # Création des matrices
        a = np.ones((200, 1), dtype=complex)
        b = np.ones((200, 1), dtype=complex)
        c = np.ones((200, 1), dtype=complex)
        d = np.ones((200, 1), dtype=complex)

        for i in range(200):
            tht = rangeAngle[i]
            a[i], b[i], c[i], d[i] = self.coefficient(tht, lambda__)

        plt.figure(1)
        plt.subplot(211)
        plt.title(f"Reflection for lambda = {lambda__}")
        plt.plot(rangeAngle, abs(c))
        plt.plot(1, reflectMin)
        plt.ylabel("Reflection")
        plt.xlabel("Angle (degrees)")
        plt.grid(True)
        plt.subplot(212)
        plt.plot(rangeAngle, np.angle(a))
        plt.ylabel("Phase")
        plt.xlabel("Angle")
        plt.title("Phase of the reflection coefficient")
        plt.tight_layout()
        plt.show()

    def __str__(self):
        
        self.angular(self.lambda__)
