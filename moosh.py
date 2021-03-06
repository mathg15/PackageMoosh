import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def Layer(Epsilon, Mu, Height):
    """
    :param Epsilon: Permittivity of the layer
    :param Mu: Permeability of the layer
    :param Height: Height of the layer
    :return: Array like [Epsilon, Mu, Height]
    """
    return np.array([[Epsilon, Mu, Height]])


def Structure(*args):
    """
    Create a multi-layers structure

    :param args: layers of the structure
    :return: array for each parameter of the layers

    Example :

    Layer1 = [1,1,100]; Layer2 = [2,1,100]

    Return : Eps = [1, 2]; Mu = [1, 1}; Height = [100, 100]
    """
    structure = args[0]
    for i in range(len(args) - 1):
        structure = np.append(structure, args[i + 1], axis=0)
    epsStruc = structure[:, 0]
    muStruc = structure[:, 1]
    heightStruc = structure[:, 2]
    return epsStruc, muStruc, heightStruc


class RandomStructure:

    def __init__(self, nbLayer, Eps):
        self.Eps = Eps
        self.nbLayer = nbLayer

    def random_layer(self):
        layers = []
        i = 0
        for i in range(self.nbLayer):
            randEps = np.random.randint(0, len(self.Eps))
            randHeight = np.random.randint(10, 2000)
            layer = np.array([self.Eps[randEps], 1, randHeight])
            layers.append(layer)
            i += 1
        return np.asarray(layers)

    def random_structure(self):
        epsStruc = self.random_layer()[:, 0]
        muStruc = self.random_layer()[:, 1]
        heightStruc = self.random_layer()[:, 2]
        return epsStruc, muStruc, heightStruc


def StructureBragg(periods):
    Eps = np.tile([1, 2], (1, periods))
    Eps = np.insert(Eps, 0, 1)
    Eps = np.append(Eps, 1)

    Mu = np.ones((2 * periods + 2))

    Height = np.tile([600 / np.sqrt(2), 600 / (4 * 2)], (1, periods))
    Height = np.insert(Height, 0, 1600)
    Height = np.append(Height, 100)

    return Eps, Mu, Height


class Reflection:

    def __init__(self, structure, polarisation):
        """
        :param structure: Parameters of the multilayer Epsilon = []; Mu = []; Height = []
        :param polarisation: TE polarisation = 0 and TM polarisation = 1
        """
        self.Epsilon, self.Mu, self.Height = structure
        self.polarisation = polarisation

    def cascade(sekf, A, B):
        t = 1 / (1 - B[0, 0] * A[1, 1])
        S = np.array([[A[0, 0] + A[0, 1] * B[0, 0] * A[1, 0] * t, A[0, 1] * B[0, 1] * t],
                      [B[1, 0] * A[1, 0] * t, B[1, 1] + A[1, 1] * B[0, 1] * B[1, 0] * t]], dtype=complex)
        return S

    def coefficient(self, angle, _lambda):

        # On consid??re que la premi??re valeur de la hauteur est 0
        self.Height[0] = 0

        # On definie k0 ?? partir de lambda
        k0 = (2 * np.pi) / _lambda

        # g est la longeur de Type
        g = len(self.Epsilon)

        if self.polarisation == 0:
            f = self.Epsilon
        elif self.polarisation == 1:
            f = self.Mu

        # D??finition de alpha et gamma en fonction de TypeEps, TypeMu, k0 et Theta
        alpha = np.sqrt(self.Epsilon[0] * self.Mu[0]) * k0 * np.sin(angle * (np.pi / 180))
        gamma = np.sqrt(self.Epsilon * self.Mu * (k0 ** 2) - np.ones(g) * (alpha ** 2))

        # On fait en sorte d'avoir un r??sultat positif en foction au cas o?? l'index serait n??gatif
        if np.real(self.Epsilon[0]) < 0 and np.real(self.Mu[0]):
            gamma[0] = -gamma[0]

        # On modifie la d??termination de la racine carr??e pour obtenir un stabilit?? parfaite
        if g > 2:
            gamma[1: g - 2] = gamma[1:g - 2] * (1 - 2 * (np.imag(gamma[1:g - 2]) < 0))
        # Condition de l'onde sortante pour le dernier milieu
        if np.real(self.Epsilon[g - 1]) < 0 and np.real(self.Mu[g - 1]) < 0 and np.real(
                np.sqrt(self.Epsilon[g - 1] * self.Mu * (k0 ** 2) - (alpha ** 2))):
            gamma[g - 1] = -np.sqrt(self.Epsilon[g - 1] * self.Mu[g - 1] * (k0 ** 2) - (alpha ** 2))
        else:
            gamma[g - 1] = np.sqrt(self.Epsilon[g - 1] * self.Mu[g - 1] * (k0 ** 2) - (alpha ** 2))

        # Definition de la matrice Scattering
        T = np.zeros((2 * g, 2, 2), dtype=complex)
        T[0] = [[0, 1], [1, 0]]

        # Cacul des matrices S
        for k in range(g - 1):
            # Matrice de diffusion des couches
            t = np.exp(1j * gamma[k] * self.Height[k])
            T[2 * k + 1] = np.array([[0, t], [t, 0]])

            # Matrice de diffusion d'interface
            b1 = gamma[k] / (f[k])
            b2 = gamma[k + 1] / (f[k + 1])
            T[2 * k + 2] = [[(b1 - b2) / (b1 + b2), (2 * b2) / (b1 + b2)],
                            [(2 * b1) / (b1 + b2), (b2 - b1) / (b1 + b2)]]

        # Matrice de diffusion pour la derni??re couche
        t = np.exp(1j * gamma[g - 1] * self.Height[g - 1])
        T[2 * g - 1] = [[0, t], [t, 0]]

        A = np.zeros((2 * g - 1, 2, 2), dtype=complex)
        A[0] = T[0]

        for j in range(len(T) - 2):
            A[j + 1] = self.cascade(A[j], T[j + 1])

        # Coefficient de reflexion de l'ensemble de la structure
        r = A[len(A) - 1][0, 0]

        # Coefficient de transmission de l'ensemble de la structure
        # tr = A[len(A) - 1][1, 0]

        # Coefficient de r??flexion de l'??nergie
        Re = np.abs(r) ** 2

        # Coefficient de transmission de l'??nergie
        # Tr = np.abs(tr)  # ** 2 * gamma[g - 1] * f[0] / (gamma[0] * f[g-1])

        return r, Re

    def angular(self, _lambda_):
        """
        :param _lambda_ : wave length
        :return: Reflection vs incident angle (0?? to 90??)
        """

        rangeAngle = np.linspace(0, 89, 200)

        # Cr??ation des matrices
        re = np.ones((200, 1), dtype=complex)
        Re = np.ones((200, 1), dtype=complex)
        # c = np.ones((200, 1), dtype=complex)
        # d = np.ones((200, 1), dtype=complex)

        for i in range(200):
            tht = rangeAngle[i]
            re[i], Re[i] = self.coefficient(tht, _lambda_)

        plt.title(f"Reflection for lambda = {_lambda_} nm")
        plt.plot(rangeAngle, abs(Re))
        plt.ylabel("Reflection")
        plt.xlabel("Angle (degrees)")
        plt.grid(True)
        plt.show()

    def spectrum(self, angle):
        """
        :param angle: incident angle
        :return: Reflection vs wave length (400nm to 800nm)
        """
        rangeLambda = np.linspace(400, 800, 1000)

        # Cr??ation des matrices
        re = np.ones((1000, 1), dtype=complex)
        Re = np.ones((1000, 1), dtype=complex)
        # c = np.ones((1000, 1), dtype=complex)
        # = np.ones((1000, 1), dtype=complex)

        for i in range(1000):
            lambda_ = rangeLambda[i]
            re[i], Re[i] = self.coefficient(angle * (np.pi / 180), lambda_)

        plt.title(f"Reflection spectrum for an incidence of {angle}??")
        plt.plot(rangeLambda, abs(Re))
        plt.ylabel("Reflection")
        plt.xlabel("Wavelength")
        plt.grid(True)
        plt.show()


class Beam:

    def __init__(self, structure, pola, angle, lambda_, pos):

        self.structure = structure
        self.pola = pola
        self.angle = angle
        self.lambda_ = lambda_
        self.pos = pos

    def cascade(self, A, B):
        # On combine deux matrices de diffusion A et B en une seule matrice de diffusion (2*2) S
        t = 1 / (1 - B[0, 0] * A[1, 1])
        S = np.array([[A[0, 0] + A[0, 1] * B[0, 0] * A[1, 0] * t, A[0, 1] * B[0, 1] * t],
                      [B[1, 0] * A[1, 0] * t, B[1, 1] + A[1, 1] * B[0, 1] * B[1, 0] * t]], dtype=complex)
        return S

    def beam(self):

        TypeEps, TypeMu, hauteur = self.structure

        _theta = self.angle * (np.pi / 180)

        # Spatial window size
        d = 70 * self.lambda_

        # Incident beam width
        w = 10 * self.lambda_

        # Number of pixels horizontally
        nx = np.floor(d / 10)

        # Number of pixels verticaly
        ny = np.floor(hauteur / 10)
        # print(ny)
        # Number of modes
        nmod = np.floor(0.83660 * (d / w))

        hauteur = hauteur / d
        l = self.lambda_ / d
        w = w / d

        if self.pola == 0:
            f = TypeEps
        elif self.pola == 1:
            f = TypeMu

        k0 = (2 * np.pi) / l

        En = np.zeros((int(np.sum(ny)), int(nx)), dtype=complex)

        # g est la longeur de Type
        g = len(TypeEps)

        # Amplitude of the different modes
        X = np.exp(-(w ** 2) * (np.pi ** 2) * np.arange(-nmod, nmod + 1) ** 2) * np.exp(
            2 * 1j * np.pi * np.arange(-nmod, nmod + 1) * self.pos)

        # Scattering matrix
        T = np.zeros((2 * g, 2, 2), dtype=complex)
        T[0] = [[0, 1], [1, 0]]

        for nm in range(int(2 * nmod)):
            alpha = np.sqrt(TypeEps[0] * TypeMu[0]) * k0 * np.sin(_theta) + 2 * np.pi * (nm - nmod - 1)
            # print("alpha :",alpha)
            gamma = np.sqrt(TypeEps * TypeMu * (k0 ** 2) - np.ones(g) * (alpha ** 2))

            if np.real(TypeEps[0]) < 0 and np.real(TypeMu[0]):
                gamma[0] = -gamma[0]

            # On modifie la d??termination de la racine carr??e pour obtenir un stabilit?? parfaite
            if g > 2:
                gamma[1: g - 2] = gamma[1:g - 2] * (1 - 2 * (np.imag(gamma[1:g - 2]) < 0))
            # Condition de l'onde sortante pour le dernier milieu
            if np.real(TypeEps[g - 1]) < 0 and np.real(TypeMu[g - 1]) < 0 and np.real(
                    np.sqrt(TypeEps[g - 1] * TypeMu * (k0 ** 2) - (alpha ** 2))):
                gamma[g - 1] = -np.sqrt(TypeEps[g - 1] * TypeMu[g - 1] - (alpha ** 2))
            else:
                gamma[g - 1] = np.sqrt(TypeEps[g - 1] * TypeMu[g - 1] * (k0 ** 2) - (alpha ** 2))

            for k in range(g - 1):
                # Matrice de diffusion des couches
                t = np.exp(1j * gamma[k] * hauteur[k])
                T[2 * k + 1] = np.array([[0, t], [t, 0]])

                # Matrice de diffusion d'interface
                b1 = gamma[k] / (f[k])
                b2 = gamma[k + 1] / (f[k + 1])
                T[2 * k + 2] = [[(b1 - b2) / (b1 + b2), (2 * b2) / (b1 + b2)],
                                [(2 * b1) / (b1 + b2), (b2 - b1) / (b1 + b2)]]

            # Matrice de diffusion pour la derni??re couche
            t = np.exp(1j * gamma[g - 1] * hauteur[g - 1])
            T[2 * g - 1] = [[0, t], [t, 0]]

            A = np.zeros((2 * g - 1, 2, 2), dtype=complex)
            A[0] = T[0]

            H = np.zeros((2 * g - 1, 2, 2), dtype=complex)
            H[0] = T[2 * g - 1]

            for j in range(len(T) - 2):
                A[j + 1] = self.cascade(A[j], T[j + 1])
                H[j + 1] = self.cascade(T[len(T) - 2 - j], H[j])

            I = np.zeros((len(T), 2, 2), dtype=complex)
            for j in range(len(T) - 1):
                I[j] = np.array([[A[j][1, 0], A[j][1, 1] * H[len(T) - j - 2][0, 1]],
                                 [A[j][1, 0] * H[len(T) - j - 2][0, 0], H[len(T) - j - 2][0, 1]]] / (
                                        1 - A[j][1, 1] * H[len(T) - j - 2][0, 0]))
            h = 0
            t = 0

            E = np.zeros((int(np.sum(ny)), 1), dtype=complex)

            for k in range(g):
                for m in range(int(ny[k])):
                    h = h + hauteur[k] / ny[k]
                    E[t, 0] = I[2 * k][0, 0] * np.exp(1j * gamma[k] * h) + I[2 * k + 1][1, 0] * np.exp(
                        1j * gamma[k] * (hauteur[k] - h))
                    t += 1
                h = 0

            E = E * np.exp(1j * alpha * np.arange(0, nx) / nx)
            En = En + X[int(nm)] * E

        V = np.abs(En)
        V = np.flip(V)
        V = V / V.max()

        norm = mcolors.Normalize(vmax=V.max(), vmin=V.min())

        plt.figure(1)
        plt.pcolormesh(V / V.max(), norm=norm, cmap='jet')
        plt.colorbar()
        plt.title(f"Light beam for lambda = {self.lambda_} \n with an incidence angle of {self.angle} degrees")
