import moosh as ms

Layer1 = ms.Layer(1, 1, 100).layer()
Layer2 = ms.Layer(2, 1, 100).layer()

structure = ms.Structure(
    Layer1,
    Layer2
).structure()

ms.Angular(structure, 0, 600).angular()
