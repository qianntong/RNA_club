from nupack import *

# https://docs.nupack.org/analysis/

model1 = Model(material='rna', celsius=37)

A = Strand('AGUCUAGGAUUCGGCGUGGGUUAA', name='A') # name is required for strands
B = Strand('UUAACCCACGCCGAAUCCUAGACUCAAAGUAGUCUAGGAUUCGGCGUG', name='B')
C = Strand('AGUCUAGGAUUCGGCGUGGGUUAACACGCCGAAUCCUAGACUACUUUG', name='C')

A.nt() # --> 24

c1 = Complex([A]) # name is optional for complexes
c2 = Complex([A, B, B, C], name='ABBC')
c3 = Complex([A, A], name='AA')

# destabilize c4 by 1 kcal/mol
c4 = Complex([A, B, C], name='ABC', bonus=+1.0)

# stabilize c5 by 10 kcal/mol
c5 = Complex([A, B], name='AB', bonus=-10.0)

t1 = Tube(strands={A: 1e-6, B: 1e-8}, name='t1') # complexes defaults to [A, B]

t2 = Tube(strands={A: 1e-6, B: 1e-8, C: 1e-12},
    complexes=SetSpec(max_size=3, include=[c2,[B, B, B, B]], exclude=[c1]),
    name='t2')

print(t1.complexes) # --> {<Complex A>, <Complex B>}
print(t2.complexes) # --> {<Complex (C+C+B)>, <Complex (B)>,
    # <Complex (A+C+B)>, <Complex (C+C+C)>, <Complex (C)>, <Complex (A+A+B)>,
    # <Complex (A+C)>, <Complex (B+B+B+B)>, <Complex (A+A)>, <Complex (A+B+B)>,
    # <Complex (B+B)>, <Complex (A+B)>, <Complex (B+B+B)>, <Complex (A+B+C)>,
    # <Complex (A+C+C)>, <Complex (A+A+A)>, <Complex (C+C)>, <Complex (A+A+C)>,
    # <Complex ABBC>, <Complex (C+B+B)>, <Complex (C+B)>}

# specify strands
a = Strand('CUGAUCGAU', name='a')
b = Strand('GAUCGUAGUC', name='b')

# specify tubes
t1 = Tube(strands={a: 1e-8, b: 1e-9}, complexes=SetSpec(max_size=3), name='t1')
t2 = Tube(strands={a: 1e-10, b: 1e-9}, complexes=SetSpec(max_size=2), name='t2')

# analyze tubes
model1 = Model()
tube_results = tube_analysis(tubes=[t1, t2], model=model1)

print(tube_results)
