from horton import *
import numpy as np
np.set_printoptions(linewidth=200)
filename = 'h2o_sto3g.fchk'
mol = IOData.from_file(filename)
# Note there can also be exp_beta
exp_mo = mol.exp_alpha
coeffs = exp_mo.coeffs
# Integals in atomic basis
olp_ab = mol.obasis.compute_overlap(mol.lf)
kin_ab = mol.obasis.compute_kinetic(mol.lf)
na_ab = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
er_ab = mol.obasis.compute_electron_repulsion(mol.lf)

'''
# manually construct fock matrix
# Integals in mo basis
#olp_mo = olp_ab.new()
#olp_mo.assign_two_index_transform(olp_ab, exp_mo)
kin_mo = kin_ab.new()
kin_mo.assign_two_index_transform(kin_ab, exp_mo)
na_mo = na_ab.new()
na_mo.assign_two_index_transform(na_ab, exp_mo)
er_mo = er_ab.new()
er_mo.assign_four_index_transform(er_ab, exp_mo)

core_ab = 0.5*kin_ab._array + na_ab._array
core_mo = coeffs.T.dot(core_ab).dot(coeffs)
coulomb_ab = np.einsum('ijij->ij', er_ab._array)
exchange_ab = np.einsum('ijji->ij', er_ab._array)
fock_ab = core_ab + 2*coulomb_ab - exchange_ab


core_mo = kin_mo._array + na_mo._array
coulomb_mo = np.einsum('ijij->ij', er_mo._array)
exchange_mo = np.einsum('ijji->ij', er_mo._array)
fock_mo = 2*np.diag(np.diag(core_mo)) + coulomb_mo - exchange_mo
print fock_mo
print exp_mo.energies, 'help'
print 2*np.diag(core_mo) + np.sum(2*coulomb_mo - exchange_mo, axis=1), 'im a robot'
print exp_mo.energies - (2*np.diag(core_mo) + np.sum(2*coulomb_mo - exchange_mo, axis=1))
'''

# construct effective hamiltonian from density matrix, then construct fock
lf = DenseLinalgFactory(mol.obasis.nbasis)
terms = [
    RTwoIndexTerm(kin_ab, 'kin'),
    RDirectTerm(er_ab, 'hartree'),
    RExchangeTerm(er_ab, 'x_hf'),
    RTwoIndexTerm(na_ab, 'ne'),
]
ham = REffHam(terms)
fock_alpha = lf.create_two_index()
'''
# get 1dms by construction
#dm_alpha = lf.create_two_index()
#dm_alpha.assign(2*mol.exp_alpha.coeffs.dot(mol.exp_alpha.coeffs.T))
#ham.reset(dm_alpha)
# or from somewhere else
# ham.cache['dm_alpha'] = mol.dm_alpha
'''
fock_beta = lf.create_two_index()
'''
# repeat for beta orbitals
#dm_beta = lf.create_two_index()
#dm_beta.assign(mol.exp_beta.coeffs.dot(mol.exp_beta.coeffs.T))
#ham.compute_fock(fock_alpha)
#temp = fock_alpha._array
#print ham.compute_energy()
'''
scf_solver = PlainSCFSolver()
occ_model = AufbauOccModel(5)
temp_exp_alpha = mol.exp_alpha.new()
temp_exp_alpha._coeffs = np.identity(mol.obasis.nbasis)
scf_solver(ham, mol.lf, olp_ab, occ_model, exp_mo)
print ham.compute_energy()
dm_one = ham.cache['dm_alpha']._array

#print'xxx'
#print dm_one - coeffs.dot(coeffs.T)
#print dm_one - 2*coeffs.dot(coeffs.T)

#ham.compute_fock(fock_alpha)
#print mol.exp_alpha.coeffs
