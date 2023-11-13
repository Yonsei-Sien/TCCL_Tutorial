from pyscf import gto,dft,scf

'''
< Unit List >

Length Unit: Angstrom
Energy Unit: kcal/mol

'''

class DFT:
        def __init__(self, atom='', spin=0, charge=0, xc='r2scan', basis='def2qzvppd'):
                self.coord=atom
                self.spin=spin
                self.charge=charge
                self.xc=xc
                self.basis=basis
                self.hf_E=0
                self.ks_E=0
        
        def data_check(self):
                if self.coord=='' or self.basis=='':
                        return False
                else:
                        return True

        def ks(self, verbose=4, max_cycle=500, level_shift=0):
                if self.data_check():
                        mol=gto.M(atom=self.coord
                                ,spin=self.spin
                                ,charge=self.charge
                                ,basis=self.basis
                                ,verbose=verbose)
                        
                        mol_sc=dft.UKS(mol)
                        mol_sc.xc=self.xc
                        mol_sc.level_shift=level_shift
                        mol_sc.conv_check=False
                        mol_sc.grids.level=4
                        mol_sc.conv_tol=1e-8
                        mol_sc.diis_space=40
                        mol_sc.max_cycle=max_cycle

                        Esc=mol_sc.kernel() * 627.509

                        self.ks_E=Esc
                        return Esc
                else:
                        print('Data Available Check False: Please Check Your Initial Input')



        def hf(self,level_shift=0,verbose=4):
                if self.data_check():
                        mol=gto.M(atom=self.coord
                                ,spin=self.spin
                                ,charge=self.charge
                                ,basis=self.basis
                                ,verbose=verbose)
                        
                        mol_sc=dft.UKS(mol)
                        mol_sc.xc=self.xc
                        mol_sc.level_shift=level_shift
                        mol_sc.conv_check=False
                        mol_sc.grids.level=4
                        mol_sc.conv_tol=1e-8
                        mol_sc.diis_space=40
                        
                        mol_hf=scf.UHF(mol)
                        mol_hf.diis=scf.ADIIS()
                        mol_hf.diis_space=40
                        mol_sc.conv_tol=1e-8
                        mol_hf.kernel()
                        Dhf=mol_hf.make_rdm1()
                        Ehf=mol_sc.energy_tot(Dhf)

                        Ehf *= 627.509

                        self.hf_E=Ehf
                        return Ehf
                else:
                        print('Data Available Check False: Please Check Your Initial Input')
        
