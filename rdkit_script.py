from rdkit import Chem

# suppl = Chem.SDMolSupplier('mk/mk_vina_gmo.sdf')
# suppl = Chem.SDMolSupplier('mk/mk_terb_vina.sdf')
# suppl = Chem.SDMolSupplier("mk/mk_pg_vina.sdf")
suppl = Chem.SDMolSupplier("mk/mk_sa_vina.sdf")


writer = Chem.SDWriter('mk/mk_sa_vina_f68.sdf')

count  = 0
for mol in suppl:
  count += 1
  if count == 68:
        # Write the 18th molecule to the output SDF file
        mol = Chem.AddHs(mol, addCoords=True)

        atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
        num_atoms = mol.GetNumAtoms()
        print("Atoms:", atoms)
        print("Number of atoms:", num_atoms)

        writer.write(mol)
        break
