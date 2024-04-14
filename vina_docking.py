from vina import Vina

def docking(protein, ligands):
    for ligand in ligands:
        v = Vina(sf_name="vina")
        v.set_receptor(protein)

        name = ligand.split(".pdbqt")[0]

        v.set_ligand_from_file(ligand)
        binding_site_position = [0, 0, 0]  # The ligand should be placed in this position as well
        box_size = [30, 30, 30]
        v.compute_vina_maps(center=binding_site_position, box_size=box_size)

        print(v.info())
        # Score the current pose
        energy = v.score()
        print('Score before minimization: %.3f (kcal/mol)' % energy[0])

        # Minimized locally the current pose
        energy_minimized = v.optimize()
        print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
        v.write_pose(name + "_minimized.pdbqt", overwrite=True)

        # Dock the ligand
        v.dock(exhaustiveness=32, n_poses=70)
        v.write_poses(name+'_vina_out.pdbqt', n_poses=70, overwrite=True)

        print(f"------------ Finished docking: {ligand} ------------")
    print(v.cite())

docking("muc1_h_stripped_model1.pdbqt", ["gmo_chemspider.pdbqt", "ter_protonated_pubchem.pdbqt", "sa_3d_pubchem.pdbqt", "pg_pubchem.pdbqt"])