{
    //const vectorField& cells = mesh.C();
    scalarField& rho_at_cells  = rho.primitiveFieldRef();
    scalarField& visc_at_cells = visc.primitiveFieldRef();
    scalarField& D_at_cells    = diffcoeff.primitiveFieldRef();
    scalarField& conc_at_cells = C.primitiveFieldRef();

    forAll(rho_at_cells,i)
    {
        rho_at_cells[i]  = rhoCoeffs.value(conc_at_cells[i]);
        visc_at_cells[i] = muCoeffs.value(conc_at_cells[i]);
        D_at_cells[i]    = DCoeffs.value(conc_at_cells[i]);
    }

    volScalarField::Boundary& rhoBf    = rho.boundaryFieldRef();
    volScalarField::Boundary& viscBf   = visc.boundaryFieldRef();
    volScalarField::Boundary& dcoeffBf = diffcoeff.boundaryFieldRef();
    volScalarField::Boundary& concBf   = C.boundaryFieldRef();

    forAll(rho.boundaryField(), patchi)
    {
        fvPatchScalarField& rhop    = rhoBf[patchi];
        fvPatchScalarField& viscp   = viscBf[patchi];
        fvPatchScalarField& dcoeffp = dcoeffBf[patchi];
        fvPatchScalarField& concp   = concBf[patchi];

        forAll(rhop, facei)
        {
            rhop[facei]    = rhoCoeffs.value(concp[facei]);
            viscp[facei]   = muCoeffs.value(concp[facei]);
            dcoeffp[facei] = DCoeffs.value(concp[facei]);
        }
    }
}
